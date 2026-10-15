"""A continuum tissue study with conservative material transfers and paid activity.

This is a different law from CellularEcology: no cells, contacts, inherited
controllers, or mechanical force solve. All simulation work remains in Julia.
"""
module LivingTissueStudy
using Interior, Random, Statistics, SHA, TOML, Serialization
using Interior: KernelAbstractions, JSON3
using KernelAbstractions: @kernel, @index
import Interior: validate, stability_limit, advance!, process_ports,
    record_process_fields, owned_material_mass
export LivingTissue, tissue_world, tissue_statistics, run_living_tissue,
    tissue_checkpoint, restore_tissue!, tissue_parameters

const CHANNELS=(:precursor,:biomass,:matrix,:spent,:excitation,:inhibition,:memory)
const MODEL_SOURCE=read(@__FILE__,String)
const MODEL_SHA=bytes2hex(sha256(MODEL_SOURCE))
const SHARED_SOURCES=Dict(name=>read(joinpath(dirname(pathof(Interior)),name),String)
    for name in ("core.jl","recurrent_material.jl"))
const SHARED_SHA=Dict(name=>bytes2hex(sha256(source)) for (name,source) in SHARED_SOURCES)

function tissue_parameters(;precursor_diffusion=.22,biomass_diffusion=.045,
        excitation_diffusion=.18,barrier=3.,taxis=0.,taxis_sensitivity=8.,growth=.45,capacity=1.2,
        turnover=.008,starvation_turnover=.04,construction=.075,erosion=.004,
        metabolic_work=.055,half_saturation=.15,tissue_saturation=.12,
        neural_gain=1.85,inhibition_gain=1.65,inhibition_time=7.,
        memory_time=45.,memory_feedback=.35,matrix_feedback=.28,
        activity_control=1.,constant_activity=.5)
    q=(;precursor_diffusion,biomass_diffusion,excitation_diffusion,barrier,taxis,taxis_sensitivity,growth,
       capacity,turnover,starvation_turnover,construction,erosion,metabolic_work,
       half_saturation,tissue_saturation,neural_gain,inhibition_gain,
       inhibition_time,memory_time,memory_feedback,matrix_feedback,activity_control,constant_activity)
    all(x->isfinite(x)&&x>=0,values(q)) || throw(ArgumentError("nonnegative finite tissue parameters required"))
    all(>(0),(capacity,half_saturation,tissue_saturation,inhibition_time,memory_time,metabolic_work)) ||
        throw(ArgumentError("tissue scales and work rate must be positive"))
    0<=activity_control<=1 || throw(ArgumentError("activity_control must lie in [0,1]"))
    0<=constant_activity<=1 || throw(ArgumentError("constant_activity must lie in [0,1]"))
    p=map(Float32,q)
    all(isfinite,values(p)) && all(>(0),(p.capacity,p.half_saturation,p.tissue_saturation,
        p.inhibition_time,p.memory_time,p.metabolic_work)) || throw(ArgumentError("parameters must remain valid Float32 values"))
    p
end

mutable struct LivingTissue{A,Q,G}
    state::A
    scratch::A
    buffers::Tuple{A,A}
    parameters::Q
    grid::G
end
tissue_world(grid;backend=CPU(),seed=911)=World(grid;backend,seed,fields=(),scratch_fields=())

function LivingTissue(w::World;initial=nothing,parameters=tissue_parameters(),seed=w.seed,pattern="ribbons")
    q=tissue_parameters(;parameters...)
    host=initial===nothing ? tissue_initial_state(w.grid,seed;pattern) : Array{Float32}(initial)
    size(host)==(w.grid.size...,7) || throw(DimensionMismatch("seven tissue channels required"))
    all(isfinite,host) && minimum(view(host,:,:,:,1:4))>=0 || throw(ArgumentError("invalid physical tissue state"))
    maximum(abs,view(host,:,:,:,5:7))<=1 || throw(ArgumentError("neural tissue state must lie in [-1,1]"))
    a=KernelAbstractions.zeros(w.backend,Float32,size(host)...);b=similar(a)
    copyto!(a,host);fill!(b,0)
    p=LivingTissue(a,b,(a,b),q,w.grid);validate(w,p);p
end

"""Deliberately seeded curved strands, not a claim of spontaneously grown anatomy."""
function tissue_initial_state(grid,seed;pattern="ribbons")
    pattern in ("ribbons","lace") || throw(ArgumentError("unknown tissue seed pattern"))
    rng=MersenneTwister(seed);dims=grid.size
    host=zeros(Float32,dims...,7)
    for k in 1:dims[3],j in 1:dims[2],i in 1:dims[1]
        x=(i-.5f0)/dims[1]-.5f0;y=(j-.5f0)/dims[2]-.5f0;z=dims[3]==1 ? 0f0 : (k-.5f0)/dims[3]-.5f0
        # A broken pair of winding ribbons lets fronts collide and seal gaps.
        ribbon=0f0
        for arm in (-1f0,1f0)
            center=arm*(.10f0+.075f0*sin(12f0*x))
            distance=(y-center)^2+(dims[3]==1 ? 0f0 : (z-.09f0*cos(10f0*x+arm))^2)
            gate=abs(x)<.27f0 && !(arm>0 && -.025f0<x<.065f0)
            gate && (ribbon=max(ribbon,exp(-distance/.0007f0)))
        end
        if pattern=="lace"
            ribbon=0f0
            for arm in -2:2
                center=.145f0*arm+.037f0*sin(15f0*x+1.7f0*arm)+.014f0*sin(37f0*x-.7f0*arm)
                distance=(y-center)^2+(dims[3]==1 ? 0f0 : (z-.08f0*sin(13f0*x+arm))^2)
                gate=abs(x)<.40f0 && cos(31f0*x+2f0*arm)>-.88f0
                gate && (ribbon=max(ribbon,exp(-distance/.000075f0)))
            end
            for arm in (-1f0,1f0)
                center=arm*(.13f0+.028f0*cos(14f0*y))
                distance=(x-center)^2+(dims[3]==1 ? 0f0 : (z+.08f0*cos(10f0*y))^2)
                ribbon=max(ribbon,.8f0*exp(-distance/.00006f0)*Float32(abs(y)<.36f0))
            end
        end
        precursor=.8f0+.3f0*cos(6.2831855f0*x)*cos(12.566371f0*y)+.12f0*cos(18.849556f0*z+5f0*x)
        host[i,j,k,1]=precursor
        host[i,j,k,2]=.42f0*ribbon
        host[i,j,k,5]=.10f0*ribbon*sin(19f0*x+11f0*y)+.012f0*(rand(rng,Float32)-.5f0)*ribbon
    end
    host
end

function validate(w::World,p::LivingTissue)
    w.grid===p.grid || throw(ArgumentError("tissue world grid was replaced"))
    size(p.state)==size(p.scratch)==(w.grid.size...,7) || throw(DimensionMismatch("tissue buffer shape"))
    p.state!==p.scratch && !Base.mightalias(p.state,p.scratch) || throw(ArgumentError("distinct tissue buffers required"))
    all(a->any(b->a===b,p.buffers),(p.state,p.scratch)) || throw(ArgumentError("tissue buffers were replaced"))
    all(a->KernelAbstractions.get_backend(a)==w.backend,p.buffers) || throw(ArgumentError("tissue backend mismatch"))
    nothing
end
function stability_limit(w,p::LivingTissue)
    q=p.parameters
    spatial=sum((inv(Float64(w.grid.spacing[d])^2) for d in 1:3 if w.grid.size[d]>1);init=0.)
    # The physical transport update is a convex combination. Excitation adds a
    # bounded relaxation with coefficient <=1, preserving [-1,1] under this bound.
    min(.08, .9/(1+2spatial*max(q.precursor_diffusion,q.biomass_diffusion+q.taxis,q.excitation_diffusion)),
        .9q.inhibition_time,.9q.memory_time)
end
@inline tissue_fraction(rate,dt)=(dt*rate)/(1f0+dt*rate)
@inline tissue_support(b,q)=b/(q.tissue_saturation+b)
@inline function tissue_taxis(s,i,j,k,h,faces,q)
    nx,ny,nz=size(s,1),size(s,2),size(s,3)
    im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
    jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
    km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
    n=s[i,j,k,1];b=s[i,j,k,2];flux=0f0
    # The same antisymmetric flux is evaluated from either side of a face.
    # Donor selection is upwind; |tanh|<=1 gives the outgoing-rate bound used
    # in stability_limit, including ordinary biomass diffusion.
    if nx>1
        for (a,face) in ((im,faces[1]),(ip,faces[2]))
            drive=tanh(q.taxis_sensitivity*(s[a,j,k,1]-n))
            flux+=face*(max(0f0,-drive)*s[a,j,k,2]-max(0f0,drive)*b)/h[1]^2
        end
    end
    if ny>1
        for (a,face) in ((jm,faces[3]),(jp,faces[4]))
            drive=tanh(q.taxis_sensitivity*(s[i,a,k,1]-n))
            flux+=face*(max(0f0,-drive)*s[i,a,k,2]-max(0f0,drive)*b)/h[2]^2
        end
    end
    if nz>1
        for (a,face) in ((km,faces[5]),(kp,faces[6]))
            drive=tanh(q.taxis_sensitivity*(s[i,j,a,1]-n))
            flux+=face*(max(0f0,-drive)*s[i,j,a,2]-max(0f0,drive)*b)/h[3]^2
        end
    end
    flux
end
@inline function tissue_laplacian(s,i,j,k,c,h,faces)
    nx,ny,nz=size(s,1),size(s,2),size(s,3)
    im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
    jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
    km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
    v=s[i,j,k,c];lap=0f0
    nx>1 && (lap+=(faces[1]*(s[im,j,k,c]-v)+faces[2]*(s[ip,j,k,c]-v))/h[1]^2)
    ny>1 && (lap+=(faces[3]*(s[i,jm,k,c]-v)+faces[4]*(s[i,jp,k,c]-v))/h[2]^2)
    nz>1 && (lap+=(faces[5]*(s[i,j,km,c]-v)+faces[6]*(s[i,j,kp,c]-v))/h[3]^2)
    lap
end

@kernel function living_tissue_kernel!(out,s,h,dt,q)
    i,j,k=@index(Global,NTuple)
    @inbounds begin
        nx,ny,nz=size(s,1),size(s,2),size(s,3)
        im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
        jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
        km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
        n=s[i,j,k,1];b=s[i,j,k,2];m=s[i,j,k,3];waste=s[i,j,k,4]
        x=s[i,j,k,5];y=s[i,j,k,6];memory=s[i,j,k,7]
        # Reuse RecurrentMaterial's symmetric face constitutive primitive.
        faces=(Interior.recurrent_face(q.barrier*m,q.barrier*s[im,j,k,3]),
            Interior.recurrent_face(q.barrier*m,q.barrier*s[ip,j,k,3]),
            Interior.recurrent_face(q.barrier*m,q.barrier*s[i,jm,k,3]),
            Interior.recurrent_face(q.barrier*m,q.barrier*s[i,jp,k,3]),
            Interior.recurrent_face(q.barrier*m,q.barrier*s[i,j,km,3]),
            Interior.recurrent_face(q.barrier*m,q.barrier*s[i,j,kp,3]))
        support=tissue_support(b,q)
        neural_faces=(faces[1]*support*tissue_support(s[im,j,k,2],q),
            faces[2]*support*tissue_support(s[ip,j,k,2],q),
            faces[3]*support*tissue_support(s[i,jm,k,2],q),
            faces[4]*support*tissue_support(s[i,jp,k,2],q),
            faces[5]*support*tissue_support(s[i,j,km,2],q),
            faces[6]*support*tissue_support(s[i,j,kp,2],q))
        # Transport first, then use the transported donor amounts for every
        # local reaction. No clipping can hide a negative transport inventory.
        n+=dt*q.precursor_diffusion*tissue_laplacian(s,i,j,k,1,h,faces)
        b+=dt*q.biomass_diffusion*tissue_laplacian(s,i,j,k,2,h,faces)
        q.taxis>0 && (b+=dt*q.taxis*tissue_taxis(s,i,j,k,h,faces,q))
        requested=dt*q.metabolic_work*b*(.25f0+abs(x))
        work=min(n,requested);n-=work;waste+=work
        paid=requested>0 ? work/requested : 0f0
        fuel=n/(q.half_saturation+n)
        active=(1f0-q.activity_control)*q.constant_activity+q.activity_control*max(0f0,x)
        growth=min(n,b*tissue_fraction(q.growth*paid*fuel*max(0f0,1f0-b/q.capacity)*(.18f0+active),dt))
        n-=growth;b+=growth
        death=b*tissue_fraction(q.turnover+q.starvation_turnover*(1f0-fuel),dt)
        b-=death;n+=death
        built=b*tissue_fraction(q.construction*paid*active,dt)
        b-=built;m+=built
        eroded=m*tissue_fraction(q.erosion*(1f0+2f0*max(0f0,-x)),dt)
        m-=eroded;n+=eroded
        # Bounded neural targets, delayed inhibition, and passive slow memory.
        # Paid excitation changes construction/growth above; matrix changes
        # thresholds here and the symmetric transport coefficients above.
        target=tanh(q.neural_gain*x-q.inhibition_gain*y-q.memory_feedback*memory-
            q.matrix_feedback*m/(q.half_saturation+m))
        vitality=tissue_support(b,q)*paid*fuel
        dx=vitality*(target-x)-(1f0-tissue_support(b,q))*x+
            q.excitation_diffusion*tissue_laplacian(s,i,j,k,5,h,neural_faces)
        out[i,j,k,1]=n;out[i,j,k,2]=b;out[i,j,k,3]=m;out[i,j,k,4]=waste
        out[i,j,k,5]=x+dt*dx
        out[i,j,k,6]=y+dt*(x-y)/q.inhibition_time
        out[i,j,k,7]=memory+dt*(x-memory)/q.memory_time
    end
end
function advance!(w::World,p::LivingTissue,dt::Real)
    Interior.check_dt(w,p,dt)
    living_tissue_kernel!(w.backend,256)(p.scratch,p.state,w.grid.spacing,Float32(dt),p.parameters;ndrange=w.grid.size)
    p.state,p.scratch=p.scratch,p.state
    w
end
record_process_fields(p::LivingTissue)=(state=p.state,parameters=p.parameters,channels=CHANNELS)
owned_material_mass(p::LivingTissue)=sum(Float64,view(Array(p.state),:,:,:,1:4))
process_ports(::LivingTissue)=ProcessPorts(owned=(:state,),drives="finite initial precursor converted to zero-potential spent material; no replenishment",
    memory="immobile paid matrix, site-bound excitation/inhibition and slower memory",
    law="symmetric face transport; bounded precursor/biomass/matrix transfers; paid delayed-inhibition neural activity")

function tissue_statistics(w,p;state=Array(p.state))
    n,b,m,spent,x,y,memory=(view(state,:,:,:,c) for c in 1:7)
    (time=w.time,material=sum(Float64,n)+sum(Float64,b)+sum(Float64,m)+sum(Float64,spent),
     precursor=sum(Float64,n),biomass=sum(Float64,b),matrix=sum(Float64,m),spent=sum(Float64,spent),
     potential=sum(Float64,n)+sum(Float64,b)+sum(Float64,m),
     excitation_std=length(x)>1 ? std(x) : 0f0,excitation_max=maximum(x),excitation_min=minimum(x),
     memory_std=length(memory)>1 ? std(memory) : 0f0,neural_maxabs=maximum(abs,view(state,:,:,:,5:7)),
     physical_min=minimum(view(state,:,:,:,1:4)),finite=all(isfinite,state),
     active_sites=count(>(.15f0),b),matrix_sites=count(>(.04f0),m),
     mean_permeability=mean(1f0 ./ (1f0 .+ 8f0*p.parameters.barrier.*m)))
end

function tissue_checkpoint(w,p,path)
    ispath(path) && throw(ArgumentError("checkpoint exists"))
    serialize(path,(format="living-tissue-v1",source_sha256=MODEL_SHA,
        shared_source_sha256=SHARED_SHA,
        dims=w.grid.size,extent=w.grid.extent,time=w.time,seed=w.seed,
        parameters=p.parameters,state=Array(p.state)))
end
function restore_tissue!(path;backend=CPU())
    c=deserialize(path)
    c.format=="living-tissue-v1" && c.source_sha256==MODEL_SHA || throw(ArgumentError("checkpoint requires its pinned study source"))
    c.shared_source_sha256==SHARED_SHA || throw(ArgumentError("shared grid/stepping or recurrent constitutive source changed"))
    w=tissue_world(Grid(c.dims...;extent=c.extent);backend,seed=c.seed)
    p=LivingTissue(w;initial=c.state,parameters=c.parameters);w.time=c.time
    w,p
end

include("living_tissue_runner.jl")
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3) || error("usage: living_tissue.jl RECIPE.toml FRESH_OUTPUT [cpu|amd|metal]")
    LivingTissueStudy.run_living_tissue(ARGS[1],ARGS[2];backend_kind=length(ARGS)==3 ? ARGS[3] : "cpu")
end
