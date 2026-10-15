"""Separate active-turnover law built from the frozen continuum-tissue process.

Finite local reserve supplies precursor; negative phases fund unbuilding and
relocation instead of merely waiting for passive erosion. Abstract inhibitory
signals communicate faster than excitation. The old study/pins are unchanged.
"""
module LivingTissueTurnover
using Interior, Random, Statistics, SHA, Serialization
using Interior: KernelAbstractions, JSON3
using KernelAbstractions: @kernel, @index
include("living_tissue.jl")
using .LivingTissueStudy
import Interior: advance!, validate, stability_limit, process_children,
    process_ports, record_process_fields, owned_material_mass
export TurnoverTissue, turnover_parameters, turnover_base_parameters,
    turnover_statistics, turnover_checkpoint, restore_turnover!, tissue_world,turnover_initial_state
const TURNOVER_SOURCE=read(@__FILE__,String)
const TURNOVER_SHA=bytes2hex(sha256(TURNOVER_SOURCE))

function turnover_parameters(;renewal=.12,target_precursor=.6,phase_work=.08,
        construction=.32,unbuilding=.65,withdrawal=.45,motility=.05,
        taxis_sensitivity=12.,inhibition_diffusion=.15,phase_threshold=.04)
    q=(;renewal,target_precursor,phase_work,construction,unbuilding,withdrawal,
        motility,taxis_sensitivity,inhibition_diffusion,phase_threshold)
    all(x->isfinite(x)&&x>=0,values(q)) || throw(ArgumentError("invalid turnover parameters"))
    phase_work>0 && phase_threshold<1 || throw(ArgumentError("positive phase work and threshold below one required"))
    p=map(Float32,q)
    all(isfinite,values(p)) && p.phase_work>0 || throw(ArgumentError("turnover parameters must remain valid Float32 values"))
    p
end
turnover_base_parameters(;kwargs...)=tissue_parameters(;precursor_diffusion=.12,
    biomass_diffusion=.0005,excitation_diffusion=.03,taxis=0.,barrier=8.,
    growth=.5,capacity=.9,turnover=.005,starvation_turnover=.02,
    construction=0.,erosion=0.,metabolic_work=.055,neural_gain=2.5,
    inhibition_gain=2.,inhibition_time=5.,memory_time=40.,memory_feedback=.3,
    matrix_feedback=.15,kwargs...)

mutable struct TurnoverTissue{P,A,Q}
    tissue::P
    reserve::A
    reserve_scratch::A
    reserve_buffers::Tuple{A,A}
    parameters::Q
end
function turnover_initial_state(grid,seed;noise_model="independent")
    rng=MersenneTwister(seed);host=zeros(Float32,grid.size...,7)
    if noise_model=="independent"
        host[:,:,:,1].=.6f0.*(1f0.+.02f0.*randn(rng,Float32,grid.size...))
        host[:,:,:,2].=.26f0.*(1f0.+.02f0.*randn(rng,Float32,grid.size...))
        host[:,:,:,5].=.02f0.*randn(rng,Float32,grid.size...)
    elseif noise_model=="bandlimited"
        nd=grid.size[3]==1 ? 2 : 3
        modes=NTuple{3,Int}[];seen=Set{NTuple{3,Int}}()
        while length(modes)<64
            mode=ntuple(d->d<=nd ? rand(rng,-12:12) : 0,3)
            radius=sum(abs2,mode)
            if 0<radius<=144 && !(mode in seen) && first(filter(!iszero,mode))>0
                push!(modes,mode);push!(seen,mode)
            end
        end
        amplitudes=randn(rng,64,3);phases=2pi.*rand(rng,64,3)
        normalizations=[sqrt(sum(abs2,view(amplitudes,:,c))/2) for c in 1:3]
        averages=[prod(sinc(modes[q][d]/grid.size[d]) for d in 1:nd) for q in 1:64]
        # Deposit the common continuous spectrum, including its exact cell
        # average and cell-center phase. Three in-place CPU inverse FFTs avoid
        # evaluating hundreds of cosines at every future 3D lattice site.
        spectrum=zeros(ComplexF64,grid.size...)
        plan=Interior.FFTW.plan_ifft!(spectrum;flags=Interior.FFTW.ESTIMATE)
        channels=(1,2,5);means=(.6,.26,0.);scales=(.012,.0052,.02)
        for c in 1:3
            fill!(spectrum,0)
            for mode in 1:64
                angle=phases[mode,c]+pi*sum(modes[mode][d]/grid.size[d] for d in 1:nd)
                coefficient=length(spectrum)*averages[mode]*amplitudes[mode,c]*cis(angle)/(2normalizations[c])
                positive=CartesianIndex(ntuple(d->mod(modes[mode][d],grid.size[d])+1,3))
                negative=CartesianIndex(ntuple(d->mod(-modes[mode][d],grid.size[d])+1,3))
                spectrum[positive]+=coefficient;spectrum[negative]+=conj(coefficient)
            end
            plan*spectrum
            for I in CartesianIndices(grid.size);host[I,channels[c]]=means[c]+scales[c]*real(spectrum[I]);end
        end
    else
        throw(ArgumentError("unknown initial noise model"))
    end
    host
end
function TurnoverTissue(w::World;parameters=turnover_parameters(),
        base_parameters=turnover_base_parameters(),reserve=4f0,initial=nothing,noise_model="independent")
    q=turnover_parameters(;parameters...)
    base_parameters.construction==base_parameters.erosion==base_parameters.taxis==0 ||
        throw(ArgumentError("base construction, erosion and unfunded taxis must be disabled"))
    if initial===nothing
        initial=turnover_initial_state(w.grid,w.seed;noise_model)
    end
    base=LivingTissue(w;initial,parameters=base_parameters)
    host=reserve isa Number ? fill(Float32(reserve),w.grid.size...) : Array{Float32}(reserve)
    size(host)==w.grid.size && all(x->isfinite(x)&&x>=0,host) || throw(ArgumentError("invalid finite reserve"))
    r=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...);s=similar(r);copyto!(r,host);fill!(s,0)
    p=TurnoverTissue(base,r,s,(r,s),q);validate(w,p);p
end
function validate(w::World,p::TurnoverTissue)
    validate(w,p.tissue)
    p.reserve!==p.reserve_scratch && !Base.mightalias(p.reserve,p.reserve_scratch) || throw(ArgumentError("distinct reserve buffers required"))
    all(a->any(b->a===b,p.reserve_buffers),(p.reserve,p.reserve_scratch)) || throw(ArgumentError("reserve buffers replaced"))
    all(a->size(a)==w.grid.size && KernelAbstractions.get_backend(a)==w.backend,p.reserve_buffers) || throw(ArgumentError("reserve grid/backend mismatch"))
    nothing
end
function stability_limit(w,p::TurnoverTissue)
    spatial=sum((inv(Float64(w.grid.spacing[d])^2) for d in 1:3 if w.grid.size[d]>1);init=0.)
    min(stability_limit(w,p.tissue),.9/(1+2spatial*max(p.parameters.motility,p.parameters.inhibition_diffusion)))
end
@inline function phase_funding(s,reserve,i,j,k,dt,q)
    n=s[i,j,k,1];b=s[i,j,k,2];m=s[i,j,k,3];x=s[i,j,k,5]
    addition=min(reserve[i,j,k],max(0f0,q.target_precursor-n)*LivingTissueStudy.tissue_fraction(q.renewal,dt))
    positive=max(0f0,x-q.phase_threshold);negative=max(0f0,-x-q.phase_threshold)
    requested=dt*q.phase_work*(b+m)*(positive+negative)
    work=min(n+addition,requested)
    paid=requested>0 ? work/requested : 0f0
    (addition,work,paid,positive,negative)
end
@inline function turnover_hopping(s,r,i,j,k,h,dt,q,barrier,funding)
    nx,ny,nz=size(s,1),size(s,2),size(s,3)
    im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
    jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
    km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
    b=s[i,j,k,2];n=s[i,j,k,1];m=s[i,j,k,3]
    donor=funding[3]*funding[5];flux=0f0
    for (a,c,d,axis) in ((im,j,k,1),(ip,j,k,1),(i,jm,k,2),(i,jp,k,2),(i,j,km,3),(i,j,kp,3))
        dims=axis==1 ? nx : axis==2 ? ny : nz
        if dims>1
            incoming=phase_funding(s,r,a,c,d,dt,q)
            neighbor_donor=incoming[3]*incoming[5]
            drive=tanh(q.taxis_sensitivity*(s[a,c,d,1]-n))
            face=Interior.recurrent_face(barrier*m,barrier*s[a,c,d,3])
            flux+=.5f0*face*(neighbor_donor*s[a,c,d,2]*(1f0-drive)-donor*b*(1f0+drive))/h[axis]^2
        end
    end
    flux
end
@kernel function active_turnover_kernel!(out,s,rout,reserve,h,dt,q,barrier)
    i,j,k=@index(Global,NTuple)
    @inbounds begin
        funding=phase_funding(s,reserve,i,j,k,dt,q)
        addition,work,paid,positive,negative=funding
        n=s[i,j,k,1]+addition-work;b=s[i,j,k,2];m=s[i,j,k,3];spent=s[i,j,k,4]+work
        q.motility>0 && (b+=dt*q.motility*turnover_hopping(s,reserve,i,j,k,h,dt,q,barrier,funding))
        built=b*LivingTissueStudy.tissue_fraction(q.construction*paid*positive,dt)
        b-=built;m+=built
        dissolved=m*LivingTissueStudy.tissue_fraction(q.unbuilding*paid*negative,dt)
        m-=dissolved;n+=dissolved
        withdrawn=b*LivingTissueStudy.tissue_fraction(q.withdrawal*paid*negative,dt)
        b-=withdrawn;n+=withdrawn
        nx,ny,nz=size(reserve)
        im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
        jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
        km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
        oldm=s[i,j,k,3]
        faces=(Interior.recurrent_face(barrier*oldm,barrier*s[im,j,k,3]),
            Interior.recurrent_face(barrier*oldm,barrier*s[ip,j,k,3]),
            Interior.recurrent_face(barrier*oldm,barrier*s[i,jm,k,3]),
            Interior.recurrent_face(barrier*oldm,barrier*s[i,jp,k,3]),
            Interior.recurrent_face(barrier*oldm,barrier*s[i,j,km,3]),
            Interior.recurrent_face(barrier*oldm,barrier*s[i,j,kp,3]))
        out[i,j,k,1]=n;out[i,j,k,2]=b;out[i,j,k,3]=m;out[i,j,k,4]=spent
        out[i,j,k,5]=s[i,j,k,5]
        out[i,j,k,6]=s[i,j,k,6]+dt*q.inhibition_diffusion*LivingTissueStudy.tissue_laplacian(s,i,j,k,6,h,faces)
        out[i,j,k,7]=s[i,j,k,7]
        rout[i,j,k]=reserve[i,j,k]-addition
    end
end
function advance!(w::World,p::TurnoverTissue,dt::Real)
    Interior.check_dt(w,p,dt)
    advance!(w,p.tissue,dt)
    active_turnover_kernel!(w.backend,256)(p.tissue.scratch,p.tissue.state,
        p.reserve_scratch,p.reserve,w.grid.spacing,Float32(dt),p.parameters,p.tissue.parameters.barrier;ndrange=w.grid.size)
    p.tissue.state,p.tissue.scratch=p.tissue.scratch,p.tissue.state
    p.reserve,p.reserve_scratch=p.reserve_scratch,p.reserve
    w
end
process_children(p::TurnoverTissue)=(p.tissue,)
process_ports(::TurnoverTissue)=ProcessPorts(owned=(:reserve,),
    drives="finite local high-potential reserve renews precursor; negative and positive phases spend additional precursor",
    memory="tissue child's neural/matrix state and remaining immobile reserve",
    law="tissue child followed by funded phase construction/unbuilding/withdrawal, conservative negative-phase hopping, and faster inhibitory signal exchange")
record_process_fields(p::TurnoverTissue)=(state=p.tissue.state,reserve=p.reserve,
    tissue_parameters=p.tissue.parameters,turnover_parameters=p.parameters,
    channels=LivingTissueStudy.CHANNELS)
owned_material_mass(p::TurnoverTissue)=owned_material_mass(p.tissue)+sum(Float64,Array(p.reserve))
function turnover_statistics(w,p;state=Array(p.tissue.state),reserve=Array(p.reserve))
    s=tissue_statistics(w,p.tissue;state);r=sum(Float64,reserve)
    nd=w.grid.size[3]==1 ? 2 : 3
    measure=prod(Float64(w.grid.spacing[d]) for d in 1:nd)
    merge(s,(reserve=r,material=s.material+r,potential=s.potential+r,
        physical_min=min(s.physical_min,minimum(reserve)),finite=s.finite&&all(isfinite,reserve),
        voxel_measure=measure,physical_material=(s.material+r)*measure,physical_potential=(s.potential+r)*measure,
        physical_precursor=s.precursor*measure,physical_biomass=s.biomass*measure,
        physical_matrix=s.matrix*measure,physical_spent=s.spent*measure,physical_reserve=r*measure))
end
function turnover_checkpoint(w,p,path)
    ispath(path)&&throw(ArgumentError("checkpoint exists"))
    serialize(path,(format="living-tissue-turnover-v1",study_source=TURNOVER_SHA,
        base_source=LivingTissueStudy.MODEL_SHA,shared_sources=LivingTissueStudy.SHARED_SHA,
        dims=w.grid.size,extent=w.grid.extent,time=w.time,seed=w.seed,parameters=p.parameters,
        base_parameters=p.tissue.parameters,state=Array(p.tissue.state),reserve=Array(p.reserve)))
end
function restore_turnover!(path;backend=CPU())
    c=deserialize(path)
    c.format=="living-tissue-turnover-v1" && c.study_source==TURNOVER_SHA &&
        c.base_source==LivingTissueStudy.MODEL_SHA && c.shared_sources==LivingTissueStudy.SHARED_SHA ||
        throw(ArgumentError("turnover checkpoint requires exact pinned sources"))
    w=tissue_world(Grid(c.dims...;extent=c.extent);backend,seed=c.seed)
    p=TurnoverTissue(w;parameters=c.parameters,base_parameters=c.base_parameters,initial=c.state,reserve=c.reserve)
    w.time=c.time;w,p
end
end
