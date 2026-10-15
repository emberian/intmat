# Volume-filling constitutive option with paid growth and recipient-limited
# exchange. The original turnover law/specimens remain separate. Solid
# occupancy is B+M; precursor, reserve and spent material are dissolved species.
if !isdefined(@__MODULE__,:LivingTissueTurnover)
    include("living_tissue_turnover.jl")
end
module LivingTissuePacking
using Interior,SHA,Serialization
using Interior:KernelAbstractions
using KernelAbstractions:@kernel,@index
using ..LivingTissueTurnover
import Interior:advance!,validate,stability_limit,process_children,process_ports,
    record_process_fields,owned_material_mass
const LT=LivingTissueTurnover
const BASE=LT.LivingTissueStudy
const PACKING_SOURCE=read(@__FILE__,String)
const PACKING_SHA=bytes2hex(sha256(PACKING_SOURCE))
export PackingTissue,packing_parameters,packing_base_parameters,packing_statistics,
    packing_checkpoint,restore_packing!

function packing_parameters(;capacity=.9,growth=1.2,basal_activity=.18,
        basal_transport=0.,basal_remodeling=0.,kwargs...)
    q=turnover_parameters(;kwargs...)
    all(x->isfinite(x)&&x>=0,(capacity,growth,basal_activity,basal_transport,basal_remodeling)) && capacity>0 && basal_activity>0 ||
        throw(ArgumentError("positive capacity and basal activity required"))
    p=merge(q,map(Float32,(;capacity,growth,basal_activity,basal_transport,basal_remodeling)))
    all(isfinite,values(p)) && p.capacity>0 && p.basal_activity>0 ||
        throw(ArgumentError("invalid Float32 packing parameters"))
    p
end
packing_base_parameters(;kwargs...)=turnover_base_parameters(;
    growth=0.,biomass_diffusion=0.,capacity=.9,kwargs...)
struct PackingTissue{P,Q}
    flow::P
    parameters::Q
end
# Both laws expose the same native tissue/reserve ports to the recording tools.
function Base.getproperty(p::PackingTissue,name::Symbol)
    name===:tissue && return getfield(p,:flow).tissue
    name===:reserve && return getfield(p,:flow).reserve
    getfield(p,name)
end
Base.propertynames(::PackingTissue)=(:flow,:parameters,:tissue,:reserve)
function PackingTissue(w;parameters=packing_parameters(),
        base_parameters=packing_base_parameters(capacity=parameters.capacity),
        reserve=12f0,initial=nothing,noise_model="bandlimited")
    q=packing_parameters(;parameters...)
    base_parameters.growth==base_parameters.biomass_diffusion==0 ||
        throw(ArgumentError("packing stage owns all growth and biomass exchange"))
    base_parameters.capacity==q.capacity || throw(ArgumentError("matching declared capacities required"))
    phase=turnover_parameters(;NamedTuple(k=>getproperty(q,k) for k in keys(turnover_parameters()))...)
    flow=TurnoverTissue(w;parameters=phase,base_parameters,reserve,initial,noise_model)
    state=Array(flow.tissue.state)
    maximum(state[:,:,:,2].+state[:,:,:,3])<=q.capacity*(1f0+16eps(Float32)) ||
        throw(ArgumentError("initial solid occupancy exceeds capacity"))
    p=PackingTissue(flow,q);validate(w,p);p
end
function validate(w,p::PackingTissue)
    validate(w,p.flow)
    p.flow.tissue.parameters.growth==p.flow.tissue.parameters.biomass_diffusion==0 ||
        throw(ArgumentError("base growth/biomass exchange must remain disabled"))
    p.flow.tissue.parameters.capacity==p.parameters.capacity || throw(ArgumentError("capacity changed"))
    nothing
end
function stability_limit(w,p::PackingTissue)
    spatial=sum((inv(Float64(w.grid.spacing[d])^2) for d in 1:3 if w.grid.size[d]>1);init=0.)
    min(stability_limit(w,p.flow.tissue),.9/(1+2spatial*max(p.parameters.motility,p.parameters.inhibition_diffusion)))
end
@inline function packing_funding(s,r,i,j,k,dt,q)
    n=s[i,j,k,1];b=s[i,j,k,2];m=s[i,j,k,3];x=s[i,j,k,5]
    addition=min(r[i,j,k],max(0f0,q.target_precursor-n)*BASE.tissue_fraction(q.renewal,dt))
    positive=max(0f0,x-q.phase_threshold);negative=max(0f0,-x-q.phase_threshold)
    requested=dt*q.phase_work*(b+m)*(q.basal_activity+positive+negative)
    work=min(n+addition,requested)
    paid=requested>0 ? work/requested : 0f0
    (addition,work,paid,positive,negative)
end
@kernel function packing_kernel!(out,s,rout,r,h,dt,q,barrier,half_saturation)
    i,j,k=@index(Global,NTuple)
    @inbounds begin
        nx,ny,nz=size(r)
        im=i==1 ? nx : i-1;ip=i==nx ? 1 : i+1
        jm=j==1 ? ny : j-1;jp=j==ny ? 1 : j+1
        km=k==1 ? nz : k-1;kp=k==nz ? 1 : k+1
        funding=packing_funding(s,r,i,j,k,dt,q)
        addition,work,paid,positive,negative=funding
        oldb=s[i,j,k,2];oldm=s[i,j,k,3];oldn=s[i,j,k,1]
        gap=max(0f0,1f0-(oldb+oldm)/q.capacity)
        donor=paid*(q.basal_transport+negative)/(1f0+q.basal_transport)
        flux=0f0;inhibition_lap=0f0
        for (a,c,d,axis) in ((im,j,k,1),(ip,j,k,1),(i,jm,k,2),(i,jp,k,2),(i,j,km,3),(i,j,kp,3))
            extent=axis==1 ? nx : axis==2 ? ny : nz
            if extent>1
                face=Interior.recurrent_face(barrier*oldm,barrier*s[a,c,d,3])
                inhibition_lap+=face*(s[a,c,d,6]-s[i,j,k,6])/h[axis]^2
                if q.motility>0
                    incoming=packing_funding(s,r,a,c,d,dt,q)
                    neighbor_donor=incoming[3]*(q.basal_transport+incoming[5])/(1f0+q.basal_transport)
                    neighbor_gap=max(0f0,1f0-(s[a,c,d,2]+s[a,c,d,3])/q.capacity)
                    drive=tanh(q.taxis_sensitivity*(s[a,c,d,1]-oldn))
                    flux+=.5f0*face*(neighbor_donor*s[a,c,d,2]*gap*(1f0-drive)-
                        donor*oldb*neighbor_gap*(1f0+drive))/h[axis]^2
                end
            end
        end
        n=oldn+addition-work;b=oldb+dt*q.motility*flux;m=oldm
        fuel=n/(half_saturation+n)
        free=max(0f0,q.capacity-(b+m))
        grown=min(n,free,b*BASE.tissue_fraction(q.growth*paid*fuel*(q.basal_activity+positive),dt))
        n-=grown;b+=grown
        built=b*BASE.tissue_fraction(q.construction*paid*positive,dt)
        b-=built;m+=built
        dissolved=m*BASE.tissue_fraction(q.unbuilding*paid*(q.basal_remodeling+negative),dt)
        m-=dissolved;n+=dissolved
        withdrawn=b*BASE.tissue_fraction(q.withdrawal*paid*negative,dt)
        b-=withdrawn;n+=withdrawn
        out[i,j,k,1]=n;out[i,j,k,2]=b;out[i,j,k,3]=m;out[i,j,k,4]=s[i,j,k,4]+work
        out[i,j,k,5]=s[i,j,k,5]
        out[i,j,k,6]=s[i,j,k,6]+dt*q.inhibition_diffusion*inhibition_lap
        out[i,j,k,7]=s[i,j,k,7]
        rout[i,j,k]=r[i,j,k]-addition
    end
end
function advance!(w::World,p::PackingTissue,dt::Real)
    Interior.check_dt(w,p,dt)
    f=p.flow
    advance!(w,f.tissue,dt)
    packing_kernel!(w.backend,256)(f.tissue.scratch,f.tissue.state,
        f.reserve_scratch,f.reserve,w.grid.spacing,Float32(dt),p.parameters,
        f.tissue.parameters.barrier,f.tissue.parameters.half_saturation;ndrange=w.grid.size)
    f.tissue.state,f.tissue.scratch=f.tissue.scratch,f.tissue.state
    f.reserve,f.reserve_scratch=f.reserve_scratch,f.reserve
    w
end
process_children(p::PackingTissue)=(p.flow,)
process_ports(::PackingTissue)=ProcessPorts(owned=(),
    drives="finite reserve plus paid basal/phase operating work; donor and solid-space limited growth",
    memory="child matrix, signal memory and immobile finite reserve",
    law="volume-filling B+M; recipient free fractions limit conservative hopping; all phase actions require actual paid fraction")
record_process_fields(p::PackingTissue)=(state=p.tissue.state,reserve=p.reserve,
    tissue_parameters=p.tissue.parameters,packing_parameters=p.parameters,channels=BASE.CHANNELS)
owned_material_mass(p::PackingTissue)=owned_material_mass(p.flow)
function packing_statistics(w,p;state=Array(p.tissue.state),reserve=Array(p.reserve))
    s=turnover_statistics(w,p.flow;state,reserve)
    maximum_solid=maximum(view(state,:,:,:,2).+view(state,:,:,:,3))
    merge(s,(capacity=p.parameters.capacity,maximum_occupancy=maximum_solid/p.parameters.capacity,
        capacity_excess=max(0f0,maximum_solid-p.parameters.capacity)))
end
function packing_checkpoint(w,p,path)
    ispath(path)&&throw(ArgumentError("checkpoint exists"))
    serialize(path,(format="living-tissue-packing-v1",packing_source=PACKING_SHA,
        turnover_source=LT.TURNOVER_SHA,base_source=BASE.MODEL_SHA,shared_sources=BASE.SHARED_SHA,
        dims=w.grid.size,extent=w.grid.extent,time=w.time,seed=w.seed,parameters=p.parameters,
        base_parameters=p.tissue.parameters,state=Array(p.tissue.state),reserve=Array(p.reserve)))
end
function restore_packing!(path;backend=CPU())
    c=deserialize(path)
    c.format=="living-tissue-packing-v1" && c.packing_source==PACKING_SHA &&
        c.turnover_source==LT.TURNOVER_SHA && c.base_source==BASE.MODEL_SHA && c.shared_sources==BASE.SHARED_SHA ||
        throw(ArgumentError("packing checkpoint requires exact pinned sources"))
    w=tissue_world(Grid(c.dims...;extent=c.extent);backend,seed=c.seed)
    p=PackingTissue(w;parameters=c.parameters,base_parameters=c.base_parameters,initial=c.state,reserve=c.reserve)
    w.time=c.time;w,p
end
end
