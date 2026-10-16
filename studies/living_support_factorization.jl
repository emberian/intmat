module LivingSupportFactorization
using Interior,KernelAbstractions
# Extra candidate scratch is study-owned; compressed production support and all
# consumers retain their original layout, candidate order and Float32 arithmetic.
@kernel function candidates!(s,e,q,g,sites,weights)
    i,c=@index(Global,NTuple)
    if e.status[1]==0 && s.active[i]==1
        t=c-1;dx=rem(t,g.candidates[1]);t=div(t,g.candidates[1]);dy=rem(t,g.candidates[2]);dz=div(t,g.candidates[2])
        lo=ntuple(d->d>g.nd ? Int32(1) : Int32(floor((s.position[i,d]-q.exchange_radius)/g.spacing[d]+.5f0)),3)
        pos=(Interior.dc_period(lo[1]+dx,g.dims[1]),Interior.dc_period(lo[2]+dy,g.dims[2]),Interior.dc_period(lo[3]+dz,g.dims[3]))
        r2=0f0
        for d in 1:Int(g.nd)
            delta=abs((pos[d]-.5f0)*g.spacing[d]-s.position[i,d]);delta=min(delta,g.extent[d]-delta);r2+=delta^2
        end
        inside=r2<q.exchange_radius^2
        sites[i,c]=inside ? Interior.dc_index(pos...,g.dims) : Int32(0)
        weights[i,c]=inside ? (1f0-r2/q.exchange_radius^2)^2 : 0f0
    end
end
@kernel function compress!(s,e,g,sites,weights,totals)
    i=@index(Global,Linear)
    if e.status[1]==0
        s.support_count[i]=0
        if s.active[i]==1
            count=Int32(0);total=0f0
            for c in 1:size(sites,2)
                site=sites[i,c]
                if site!=0
                    count+=Int32(1);weight=weights[i,c]
                    s.sites[i,count]=site;s.support[i,count]=weight;total+=weight
                end
            end
            if count==0
                pos=ntuple(d->d>g.nd ? Int32(1) : min(g.dims[d],Int32(floor(s.position[i,d]/g.spacing[d]))+1),3)
                count=Int32(1);s.sites[i,1]=Interior.dc_index(pos...,g.dims);s.support[i,1]=1f0;total=1f0
            end
            totals[i]=total;s.support_count[i]=count
        end
    end
end
@kernel function normalize!(s,e,totals)
    i,c=@index(Global,NTuple)
    if e.status[1]==0 && s.active[i]==1 && c<=s.support_count[i]
        s.support[i,c]/=totals[i]
    end
end
struct Carrier{P,S,W,T}
    process::P
    sites::S
    weights::W
    totals::T
end
function Carrier(p)
    hasproperty(p.state,:sites)||error("materialized support required")
    sites=similar(p.state.sites);weights=similar(p.state.support);totals=similar(p.state.support_count,Float32)
    Carrier(p,sites,weights,totals)
end
Base.getproperty(p::Carrier,k::Symbol)=k in (:process,:sites,:weights,:totals) ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::Carrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::Carrier,k,args...;ndrange)
    if k===Interior.dc_support!
        s,e,q,g=args
        Interior.device_cellular_launch!(p.process,candidates!,s,e,q,g,p.sites,p.weights;ndrange=size(p.sites))
        Interior.device_cellular_launch!(p.process,compress!,s,e,g,p.sites,p.weights,p.totals;ndrange=Int(g.capacity))
        Interior.device_cellular_launch!(p.process,normalize!,s,e,p.totals;ndrange=size(p.sites))
    else
        Interior.device_cellular_launch!(p.process,k,args...;ndrange)
    end
end
end
