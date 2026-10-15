"""Experimental conservative edge-exponential transport; no production integration.
Even periodic grids only. Permeability is frozen during each transport interval.
"""
module PairTransportEcology
using Interior, Random, Statistics, LinearAlgebra
using Interior: KernelAbstractions
using Interior.KernelAbstractions: @kernel,@index
include(joinpath(@__DIR__,"..","scripts","device_habitat_pilot_fixture.jl"))
export pair_transport!, transport_fixture, transport_step!, compare_transport

# Each matching owns disjoint pairs, so there are no floating atomics or global solve.
@kernel function pair_matching!(e,q,g,dt,axis,parity)
    site,k=@index(Global,NTuple)
    if e.status[1]==0
        pos=Interior.dc_coords(site,g.dims)
        if (pos[axis]-1)%2==parity
            next=ntuple(d->d==axis ? Interior.dc_period(pos[d]+1,g.dims[d]) : pos[d],3)
            j=Interior.dc_index(next...,g.dims)
            x=e.fields[site,k];y=e.fields[j,k]
            a=Float64(q.diffusion[k])*(Float64(e.permeability[site])+Float64(e.permeability[j]))*.5/Float64(g.spacing[axis])^2
            # expm1 retains tiny exchanges; equal and opposite flux in real arithmetic.
            flux=.5*(-expm1(-2a*Float64(dt)))*(Float64(y)-Float64(x))
            xx=Float32(Float64(x)+flux);yy=Float32(Float64(y)-flux)
            if !(isfinite(xx)&&isfinite(yy)&&xx>=0&&yy>=0)
                KernelAbstractions.@atomic e.status[1]|=Int32(8)
            else
                e.fields[site,k]=xx;e.fields[j,k]=yy
            end
        end
    end
end
@kernel function copy_transport!(e)
    i,k=@index(Global,NTuple)
    e.status[1]==0&&(e.scratch[i,k]=e.fields[i,k])
end
function pair_transport!(p,e,q,g,dt;cycles=1)
    cycles>0||error("positive splitting cycles required")
    axes=filter(d->g.dims[d]>1,1:3)
    all(d->iseven(g.dims[d]),axes)||error("pair transport requires even periodic dimensions; odd lattices are not silently approximated")
    matches=[(d,parity) for d in axes for parity in (0,1)]
    V=prod(Int.(g.dims))
    for _ in 1:cycles
        # Palindromic composition is second order for held spatial permeability.
        for (d,parity) in matches
            Interior.device_cellular_launch!(p,pair_matching!,e,q,g,Float32(dt/(2cycles)),d,parity;ndrange=(V,3))
        end
        for (d,parity) in reverse(matches)
            Interior.device_cellular_launch!(p,pair_matching!,e,q,g,Float32(dt/(2cycles)),d,parity;ndrange=(V,3))
        end
    end
    Interior.device_cellular_launch!(p,copy_transport!,e;ndrange=(V,3))
end

"""CPU reference factorization: cyclic M-matrix line solves, first order in time.
No CUDA/ROCm claim: device batched tridiagonal execution remains future work.
"""
function line_transport!(e,q,g,dt;cycles=1)
    e.fields isa Array||error("reference line solver is CPU-only; device batched implementation is not supplied")
    e.status[1]==0||return
    cycles>0||error("positive line cycles required")
    dims=Int.(g.dims);nmax=maximum(dims)
    diag=zeros(nmax);off=zeros(max(0,nmax-1));rhs=zeros(nmax);z=zeros(nmax);ids=zeros(Int,nmax)
    for _ in 1:cycles,axis in 1:3
        n=dims[axis];n==1&&continue
        others=filter(!=(axis),1:3)
        for b in 1:dims[others[2]],a in 1:dims[others[1]]
            for i in 1:n
                pos=ntuple(d->d==axis ? i : d==others[1] ? a : b,3)
                ids[i]=Interior.dc_index(pos...,g.dims)
            end
            for k in 1:3
                theta=Float64(dt)*q.diffusion[k]/(cycles*Float64(g.spacing[axis])^2)
                theta==0&&continue
                for i in 1:n
                    left=i==1 ? n : i-1;right=i==n ? 1 : i+1
                    l=(Float64(e.permeability[ids[i]])+e.permeability[ids[left]])*.5
                    r=(Float64(e.permeability[ids[i]])+e.permeability[ids[right]])*.5
                    diag[i]=1+theta*(l+r);rhs[i]=e.fields[ids[i],k]
                    i<n&&(off[i]=-theta*r)
                end
                if n==2
                    # Two periodic faces connect the same sites, as in native stencil.
                    rate=theta*(Float64(e.permeability[ids[1]])+e.permeability[ids[2]])
                    mean=(rhs[1]+rhs[2])*.5;contrast=(rhs[1]-rhs[2])*.5/(1+2rate)
                    rhs[1]=mean+contrast;rhs[2]=mean-contrast
                else
                    wrap=-theta*(Float64(e.permeability[ids[1]])+e.permeability[ids[n]])*.5
                    gamma=-diag[1];diag[1]-=gamma;diag[n]-=wrap*wrap/gamma
                    factor=ldlt!(SymTridiagonal(view(diag,1:n),view(off,1:n-1)))
                    ldiv!(factor,view(rhs,1:n));fill!(z,0.);z[1]=gamma;z[n]=wrap
                    ldiv!(factor,view(z,1:n))
                    correction=(rhs[1]+wrap*rhs[n]/gamma)/(1+z[1]+wrap*z[n]/gamma)
                    for i in 1:n;rhs[i]-=correction*z[i];end
                end
                for i in 1:n
                    v=Float32(rhs[i]);isfinite(v)&&v>=0||error("line solve positivity/finite failure")
                    e.fields[ids[i],k]=v
                end
            end
        end
    end
    copyto!(e.scratch,e.fields)
end

struct TransportCarrier{P}
    process::P
    mode::Symbol
    cycles::Int
end
Base.getproperty(p::TransportCarrier,k::Symbol)=k in (:process,:mode,:cycles) ? getfield(p,k) : getproperty(getfield(p,:process),k)
Interior.device_support_descriptor(p::TransportCarrier)=Interior.device_support_descriptor(p.process)
function Interior.device_cellular_launch!(p::TransportCarrier,kernel,args...;ndrange)
    if kernel===Interior.dc_diffuse!
        e,q,g,dt=args
        if p.mode===:lines
            line_transport!(e,q,g,dt;cycles=p.cycles)
        elseif p.mode===:pairs
            pair_transport!(p.process,e,q,g,dt;cycles=p.cycles)
        else
            # Same macroscopic biology clock; explicit transport alone subcycles.
            S=sum(inv(Float64(g.spacing[d])^2) for d in 1:3 if g.dims[d]>1;init=0.)
            count=max(1,ceil(Int,Float64(dt)*maximum(q.diffusion)*S/.45))
            for j in 1:count
                Interior.device_cellular_launch!(p.process,Interior.dc_diffuse!,e,q,g,Float32(dt/count);ndrange)
                j<count&&Interior.device_cellular_launch!(p.process,Interior.dc_commit_resources!,e;ndrange)
            end
        end
    else
        Interior.device_cellular_launch!(p.process,kernel,args...;ndrange)
    end
end
Interior.device_cellular_mechanics!(p::TransportCarrier,dt,mode::Val)=Interior.device_cellular_mechanics!(p.process,dt,mode)

function transport_fixture(;size=16,founders=8,capacity=64,seed=911,diffusion=1.,mode=:pairs,cycles=1)
    mode in (:pairs,:explicit,:lines)||error("unknown transport strategy")
    w,p=device_habitat_pilot_fixture(CPU();size,founders,capacity,seed,support_strategy=:materialized,energy_type=Float64)
    c=p.cells;q=merge(c.parameters,(diffusion=ntuple(_->Float32(diffusion),3),))
    c2=Interior.DeviceCellular(c.state,c.environment,c.backend,q,c.geometry,c.attempted_steps,
        Interior.device_cellular_contract(c.state,c.environment,q,c.geometry))
    h=p.habitat
    all(iszero,(q.source_rates...,q.sink_rates...,h.feedstock.parameters.source_rates...,h.feedstock.parameters.sink_rates...))||error("closed-only transport study")
    feed=merge(h.feedstock,(parameters=merge(h.feedstock.parameters,(diffusion=(Float32(diffusion),Float32(diffusion),0f0),)),))
    argument=Interior.device_habitat_argument(c2,h.arrays,h.parameters)
    ownership=(cells=c2,arrays=h.arrays,feedstock=feed,parameters=h.parameters,argument,array_shapes=map(Base.size,h.arrays),feed_shapes=map(Base.size,feed.arrays))
    h2=Interior.DeviceCatalyticHabitat(c2,h.arrays,feed,h.parameters,argument,ownership)
    # Declared initial concentration wave exercises transport immediately; no
    # painted compartments or unaccounted construction source is introduced.
    for site in axes(c2.environment.fields,1)
        xyz=Interior.dc_coords(site,c2.geometry.dims)
        wave=.5f0*cos(2f0*Float32(pi)*(xyz[1]-1)/size)
        for k in 1:2
            c2.environment.fields[site,k]*=1f0+(k==1 ? wave : -wave)
            h2.feedstock.arrays.fields[site,k]*=1f0+(k==1 ? wave : -wave)
        end
    end
    process=Interior.DeviceConstructiveEcology(c2,h2)
    (world=w,process,carrier=TransportCarrier(c2,mode,cycles))
end
nondiffusive_limit(q)=min(.25,Float64(q.response_time)/2,.1*min(q.response_time,q.adaptation_time))
function transport_step!(r,dt)
    p=r.process.cells;delta=Float32(dt)
    # Transport CFL is handled internally; preserve all nontransport bounds.
    nd=Float64(delta);q=p.parameters
    nd>0&&isfinite(nd)&&nd<=nondiffusive_limit(q)||error("neural clock exceeds study bound")
    Interior.device_cellular_validate(r.world,p);Interior.device_cellular_status!(p)
    p.attempted_steps+=1
    Interior.device_cellular_enqueue!(r.world,r.carrier,delta,:projected,r.process.habitat,UInt32(p.attempted_steps))
    Interior.device_cellular_status!(p);r
end
function compare_transport(;size=24,founders=8,capacity=64,diffusion=1.,dt=.01,steps=8,cycles=(1,2,4),strategies=(:explicit,:pairs,:lines))
    modes=[(:explicit,1);[(:pairs,c) for c in cycles];[(:lines,c) for c in cycles]]
    modes=filter(x->x[1] in strategies,modes)
    !isempty(modes)&&first(modes)[1]===:explicit||error("explicit control must be first")
    for (mode,c) in modes
        r=transport_fixture(;size,founders,capacity,diffusion,mode,cycles=c);transport_step!(r,dt)
    end
    rows=NamedTuple[];reference=nothing
    for (mode,c) in modes
        r=transport_fixture(;size,founders,capacity,diffusion,mode,cycles=c)
        initial=Interior.device_habitat_audit(Interior.device_habitat_snapshot(r.process))
        seconds=@elapsed for _ in 1:steps;transport_step!(r,dt);end
        x=Interior.device_habitat_snapshot(r.process);audit=Interior.device_habitat_audit(x)
        reference===nothing&&(reference=x)
        error=sqrt(sum(abs2,Float64.(x.cells.environment.fields).-reference.cells.environment.fields)/sum(abs2,Float64.(reference.cells.environment.fields)))
        rawerror=sqrt(sum(abs2,Float64.(x.feedstock.fields).-reference.feedstock.fields)/sum(abs2,Float64.(reference.feedstock.fields)))
        bounderror=sqrt(sum(abs2,Float64.(x.arrays.bound).-reference.arrays.bound)/max(eps(),sum(abs2,Float64.(reference.arrays.bound))))
        push!(rows,(mode,cycles=c,seconds,simulated_time_per_wall_second=steps*dt/seconds,field_relative_l2=error,raw_relative_l2=rawerror,bound_relative_l2=bounderror,
            material_residual=audit.material-initial.material,potential_residual=audit.potential+audit.work+audit.heat-initial.potential,audit))
    end
    (recipe=(size,founders,capacity,diffusion,dt,steps),rows)
end
end
