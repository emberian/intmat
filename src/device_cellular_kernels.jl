# Device algorithms use capacity-major arrays. No heap allocation or host callbacks.
@inline dc_wrap(x,L)=x-floor(x/L)*L
@inline dc_squash(x)=x/(1f0+abs(x))
@inline function dc_delta(s,g,i,j,d)
    d>g.nd && return 0f0
    a=min(i,j);b=max(i,j);x=s.position[b,d]-s.position[a,d];L=g.extent[d]
    delta=x-floor(x/L+.5f0)*L
    i==a ? delta : -delta
end
@inline function dc_distance(s,g,i,j)
    a=min(i,j);b=max(i,j)
    sqrt(dc_delta(s,g,a,b,1)^2+dc_delta(s,g,a,b,2)^2+dc_delta(s,g,a,b,3)^2)
end
@inline dc_index(i,j,k,dims)=i+dims[1]*((j-1)+dims[2]*(k-1))
@inline function dc_coords(site,dims)
    x=site-1;i=rem(x,dims[1])+1;x=div(x,dims[1]);j=rem(x,dims[2])+1
    (i,j,div(x,dims[2])+1)
end
@inline dc_period(i,n)=mod(i-1,n)+1
@inline function dc_bin(s,g,i)
    ntuple(d->d>g.nd ? Int32(1) : min(g.bins[d],Int32(floor(s.position[i,d]/g.extent[d]*g.bins[d]))+1),3)
end
@inline function dc_neighbor_bin(key,offset,n)
    dc_period(key+(n>2 ? offset-1 : offset),n)
end
# Exact periodic cover using the actual rebuilt maximum radius. Existing bonds
# are enumerated separately and are never constrained by this capture stencil.
@inline function dc_search_range(s,e,g,i,d)
    d>g.nd && return Int32(0)
    radius=s.radii[i]+reinterpret(Float32,e.metadata[7])+.15f0
    width=g.extent[d]/Float32(g.bins[d])
    # Cover Float32 rounding in both position→bin mapping operations as well
    # as the radius/width quotient, including exact bin-boundary contacts.
    padding=8f0*eps(Float32)*Float32(g.bins[d])
    Int32(min(Float32(g.bins[d]),ceil(radius/width+padding)))
end
@inline dc_search_count(n,r)=r>=div(n,2) ? n : 2r+1
@inline dc_neighbor_bin(key,offset,n,r)=dc_period(key+(r>=div(n,2) ? offset : offset-r),n)
@inline function dc_bond_valid(s,i,k)
    j=s.bonds[i,k]
    0<j<=length(s.active) && s.active[j]==1 && s.ids[j]==s.bond_ids[i,k]
end
@inline function dc_mix(x::UInt32)
    x=(x⊻(x>>16))*0x7feb352d;x=(x⊻(x>>15))*0x846ca68b;x⊻(x>>16)
end
@inline dc_random(seed,id,step,stream)=Float32((dc_mix(seed⊻UInt32(id)*0x9e3779b9⊻step*0x85ebca6b⊻UInt32(stream)*0xc2b2ae35)>>8)+UInt32(1))*5.9604645f-8
@inline function dc_normal(seed,id,step,stream)
    u=max(5.9604645f-8,dc_random(seed,id,step,2stream))
    sqrt(-2f0*log(u))*cos(6.2831855f0*dc_random(seed,id,step,2stream+1))
end

@kernel function dc_resource_quota!(e,q,dt)
    k=@index(Global,Linear)
    if e.status[1]==0
        e.resource_ledger[k,5]=0;e.resource_ledger[k,6]=0
        source=dt*q.source_rates[k]*q.source_totals[k]
        sink=dt*q.sink_rates[k]*q.sink_totals[k]
        available=max(0f0,e.resource_ledger_low[k,1]<0 ? prevfloat(e.resource_ledger[k,1]) : e.resource_ledger[k,1])
        e.quota[k,1]=source>0 ? min(1f0,available/(source*(1f0+8eps(Float32)))) : 0f0
        available=max(0f0,e.resource_ledger_low[k,2]<0 ? prevfloat(e.resource_ledger[k,2]) : e.resource_ledger[k,2])
        e.quota[k,2]=sink>0 ? min(1f0,available/(sink*(1f0+8eps(Float32)))) : 0f0
    end
end
@kernel function dc_diffuse!(e,q,g,dt)
    site,k=@index(Global,NTuple)
    if e.status[1]==0
        pos=dc_coords(site,g.dims);value=e.fields[site,k];flux=zero(value)
        for d in 1:3
            if g.dims[d]>1
                for direction in (-1,1)
                    neighbor=ntuple(a->a==d ? dc_period(pos[a]+direction,g.dims[a]) : pos[a],3)
                    j=dc_index(neighbor...,g.dims)
                    face=(e.permeability[site]+e.permeability[j])*.5f0
                    flux+=q.diffusion[k]*face*(e.fields[j,k]-value)/g.spacing[d]^2
                end
            end
        end
        e.scratch[site,k]=value+dt*flux
    end
end
@kernel function dc_geometry!(s,e,q,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        s.next[i]=0
        if s.active[i]==1
            volume=(s.mass[i]+s.metabolites[i,1]+s.metabolites[i,2])/q.density
            radius=Float32(g.nd==2 ? sqrt(volume/3.1415927f0) : cbrt(3f0*volume/12.566371f0))
            s.radii[i]=radius
            if !(isfinite(radius)&&0<radius<=q.radius_bound)
                KernelAbstractions.@atomic e.status[1]|=Int32(1)
            else
                KernelAbstractions.@atomic max(e.metadata[7],reinterpret(Int32,radius))
                key=dc_bin(s,g,i);bin=dc_index(key...,g.bins);s.cell_bin[i]=bin
                old=KernelAbstractions.@atomicswap e.heads[bin]=Int32(i)
                s.next[i]=old
            end
        else
            s.radii[i]=0
        end
    end
end
@inline function dc_materialized_support!(s,e,q,g,i)
    if e.status[1]==0
        s.support_count[i]=0
        if s.active[i]==1
            lo=ntuple(d->d>g.nd ? Int32(1) : Int32(floor((s.position[i,d]-q.exchange_radius)/g.spacing[d]+.5f0)),3)
            count=Int32(0);total=0f0
            for dz in 0:g.candidates[3]-1,dy in 0:g.candidates[2]-1,dx in 0:g.candidates[1]-1
                pos=(dc_period(lo[1]+dx,g.dims[1]),dc_period(lo[2]+dy,g.dims[2]),dc_period(lo[3]+dz,g.dims[3]))
                r2=0f0
                for d in 1:Int(g.nd)
                    delta=abs((pos[d]-.5f0)*g.spacing[d]-s.position[i,d]);delta=min(delta,g.extent[d]-delta)
                    r2+=delta^2
                end
                if r2<q.exchange_radius^2
                    weight=(1f0-r2/q.exchange_radius^2)^2
                    count+=1;s.sites[i,count]=dc_index(pos...,g.dims);s.support[i,count]=weight;total+=weight
                end
            end
            if count==0
                pos=ntuple(d->d>g.nd ? Int32(1) : min(g.dims[d],Int32(floor(s.position[i,d]/g.spacing[d]))+1),3)
                count=Int32(1);s.sites[i,1]=dc_index(pos...,g.dims);s.support[i,1]=1f0;total=1f0
            end
            for edge in 1:count;s.support[i,edge]/=total;end
            s.support_count[i]=count
        end
    end
end
@inline function dc_axis_support!(s,e,q,g,i)
    if e.status[1]==0
        s.support_count[i]=0;s.support_fallback[i]=0;s.support_normalization[i]=0f0
        if s.active[i]==1
            lo=ntuple(d->d>g.nd ? Int32(1) : Int32(floor((s.position[i,d]-q.exchange_radius)/g.spacing[d]+.5f0)),3)
            for d in 1:3,offset in 0:g.candidates[d]-1
                pos=dc_period(lo[d]+offset,g.dims[d]);s.axis_indices[i,offset+1,d]=pos
                delta=d>g.nd ? 0f0 : abs((pos-.5f0)*g.spacing[d]-s.position[i,d])
                delta=d>g.nd ? 0f0 : min(delta,g.extent[d]-delta)
                s.axis_distance2[i,offset+1,d]=delta^2
            end
            total=0f0;count=Int32(0)
            for z in 1:g.candidates[3],y in 1:g.candidates[2],x in 1:g.candidates[1]
                r2=(s.axis_distance2[i,x,1]+s.axis_distance2[i,y,2])+s.axis_distance2[i,z,3]
                if r2<q.exchange_radius^2
                    total+=(1f0-r2/q.exchange_radius^2)^2;count+=Int32(1)
                end
            end
            if count==0
                pos=ntuple(d->d>g.nd ? Int32(1) : min(g.dims[d],Int32(floor(s.position[i,d]/g.spacing[d]))+1),3)
                s.support_fallback[i]=dc_index(pos...,g.dims);count=Int32(1);total=1f0
            end
            s.support_count[i]=count;s.support_normalization[i]=total
        end
    end
end
@kernel function dc_support!(s,e,q,g)
    i=@index(Global,Linear)
    if hasproperty(s,:sites)
        dc_materialized_support!(s,e,q,g,i)
    else
        dc_axis_support!(s,e,q,g,i)
    end
end
@kernel function dc_sense!(s,e,q,g)
    i=@index(Global,Linear)
    if e.status[1]==0 && s.active[i]==1
        h=Int(g.hidden);base=2h
        for j in 1:h;s.inputs[i,j]=s.hidden[i,j];s.inputs[i,h+j]=0;end
        for j in 1:CELLULAR_SENSORS;s.inputs[i,base+j]=0;end
        concentration1=0f0;concentration2=0f0
        for (site,edgeweight) in device_support_edges(device_support_descriptor(s,g),i)
            weight=edgeweight/g.volume
            concentration1+=weight*Float32(e.fields[site,1]);concentration2+=weight*Float32(e.fields[site,2])
            pos=dc_coords(site,g.dims)
            for d in 1:3
                if g.dims[d]>1
                    plus=ntuple(a->a==d ? dc_period(pos[a]+1,g.dims[a]) : pos[a],3)
                    minus=ntuple(a->a==d ? dc_period(pos[a]-1,g.dims[a]) : pos[a],3)
                    a=dc_index(plus...,g.dims);b=dc_index(minus...,g.dims)
                    for k in 1:2
                        s.inputs[i,base+2+3(k-1)+d]+=weight*(e.fields[a,k]-e.fields[b,k])/(2g.spacing[d])
                    end
                end
            end
        end
        s.inputs[i,base+1]=dc_squash(concentration1);s.inputs[i,base+2]=dc_squash(concentration2)
        for j in 3:8;s.inputs[i,base+j]=dc_squash(s.inputs[i,base+j]*q.exchange_radius);end
        degree=0;tension=0f0
        for k in 1:Int(g.degree)
            if dc_bond_valid(s,i,k)
                j=s.bonds[i,k];degree+=1
                tension+=max(0f0,dc_distance(s,g,i,j)/s.rest[i,k]-1f0)
                gate=(s.actions[i,9]+1f0)*(s.actions[j,9]+1f0)*.25f0
                for a in 1:h;s.inputs[i,h+a]+=gate*s.hidden[j,a];end
            end
        end
        if degree>0
            for a in 1:h;s.inputs[i,h+a]/=Float32(degree);end
        end
        pressure=0f0;key=dc_bin(s,g,i)
        range=ntuple(d->dc_search_range(s,e,g,i,d),3)
        for dz in 0:dc_search_count(g.bins[3],range[3])-1,dy in 0:dc_search_count(g.bins[2],range[2])-1,dx in 0:dc_search_count(g.bins[1],range[1])-1
            b=(dc_neighbor_bin(key[1],dx,g.bins[1],range[1]),dc_neighbor_bin(key[2],dy,g.bins[2],range[2]),dc_neighbor_bin(key[3],dz,g.bins[3],range[3]))
            j=e.heads[dc_index(b...,g.bins)]
            while j>0
                if j!=i
                    pressure+=max(0f0,s.radii[i]+s.radii[j]-dc_distance(s,g,i,j))/max(s.radii[i],eps(Float32))
                end
                j=s.next[j]
            end
        end
        s.inputs[i,base+9]=s.metabolites[i,1]/max(s.mass[i]*q.store_capacity,eps(Float32))
        s.inputs[i,base+10]=s.metabolites[i,2]/max(s.mass[i]*q.store_capacity,eps(Float32))
        s.inputs[i,base+11]=dc_squash(Float32(s.energy[i])/max(s.mass[i],eps(Float32)))
        s.inputs[i,base+12]=dc_squash(s.mass[i]/q.division_mass)
        s.inputs[i,base+13]=Float32(degree)/max(1,g.degree);s.inputs[i,base+14]=dc_squash(tension)
        s.inputs[i,base+15]=dc_squash(pressure);s.inputs[i,base+19]=dc_squash(s.age[i]/20f0);s.inputs[i,base+20]=1
    end
end
@kernel function dc_controller!(s,e,g)
    i,row=@index(Global,NTuple)
    if e.status[1]==0 && s.active[i]==1
        value=0f0
        for col in 1:2Int(g.hidden)+CELLULAR_SENSORS
            value+=(s.weights[i,col,row]+s.marks[i,col,row])*s.inputs[i,col]
        end
        s.outputs[i,row]=tanh(value)
    end
end
@kernel function dc_regulate!(s,e,q,g,dt)
    i=@index(Global,Linear)
    if e.status[1]==0 && s.active[i]==1
        h=Int(g.hidden);total=0f0
        for j in 1:4;value=exp(2f0*s.outputs[i,h+3+j]);s.enzymes[i,j]=value;total+=value;end
        for j in 1:4;s.enzymes[i,j]/=total;end
        room=max(0f0,q.store_capacity*s.mass[i]-s.metabolites[i,1]-s.metabolites[i,2]);requestsum=zero(eltype(s.requests))
        for k in 1:2
            concentration=zero(eltype(e.fields))
            for (site,weight) in device_support_edges(device_support_descriptor(s,g),i);concentration+=weight*e.fields[site,k]/g.volume;end
            amount=dt*q.uptake_rate*s.mass[i]*s.enzymes[i,k]^q.enzyme_exponent*concentration/(q.half_saturation+concentration)
            s.requests[i,k]=amount;requestsum+=amount
        end
        if requestsum>room
            for k in 1:2;s.requests[i,k]*=room/requestsum;end
        end
        alpha=-expm1(-dt/q.response_time);beta=-expm1(-dt/q.adaptation_time)
        for row in 1:h
            target=tanh(atanh(clamp(s.outputs[i,row],-.999999f0,.999999f0))-q.adaptation_strength*s.adaptation[i,row])
            s.hidden[i,row]+=alpha*(target-s.hidden[i,row]);s.adaptation[i,row]+=beta*(s.hidden[i,row]-s.adaptation[i,row])
            activity=s.hidden[i,row]
            for col in 1:2h+CELLULAR_SENSORS
                mark=s.marks[i,col,row]
                s.marks[i,col,row]=clamp(mark+dt*(q.plasticity_rate*(activity*s.inputs[i,col]-activity^2*mark)-q.plasticity_decay*mark),-.3f0,.3f0)
            end
        end
        for row in 1:CELLULAR_ACTIONS;s.actions[i,row]=s.outputs[i,h+row];end
        for (site,weight) in device_support_edges(device_support_descriptor(s,g),i),k in 1:2
            proposal=s.requests[i,k]*weight
            KernelAbstractions.@atomic e.demand[site,k]+=proposal
        end
    end
end
@kernel function dc_allocate_sites!(e)
    site,k=@index(Global,NTuple)
    if e.status[1]==0
        demand=e.demand[site,k];before=e.fields[site,k];take=min(before,demand);after=max(0f0,before-take)
        before-after>take && (after=nextfloat(after))
        e.fields[site,k]=after;e.scale[site,k]=demand>0 ? (before-after)/demand : 0f0
    end
end
@kernel function dc_gather!(s,e,g)
    i=@index(Global,Linear)
    if e.status[1]==0 && s.active[i]==1
        for k in 1:2
            grant=zero(eltype(s.metabolites))
            for (site,weight) in device_support_edges(device_support_descriptor(s,g),i);grant+=s.requests[i,k]*weight*e.scale[site,k];end
            s.grants[i,k]=grant;s.metabolites[i,k]+=grant;dc_ledger_add!(s.ledger,s.ledger_low,i,2+k,grant)
        end
    end
end
@inline function dc_exchange_flux(metabolites,mass,actions,rate,i,j,k,dt)
    a=min(i,j);b=max(i,j);va=mass[a];vb=mass[b]
    openness=(actions[a,9]+1f0)*(actions[b,9]+1f0)*.25f0
    fraction=-expm1(-rate*dt*openness*(inv(va)+inv(vb)))
    canonical=fraction*(metabolites[a,k]*vb-metabolites[b,k]*va)/(va+vb)
    i==a ? canonical : -canonical
end
@kernel function dc_exchange_scale!(s,e,q,g,dt)
    i=@index(Global,Linear)
    @inbounds begin
    if e.status[1]==0 && s.active[i]==1
        for k in 1:2
            outgoing=zero(eltype(s.metabolites))
            for b in 1:Int(g.degree)
                dc_bond_valid(s,i,b) && (outgoing+=max(0f0,dc_exchange_flux(s.metabolites,s.mass,s.actions,q.junction_rate,i,s.bonds[i,b],k,dt)))
            end
            s.exchange_scale[i,k]=outgoing>0 ? min(1f0,s.metabolites[i,k]/outgoing) : 1f0
        end
    end
    end
end
@kernel function dc_exchange_commit!(s,e,q,g,dt)
    i=@index(Global,Linear)
    @inbounds begin
    if e.status[1]==0
        for k in 1:2
            next=s.metabolites[i,k]
            if s.active[i]==1
                change=zero(eltype(s.metabolites));transferred=zero(eltype(s.metabolites))
                for b in 1:Int(g.degree)
                    if dc_bond_valid(s,i,b)
                        j=s.bonds[i,b];flux=dc_exchange_flux(s.metabolites,s.mass,s.actions,q.junction_rate,i,j,k,dt)
                        actual=flux*(flux>0 ? s.exchange_scale[i,k] : s.exchange_scale[j,k])
                        change-=actual;transferred+=abs(actual)
                    end
                end
                next=max(0f0,next+change);dc_ledger_add!(s.ledger,s.ledger_low,i,5,.5f0*transferred)
            end
            s.metabolites_next[i,k]=next
        end
    end
    end
end
# Legacy Float32 operations remain unchanged. Wider batteries record actual endpoints.
@inline dc_energy_take(energy::Float32,request)=(energy-min(energy,request),min(energy,request))
@inline function dc_energy_take(energy::Float64,request)
    after=max(0.0,energy-Float64(request))
    after,energy-after
end
@inline function dc_energy_overflow(energy::Float32,capacity)
    overflow=max(0f0,energy-capacity)
    energy-overflow,overflow
end
@inline function dc_energy_overflow(energy::Float64,capacity)
    after=min(energy,Float64(capacity))
    after,energy-after
end
@inline dc_energy_control(energy)=Float32(energy)
@inline dc_energy_credit(energy::Float32,released,efficiency)=(energy+efficiency*released,(1f0-efficiency)*released)
@inline function dc_energy_credit(energy::Float64,released,efficiency)
    after=energy+Float64(efficiency*released)
    after-energy>Float64(released) && (after=prevfloat(after))
    after,Float64(released)-(after-energy)
end
@inline function device_anabolic_build!(::Nothing,i,a,b,mass,energy::Float64,growth_request,growth_work,potential_a,potential_b)
    build=min(2f0*min(a,b),growth_request)
    growth_work>0 && (build=min(build,Float32(energy/Float64(growth_work))))
    for _ in 1:4
        Float64(build*growth_work)<=energy && break
        build=max(0f0,prevfloat(build))
    end
    Float64(build*growth_work)>energy && (build=0f0)
    a-=build*.5f0;b-=build*.5f0;mass+=build
    energy,work=dc_energy_take(energy,build*growth_work)
    (a,b,mass,energy,work,build*.5f0*(potential_a+potential_b),build)
end
@inline device_anabolic_planned_work(planned,::Nothing,i)=planned
@inline function device_anabolic_build!(::Nothing,i,a,b,mass,energy,growth_request,growth_work,potential_a,potential_b)
    build=min(2f0*min(a,b),growth_request)
    growth_work>0 && (build=min(build,energy/growth_work))
    a-=build*.5f0;b-=build*.5f0;mass+=build;work=build*growth_work;energy-=work
    (a,b,mass,energy,work,build*.5f0*(potential_a+potential_b),build)
end
# Proportional withdrawal is defined by representable donor endpoints. In
# particular, `stock*a/stock` can exceed `a` by an ulp at full depletion. Never
# pass that negative remainder into construction, nor credit nominal rather
# than actual fuel. Keep this type-generic for Float32-only device backends.
@inline function dc_inventory_take_down(before::T,want) where {T<:AbstractFloat}
    after=T(max(zero(before),before-want))
    taken=before-after
    if taken>want
        after=nextfloat(after);taken=before-after
    end
    after,taken
end
@inline function dc_proportional_debit(a::T,b::T,request) where {T<:AbstractFloat}
    stock=a+b
    request<=0 && return (a,b,zero(a),zero(b))
    request>=stock && return (zero(a),zero(b),a,b)
    nexta,da=dc_inventory_take_down(a,min(a,request*(a/stock)))
    remaining=max(zero(request),request-da)
    da+remaining>request && (remaining=prevfloat(remaining))
    nextb,db=dc_inventory_take_down(b,min(b,remaining,request*(b/stock)))
    nexta,nextb,da,db
end
@kernel function dc_metabolism!(s,e,q,g,anabolism,dt)
    i=@index(Global,Linear)
    if e.status[1]==0
        s.waste[i]=0
        if s.active[i]==1
            mass=s.mass[i];a=s.metabolites[i,1];b=s.metabolites[i,2];energy=s.energy[i]
            stock=a+b;potential=a*q.potentials[1]+b*q.potentials[2]
            reserve=q.energy_capacity*mass*(s.actions[i,11]+1f0)*.5f0
            motor=0f0
            for d in 1:Int(g.nd);motor+=s.actions[i,d]^2;end
            planned=dt*((q.basal_rate+q.growth_work*q.growth_rate*s.enzymes[i,3])*mass+q.neural_cost*g.hidden+q.motor_cost*q.speed^2*motor)
            planned=device_anabolic_planned_work(planned,anabolism,i)
            demand=max(0f0,reserve+planned-dc_energy_control(energy))
            fuel=stock>0 && q.efficiency*potential>0 ? demand*stock/(q.efficiency*potential) : 0f0
            oxidize=min(stock,dt*q.catabolic_rate*mass*s.enzymes[i,4],fuel)
            if oxidize>0
                a,b,da,db=dc_proportional_debit(a,b,oxidize)
                oxidize=da+db
                released=da*q.potentials[1]+db*q.potentials[2];energy,respiration_heat=dc_energy_credit(energy,released,q.efficiency)
                dc_ledger_add!(s.ledger,s.ledger_low,i,2,respiration_heat);dc_ledger_add!(s.ledger,s.ledger_low,i,7,oxidize);s.waste[i]+=oxidize
            end
            maintenance=dt*(q.basal_rate*mass+q.neural_cost*g.hidden);paid=min(energy,maintenance)
            energy,paid=dc_energy_take(energy,maintenance);dc_ledger_add!(s.ledger,s.ledger_low,i,1,paid)
            if paid<maintenance
                lost=min(mass,dt*q.starvation_rate*mass*(1f0-Float32(paid/maintenance)));mass-=lost;s.waste[i]+=lost;dc_ledger_add!(s.ledger,s.ledger_low,i,8,lost)
            end
            stock=a+b;excess=max(0f0,stock-q.store_capacity*mass)
            if excess>0
                a,b,da,db=dc_proportional_debit(a,b,excess)
                dc_ledger_add!(s.ledger,s.ledger_low,i,2,da*q.potentials[1]+db*q.potentials[2]);s.waste[i]+=da+db
            end
            growth_request=dt*q.growth_rate*mass*s.enzymes[i,3]
            a,b,mass,energy,work,heat,build=device_anabolic_build!(anabolism,i,a,b,mass,energy,growth_request,q.growth_work,q.potentials[1],q.potentials[2])
            dc_ledger_add!(s.ledger,s.ledger_low,i,1,work);dc_ledger_add!(s.ledger,s.ledger_low,i,2,heat);dc_ledger_add!(s.ledger,s.ledger_low,i,6,build)
            energy,overflow=dc_energy_overflow(energy,q.energy_capacity*mass);dc_ledger_add!(s.ledger,s.ledger_low,i,2,overflow);s.age[i]+=dt
            if mass<q.minimum_mass
                s.waste[i]+=mass+a+b;dc_ledger_add!(s.ledger,s.ledger_low,i,2,a*q.potentials[1]+b*q.potentials[2]+energy)
                s.active[i]=0;mass=0;a=0;b=0;energy=0
                KernelAbstractions.@atomic e.metadata[4]+=Int32(1)
            end
            s.mass[i]=mass;s.metabolites[i,1]=max(0f0,a);s.metabolites[i,2]=max(0f0,b);s.energy[i]=max(0f0,energy)
            !(isfinite(mass)&&isfinite(energy)) && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
        end
    end
end
@kernel function dc_return_waste!(s,e,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        remaining=s.waste[i];edge=Int32(0)
        for (site,weight) in device_support_edges(device_support_descriptor(s,g),i)
            edge+=Int32(1)
            # Float32 support weights need not sum to exactly one. Wide physical
            # inventory closes the scatter at its final edge; no mass is minted
            # or discarded by the interpolation weights themselves.
            amount=eltype(s.waste)===Float64 ?
                (edge==s.support_count[i] ? remaining : min(remaining,s.waste[i]*Float64(weight))) : s.waste[i]*weight
            remaining-=amount
            KernelAbstractions.@atomic e.returns[site]+=amount
        end
    end
end
@kernel function dc_commit_waste!(e)
    i=@index(Global,Linear)
    if e.status[1]==0
        before=e.fields[i,3];request=e.returns[i];after=before+request
        if !isfinite(after)
            KernelAbstractions.@atomic e.status[1]|=Int32(8)
        else
            after-before>request && (after=prevfloat(after))
            e.fields[i,3]=after
            hi,lo=dc_inventory_pair_add(e.inert[i],e.inert_low[i],request-(after-before))
            e.inert[i]=hi;e.inert_low[i]=lo
        end
    end
end
@inline function dc_forces_at!(s,e,q,g,dt,i)
    if e.status[1]==0 && s.active[i]==1
        fx=0f0;fy=0f0;fz=0f0;rows=0f0;key=dc_bin(s,g,i)
        range=ntuple(d->dc_search_range(s,e,g,i,d),3)
        for dz in 0:dc_search_count(g.bins[3],range[3])-1,dy in 0:dc_search_count(g.bins[2],range[2])-1,dx in 0:dc_search_count(g.bins[1],range[1])-1
            bin=(dc_neighbor_bin(key[1],dx,g.bins[1],range[1]),dc_neighbor_bin(key[2],dy,g.bins[2],range[2]),dc_neighbor_bin(key[3],dz,g.bins[3],range[3]))
            j=e.heads[dc_index(bin...,g.bins)]
            while j>0
                if j!=i
                    distance=dc_distance(s,g,i,j);overlap=max(0f0,s.radii[i]+s.radii[j]-distance)
                    if overlap>0
                        denominator=max(distance,1f-4)
                        ux=dc_delta(s,g,i,j,1)/denominator;uy=dc_delta(s,g,i,j,2)/denominator;uz=dc_delta(s,g,i,j,3)/denominator
                        if distance==0
                            axis=mod(s.ids[i]+s.ids[j],g.nd)+1;sign=i<j ? 1f0 : -1f0
                            ux=axis==1 ? sign : 0f0;uy=axis==2 ? sign : 0f0;uz=axis==3 ? sign : 0f0
                        end
                        magnitude=q.repulsion*overlap;fx-=magnitude*ux;fy-=magnitude*uy;fz-=magnitude*uz
                        rows+=q.repulsion*(1f0+overlap/denominator)
                    end
                end
                j=s.next[j]
            end
        end
        for b in 1:Int(g.degree)
            if dc_bond_valid(s,i,b)
                j=s.bonds[i,b];distance=dc_distance(s,g,i,j);extension=max(0f0,distance-s.rest[i,b])
                if distance>0 && extension>0
                    magnitude=q.stiffness*extension/distance
                    fx+=magnitude*dc_delta(s,g,i,j,1);fy+=magnitude*dc_delta(s,g,i,j,2);fz+=magnitude*dc_delta(s,g,i,j,3)
                    rows+=q.stiffness
                end
            end
        end
        s.force[i,1]=fx;s.force[i,2]=fy;s.force[i,3]=fz;s.rows[i]=rows
        motor=0f0
        for d in 1:Int(g.nd);motor+=q.speed^2*s.actions[i,d]^2;end
        speed=sqrt(fx*fx+fy*fy+fz*fz)+sqrt(motor)
        if !(isfinite(fx)&&isfinite(fy)&&isfinite(fz)&&isfinite(rows)&&isfinite(speed))
            KernelAbstractions.@atomic e.status[1]|=Int32(8)
        elseif 2f0*rows*dt>.45f0 || speed*dt>.2f0*s.radii[i]
            KernelAbstractions.@atomic e.status[1]|=Int32(2)
        end
    end
end
@inline function dc_move_at!(s,e,q,g,dt,i)
    if e.status[1]==0
        for d in 1:3;s.position_next[i,d]=s.position[i,d];end
    end
    if e.status[1]==0 && s.active[i]==1
        v2=0f0;f2=0f0
        for d in 1:Int(g.nd);v2+=q.speed^2*s.actions[i,d]^2;f2+=s.force[i,d]^2;end
        linear=dt*sqrt(f2*v2);quadratic=(dt*q.motor_cost+dt^2*s.rows[i])*v2
        cost=linear+quadratic;factor=1f0
        if cost>s.energy[i]
            denominator=linear+sqrt(linear^2+4f0*quadratic*s.energy[i])
            factor=denominator>0 ? min(1f0,2f0*s.energy[i]/denominator) : 0f0
        end
        factor=Float32(factor)
        if eltype(s.energy)===Float64
            for _ in 1:4
                linear*factor+quadratic*factor^2<=s.energy[i] && break
                factor=max(0f0,prevfloat(factor))
            end
            linear*factor+quadratic*factor^2>s.energy[i] && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
        end
        s.energy[i],paid=dc_energy_take(s.energy[i],linear*factor+quadratic*factor^2);dc_ledger_add!(s.ledger,s.ledger_low,i,1,paid)
        for d in 1:Int(g.nd)
            s.position_next[i,d]=dc_wrap(s.position[i,d]+dt*(s.force[i,d]+factor*q.speed*s.actions[i,d]),g.extent[d])
        end
    end
end
@kernel function dc_forces!(s,e,q,g,dt)
    i=@index(Global,Linear)
    dc_forces_at!(s,e,q,g,dt,i)
end
@kernel function dc_move!(s,e,q,g,dt)
    i=@index(Global,Linear)
    dc_move_at!(s,e,q,g,dt,i)
end
@kernel function dc_propose_bonds!(s,e,q,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        for k in 1:Int(g.degree);s.candidates[i,k]=0;s.candidate_rest[i,k]=0;end
        if s.active[i]==1 && s.actions[i,8]>=0
            count=0
            for k in 1:Int(g.degree)
                if dc_bond_valid(s,i,k)
                    j=s.bonds[i,k];distance=dc_distance(s,g,i,j)
                    if s.actions[j,8]>=0 && distance<=(1f0+q.break_strain)*s.rest[i,k]
                        count+=1;s.candidates[i,count]=j;s.candidate_rest[i,count]=s.rest[i,k]
                    end
                end
            end
            # Existing bonds retain their slots; new proposals are the nearest
            # available candidates, with stable-ID tie-breaking.
            persistent=count;key=dc_bin(s,g,i)
            range=ntuple(d->dc_search_range(s,e,g,i,d),3)
        for dz in 0:dc_search_count(g.bins[3],range[3])-1,dy in 0:dc_search_count(g.bins[2],range[2])-1,dx in 0:dc_search_count(g.bins[1],range[1])-1
                bin=(dc_neighbor_bin(key[1],dx,g.bins[1],range[1]),dc_neighbor_bin(key[2],dy,g.bins[2],range[2]),dc_neighbor_bin(key[3],dz,g.bins[3],range[3]))
                j=e.heads[dc_index(bin...,g.bins)]
                while j>0
                    if j!=i && s.actions[j,8]>=0
                        distance=dc_distance(s,g,i,j)
                        if distance<=s.radii[i]+s.radii[j]+.15f0
                            exists=false
                            for k in 1:count;exists|=s.candidates[i,k]==j;end
                            if !exists
                                insertion=count+1
                                for k in persistent+1:count
                                    previous=s.candidates[i,k];other=dc_distance(s,g,i,previous)
                                    if distance<other || (distance==other && s.ids[j]<s.ids[previous])
                                        insertion=k;break
                                    end
                                end
                                if insertion<=g.degree
                                    newcount=min(Int(g.degree),count+1)
                                    for k in newcount:-1:insertion+1
                                        s.candidates[i,k]=s.candidates[i,k-1];s.candidate_rest[i,k]=s.candidate_rest[i,k-1]
                                    end
                                    s.candidates[i,insertion]=j;s.candidate_rest[i,insertion]=max(distance,s.radii[i]+s.radii[j]);count=newcount
                                end
                            end
                        end
                    end
                    j=s.next[j]
                end
            end
        end
    end
end
@kernel function dc_commit_bonds!(s,e,q,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        count=0
        # Candidate arrays are read-only throughout this commit.
        for k in 1:Int(g.degree)
            j=s.candidates[i,k]
            if j>0 && s.active[i]==1 && s.active[j]==1
                mutual=false
                for a in 1:Int(g.degree);mutual|=s.candidates[j,a]==i;end
                if mutual
                    count+=1;s.bonds[i,count]=j;s.bond_ids[i,count]=s.ids[j];s.rest[i,count]=s.candidate_rest[i,k]
                end
            end
        end
        for k in count+1:Int(g.degree);s.bonds[i,k]=0;s.bond_ids[i,k]=0;s.rest[i,k]=0;end
    end
end
@inline function dc_heap_less(s,e,a,b)
    e.heap_keys[a]<e.heap_keys[b] || (e.heap_keys[a]==e.heap_keys[b] && s.ids[e.heap[a]]<s.ids[e.heap[b]])
end
@inline function dc_lifecycle_impl!(s,e,q,g,step,highwater,::Val{dense}) where {dense}
    if e.status[1]==0
        limit=highwater===nothing || dense ? Int(g.capacity) : Int(highwater[1])
        if !(0<=limit<=g.capacity)
            e.status[1]|=Int32(8)
            return nothing
        end
        count=0;freecount=0;alive=0
        for i in 1:limit
            s.new_parent[i]=0
            if s.active[i]==0
                freecount+=1;e.free[freecount]=i
            else
                alive+=1
                if s.mass[i]>=q.division_mass && s.actions[i,10]>0
                    count+=1;e.heap[count]=i;e.heap_keys[count]=dc_mix(q.seed⊻UInt32(s.ids[i])*0x9e3779b9⊻step*0x85ebca6b)
                    slot=count
                    while slot>1
                        parent=div(slot,2)
                        if dc_heap_less(s,e,slot,parent)
                            a=e.heap[slot];b=e.heap_keys[slot];e.heap[slot]=e.heap[parent];e.heap_keys[slot]=e.heap_keys[parent];e.heap[parent]=a;e.heap_keys[parent]=b
                            slot=parent
                        else
                            break
                        end
                    end
                end
            end
        end
        if !dense
            for i in limit+1:min(Int(g.capacity),limit+count)
                freecount+=1;e.free[freecount]=i
            end
        end
        births=min(count,freecount);e.metadata[5]+=count-births
        if Int64(e.metadata[1])+births>=typemax(Int32)
            e.status[1]|=Int32(4)
        else
            heapcount=count
            for birth in 1:births
                parent=e.heap[1];child=e.free[birth];s.new_parent[child]=parent
                s.ids[child]=e.metadata[1];e.metadata[1]+=1;s.active[child]=1
                e.heap[1]=e.heap[heapcount];e.heap_keys[1]=e.heap_keys[heapcount];heapcount-=1;slot=1
                while 2slot<=heapcount
                    next=2slot
                    next<heapcount && dc_heap_less(s,e,next+1,next) && (next+=1)
                    if dc_heap_less(s,e,next,slot)
                        a=e.heap[slot];b=e.heap_keys[slot];e.heap[slot]=e.heap[next];e.heap_keys[slot]=e.heap_keys[next];e.heap[next]=a;e.heap_keys[next]=b;slot=next
                    else
                        break
                    end
                end
            end
            e.metadata[2]=alive+births;e.metadata[3]+=births
            if highwater!==nothing && births>0
                highwater[1]=max(highwater[1],e.free[births])
            end
        end
    end
end
# The dense oracle and production prefix share every numerical operation.
@kernel function dc_lifecycle!(s,e,q,g,step)
    thread=@index(Global,Linear)
    thread==1 && dc_lifecycle_impl!(s,e,q,g,step,nothing,Val(true))
end
@kernel function dc_lifecycle_dense!(s,e,q,g,step,highwater)
    thread=@index(Global,Linear)
    thread==1 && dc_lifecycle_impl!(s,e,q,g,step,highwater,Val(true))
end
@kernel function dc_lifecycle_prefix!(s,e,q,g,step,highwater)
    thread=@index(Global,Linear)
    thread==1 && dc_lifecycle_impl!(s,e,q,g,step,highwater,Val(false))
end
@kernel function dc_divide!(s,e,q,g,step)
    child=@index(Global,Linear)
    if e.status[1]==0
        parent=s.new_parent[child]
        if parent>0
            id=s.ids[child];h=Int(g.hidden)
            s.parents[child]=s.ids[parent];s.lineages[child]=s.lineages[parent];s.generations[child]=s.generations[parent]+1
            s.mass[parent]*=.5f0;s.mass[child]=s.mass[parent];s.energy[parent]*=.5f0;s.energy[child]=s.energy[parent]
            for k in 1:2;s.metabolites[parent,k]*=.5f0;s.metabolites[child,k]=s.metabolites[parent,k];end
            s.age[parent]=0;s.age[child]=0
            for row in 1:h+CELLULAR_ACTIONS,col in 1:2h+CELLULAR_SENSORS
                value=s.weights[parent,col,row];stream=col+(row-1)*(2h+CELLULAR_SENSORS)
                if dc_random(q.seed,id,step,stream)<q.mutation_probability
                    value+=q.mutation*dc_normal(q.seed,id,step,stream+10000)
                end
                !isfinite(value) && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
                s.weights[child,col,row]=value;s.marks[child,col,row]=s.marks[parent,col,row]
            end
            for k in 1:h
                value=s.hidden[parent,k];shift=q.developmental_partition*dc_normal(q.seed,id,step,k+20000)
                s.hidden[child,k]=clamp(value+shift,-1f0,1f0);s.hidden[parent,k]=clamp(value-shift,-1f0,1f0)
                s.adaptation[child,k]=s.adaptation[parent,k]
            end
            for k in 1:CELLULAR_ACTIONS;s.actions[child,k]=s.actions[parent,k];end
            dx=dc_normal(q.seed,id,step,30001);dy=dc_normal(q.seed,id,step,30002);dz=g.nd==3 ? dc_normal(q.seed,id,step,30003) : 0f0
            length=sqrt(dx*dx+dy*dy+dz*dz)
            if length==0;dx=1f0;length=1f0;end
            direction=(dx/length,dy/length,dz/length)
            for d in 1:3
                center=s.position[parent,d];shift=.05f0*q.exchange_radius*direction[d]
                s.position[parent,d]=d>g.nd ? 0f0 : dc_wrap(center-shift,g.extent[d])
                s.position[child,d]=d>g.nd ? 0f0 : dc_wrap(center+shift,g.extent[d])
            end
            volume=(s.mass[parent]+s.metabolites[parent,1]+s.metabolites[parent,2])/q.density
            radius=Float32(g.nd==2 ? sqrt(volume/3.1415927f0) : cbrt(3f0*volume/12.566371f0))
            s.radii[parent]=radius;s.radii[child]=radius
            for b in 1:Int(g.degree);s.bonds[child,b]=0;s.bond_ids[child,b]=0;s.rest[child,b]=0;end
            if s.actions[parent,8]>=0 && g.degree>0
                slot=0
                for b in 1:Int(g.degree)
                    if s.bonds[parent,b]==0;slot=b;break;end
                end
                if slot>0
                    s.bonds[parent,slot]=child;s.bond_ids[parent,slot]=id;s.rest[parent,slot]=2f0*radius
                    s.bonds[child,1]=parent;s.bond_ids[child,1]=s.ids[parent];s.rest[child,1]=2f0*radius
                end
            end
        end
    end
end
@kernel function dc_porosity_scatter!(s,e,g)
    i=@index(Global,Linear)
    if e.status[1]==0 && s.active[i]==1
        mass=s.mass[i]+s.metabolites[i,1]+s.metabolites[i,2]
        for (site,weight) in device_support_edges(device_support_descriptor(s,g),i)
            amount=mass*weight/g.volume
            KernelAbstractions.@atomic e.occupancy[site]+=amount
        end
    end
end
@kernel function dc_porosity_commit!(e,q)
    site=@index(Global,Linear)
    e.status[1]==0 && (e.permeability[site]=inv(1f0+q.packing_barrier*e.occupancy[site]))
end
@kernel function dc_complete_step!(e)
    i=@index(Global,Linear)
    if i==1 && e.status[1]==0;e.metadata[6]+=1;end
end
@inline dc_inventory_pair_add(hi,lo,x)=dc_pair_add(hi,lo,x)
@inline function dc_inventory_pair_add(hi,lo,x::Float64)
    high=Float32(x);low=Float32(x-Float64(high))
    hi,lo=dc_pair_add(hi,lo,high)
    dc_pair_add(hi,lo,low)
end
@inline function dc_pair_add(hi,lo,x)
    sum=hi+x;b=sum-hi;error=(hi-(sum-b))+(x-b)
    correction=lo+error;next=sum+correction
    next,correction-(next-sum)
end
@inline function dc_ledger_add!(ledger,ledger_low,i,k,amount::Float64)
    high=Float32(amount);low=Float32(amount-Float64(high))
    dc_ledger_add!(ledger,ledger_low,i,k,high)
    dc_ledger_add!(ledger,ledger_low,i,k,low)
    nothing
end
@inline function dc_ledger_add!(ledger,ledger_low,i,k,amount)
    hi,lo=dc_pair_add(ledger[i,k],ledger_low[i,k],amount)
    ledger[i,k]=hi;ledger_low[i,k]=lo
    nothing
end
@inline function dc_resource_add!(e,k,column,amount)
    before=e.resource_ledger[k,column];low=e.resource_ledger_low[k,column]
    hi,lo=dc_pair_add(before,low,amount)
    if !(isfinite(hi)&&isfinite(lo))
        KernelAbstractions.@atomic e.status[1]|=Int32(8)
    elseif (amount!=0 && hi==before && lo==low) || (column<=2 && hi+lo<0)
        KernelAbstractions.@atomic e.status[1]|=Int32(16)
    else
        e.resource_ledger[k,column]=hi;e.resource_ledger_low[k,column]=lo
    end
end
@kernel function dc_clear_grid!(e)
    i=@index(Global,Linear)
    if e.status[1]==0
        e.demand[i,1]=0;e.demand[i,2]=0;e.returns[i]=0
    end
end
@kernel function dc_clear_heads!(e)
    i=@index(Global,Linear)
    if e.status[1]==0
        e.heads[i]=0
        i==1 && (e.metadata[7]=0)
    end
end
@kernel function dc_clear_occupancy!(e)
    i=@index(Global,Linear)
    e.status[1]==0 && (e.occupancy[i]=0)
end
@kernel function dc_commit_resources!(e)
    i,k=@index(Global,NTuple)
    e.status[1]==0 && (e.fields[i,k]=e.scratch[i,k])
end
@kernel function dc_commit_metabolites!(s,e)
    i=@index(Global,Linear)
    if e.status[1]==0
        s.metabolites[i,1]=s.metabolites_next[i,1];s.metabolites[i,2]=s.metabolites_next[i,2]
    end
end
@kernel function dc_commit_positions!(s,e)
    i=@index(Global,Linear)
    if e.status[1]==0
        for d in 1:3;s.position[i,d]=s.position_next[i,d];end
    end
end
@kernel function dc_validate_complete!(s,e,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        valid=true
        if i<=size(e.fields,1)
            for k in 1:3;valid&=isfinite(e.fields[i,k]) && e.fields[i,k]>=0;end
            valid&=isfinite(e.inert[i]) && isfinite(e.inert_low[i])
        end
        if i<=g.capacity
            for k in 1:9;valid&=isfinite(s.ledger[i,k]) && isfinite(s.ledger_low[i,k]);end
            if s.active[i]==1
                valid&=isfinite(s.mass[i]) && s.mass[i]>0 && isfinite(s.energy[i]) && s.energy[i]>=0
                for k in 1:2;valid&=isfinite(s.metabolites[i,k]) && s.metabolites[i,k]>=0;end
                for d in 1:3;valid&=isfinite(s.position[i,d]);end
                for k in 1:Int(g.hidden);valid&=isfinite(s.hidden[i,k]) && isfinite(s.adaptation[i,k]);end
                for k in 1:CELLULAR_ACTIONS;valid&=isfinite(s.actions[i,k]);end
            end
        end
        !valid && (KernelAbstractions.@atomic e.status[1]|=Int32(8))
    end
end
@kernel function dc_observe_cells!(s,e,q,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        KernelAbstractions.@atomic e.observation[5]+=s.ledger[i,1]+s.ledger_low[i,1]
        KernelAbstractions.@atomic e.observation[6]+=s.ledger[i,2]+s.ledger_low[i,2]
        KernelAbstractions.@atomic e.observation[8]+=s.ledger[i,6]+s.ledger_low[i,6]
        KernelAbstractions.@atomic e.observation[9]+=s.ledger[i,7]+s.ledger_low[i,7]
        if s.active[i]==1
            material=s.mass[i]+s.metabolites[i,1]+s.metabolites[i,2]
            potential=Float32(s.energy[i])+s.metabolites[i,1]*q.potentials[1]+s.metabolites[i,2]*q.potentials[2]
            KernelAbstractions.@atomic e.observation[3]+=material
            KernelAbstractions.@atomic e.observation[4]+=potential
            bonds=0f0
            for b in 1:Int(g.degree);dc_bond_valid(s,i,b) && i<s.bonds[i,b] && (bonds+=1f0);end
            KernelAbstractions.@atomic e.observation[7]+=bonds
            KernelAbstractions.@atomic min(e.observation_min[2],reinterpret(Int32,Float32(s.energy[i])))
            !(isfinite(material)&&isfinite(potential)) && (KernelAbstractions.@atomic e.observation[10]+=1f0)
        end
    end
end
@kernel function dc_observe_grid!(e,q)
    i=@index(Global,Linear)
    if e.status[1]==0
        material=e.inert[i]+e.inert_low[i];potential=0f0;minvalue=Inf32
        for k in 1:3
            material+=e.fields[i,k];potential+=e.fields[i,k]*q.potentials[k];minvalue=min(minvalue,e.fields[i,k])
        end
        KernelAbstractions.@atomic e.observation[1]+=material
        KernelAbstractions.@atomic e.observation[2]+=potential
        KernelAbstractions.@atomic min(e.observation_min[1],reinterpret(Int32,minvalue))
        !(isfinite(material)&&isfinite(potential)) && (KernelAbstractions.@atomic e.observation[10]+=1f0)
    end
end

