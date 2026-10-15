# Normalize imported scaffold after Float64→Float32 conversion without losing
# unsupported bound material. This is the same compulsory collapse as a step.
@kernel function dh_initial_collapse!(a,e,q)
    site=@index(Global,Linear)
    if e.status[1]==0
        for k in 1:2
            enzyme=reinterpret(Float32,a.enzyme_bits[site,k]);other=reinterpret(Float32,a.enzyme_bits[site,3-k])
            capacity=dh_float_down(Float64(q.capacity_per_enzyme)*((1-Float64(q.cross_binding))*enzyme+Float64(q.cross_binding)*other))
            before=e.fields[site,k]
            free,bound,residue=dh_sorption_collapse(before,a.bound[site,k],capacity)
            e.fields[site,k]=free;a.bound[site,k]=bound
            dh_add!(a.grid_ledger,a.grid_ledger_low,site,k+8,max(0.,Float64(free)-Float64(before)))
            dh_inert_add!(a.grid_inert,a.grid_inert_low,site,residue)
            dh_add!(a.grid_ledger,a.grid_ledger_low,site,11,residue*q.food_potentials[k])
        end
    end
end
# Site-local chemistry and capacity-major inherited regulation kernels.
@kernel function dh_prepare!(a,feed,s,descriptor,e,q,cq,g,dt)
    i=@index(Global,Linear)
    if e.status[1]==0
        a.requests[i,1]=0;a.requests[i,2]=0;a.last_secretion[i,1]=0;a.last_secretion[i,2]=0
        if s.active[i]==1
            env1=0f0;env2=0f0;enzyme1=0f0;enzyme2=0f0
            for (site,support_weight) in device_support_edges(descriptor,i)
                weight=support_weight/g.volume
                env1+=weight*feed.fields[site,1];env2+=weight*feed.fields[site,2]
                enzyme1+=weight*reinterpret(Float32,a.enzyme_bits[site,1]);enzyme2+=weight*reinterpret(Float32,a.enzyme_bits[site,2])
            end
            env1/=1f0+env1;env2/=1f0+env2;enzyme1/=1f0+enzyme1;enzyme2/=1f0+enzyme2
            x=0f0;y=0f0;z=0f0;h=Int(g.hidden)
            for col in 1:h+10
                j=col-h
                input=col<=h ? s.hidden[i,col] : j==1 ? env1 : j==2 ? env2 : j==3 ? enzyme1 : j==4 ? enzyme2 :
                    j==5 ? s.metabolites[i,1]/max(s.mass[i]*cq.store_capacity,eps(Float32)) :
                    j==6 ? s.metabolites[i,2]/max(s.mass[i]*cq.store_capacity,eps(Float32)) :
                    j==7 ? Float32(s.energy[i])/(s.mass[i]+Float32(s.energy[i])) : j==8 ? s.mass[i]/(cq.division_mass+s.mass[i]) :
                    j==9 ? s.age[i]/(20f0+s.age[i]) : 1f0
                x+=a.weights[i,col,1]*input;y+=a.weights[i,col,2]*input;z+=a.weights[i,col,3]*input
            end
            investment=inv(1f0+exp(-clamp(q.regulation ? x : q.investment_bias,-40f0,40f0)))
            allocation=q.regulation ? inv(1f0+exp(-clamp(y-z,-40f0,40f0))) : .5f0
            a.expression[i,1]=investment;a.expression[i,2]=allocation;a.expression[i,3]=1f0-allocation
            a.requests[i,1]=dt*q.secretion_rate*s.mass[i]*investment*allocation^q.allocation_exponent
            a.requests[i,2]=dt*q.secretion_rate*s.mass[i]*investment*(1f0-allocation)^q.allocation_exponent
            !(isfinite(a.requests[i,1])&&isfinite(a.requests[i,2])) && dh_fail!(e.status,8)
        end
    end
end
@kernel function dh_inherit!(a,s,e,q,g,step)
    child=@index(Global,Linear)
    if e.status[1]==0
        parent=s.new_parent[child]
        if parent>0
            for row in 1:3,col in 1:Int(g.hidden)+10
                value=a.weights[parent,col,row];stream=col+(row-1)*(Int(g.hidden)+10)
                if dc_random(q.seed⊻0x831c47e2,s.ids[child],step,stream)<q.mutation_probability
                    value+=q.mutation*dc_normal(q.seed⊻0x831c47e2,s.ids[child],step,stream+10000)
                end
                !isfinite(value) && dh_fail!(e.status,8)
                a.weights[child,col,row]=value
            end
            for k in 1:3;a.expression[child,k]=0;end
            for k in 1:2;a.requests[child,k]=0;a.last_secretion[child,k]=0;end
        end
    end
end
@kernel function dh_permeability!(a,e,g,q)
    site=@index(Global,Linear)
    if e.status[1]==0
        enzyme=Float64(reinterpret(Float32,a.enzyme_bits[site,1]))+Float64(reinterpret(Float32,a.enzyme_bits[site,2]))
        factor=Float32(inv(1.0+Float64(q.transport_barrier)*enzyme/Float64(g.volume)))
        a.permeability[site]=factor;e.permeability[site]*=factor
        !(isfinite(factor)&&0<=factor<=1) && dh_fail!(e.status,64)
    end
end
@kernel function dh_react_sorb!(a,feed,e,q,fq,g,dt)
    site=@index(Global,Linear)
    if e.status[1]==0
        inert=0.;heat=0.
        for k in 1:2
            enzyme=reinterpret(Float32,a.enzyme_bits[site,k]);raw=feed.fields[site,k]
            request=Float64(dt)*q.catalytic_rate*Float64(enzyme)*Float64(raw)/(Float64(q.half_saturation)*g.volume+Float64(raw))
            next,taken=dh_take_down(raw,min(Float64(raw),request))
            food,delivered=dh_add_down(e.fields[site,k],taken)
            feed.fields[site,k]=next;e.fields[site,k]=food
            inert+=taken-delivered;heat+=taken*fq.potentials[k]-delivered*q.food_potentials[k]
            dh_add!(a.grid_ledger,a.grid_ledger_low,site,k,delivered)
            a.release[site,k]=Float32(delivered/(Float64(dt)*g.volume))
            remaining=enzyme*exp(-q.decay_rate*dt);decayed=Float64(enzyme)-Float64(remaining)
            a.enzyme_bits[site,k]=reinterpret(UInt32,remaining)
            waste,received=dh_add_down(e.fields[site,3],decayed);e.fields[site,3]=waste
            inert+=decayed-received;dh_add!(a.grid_ledger,a.grid_ledger_low,site,k+2,decayed)
        end
        for k in 1:2
            enzyme=reinterpret(Float32,a.enzyme_bits[site,k]);other=reinterpret(Float32,a.enzyme_bits[site,3-k])
            capacity=dh_float_down(Float64(q.capacity_per_enzyme)*((1-Float64(q.cross_binding))*enzyme+Float64(q.cross_binding)*other))
            free=e.fields[site,k];bound=a.bound[site,k]
            initial_free=free
            free,bound,residue=dh_sorption_collapse(free,bound,capacity)
            dh_add!(a.grid_ledger,a.grid_ledger_low,site,k+8,max(0.,Float64(free)-Float64(initial_free)))
            inert+=residue;heat+=residue*q.food_potentials[k]
            target=dh_sorption_target(bound,Float64(free)+Float64(bound),capacity,q.adsorption_rate,q.desorption_rate,g.volume,dt)
            before_free=free;before_bound=bound
            free,bound,residue=dh_sorption_exchange(free,bound,target)
            if bound>before_bound
                dh_add!(a.grid_ledger,a.grid_ledger_low,site,k+4,Float64(bound)-Float64(before_bound))
            elseif free>before_free
                dh_add!(a.grid_ledger,a.grid_ledger_low,site,k+6,Float64(free)-Float64(before_free))
            end
            e.fields[site,k]=free;a.bound[site,k]=bound
            inert+=residue;heat+=residue*q.food_potentials[k]
        end
        if !(isfinite(inert)&&inert>=0&&isfinite(heat)&&heat>=0)
            dh_fail!(e.status,64)
        else
            dh_inert_add!(a.grid_inert,a.grid_inert_low,site,inert)
            dh_add!(a.grid_ledger,a.grid_ledger_low,site,11,heat)
        end
    end
end
@kernel function dh_validate_complete!(a,feed,e,g)
    i=@index(Global,Linear)
    if e.status[1]==0
        if i<=g.capacity
            good=isfinite(a.cell_inert[i])&&isfinite(a.cell_inert_low[i])&&Float64(a.cell_inert[i])+a.cell_inert_low[i]>=0
            for k in 1:10
                value=Float64(a.cell_ledger[i,k])+a.cell_ledger_low[i,k]
                good&=isfinite(value)&&(k>=9||value>=0)
            end
            for row in 1:3,col in 1:Int(g.hidden)+10;good&=isfinite(a.weights[i,col,row]);end
            !good && dh_fail!(e.status,64)
        end
        if i<=length(a.grid_inert)
            good=isfinite(a.grid_inert[i])&&isfinite(a.grid_inert_low[i])&&Float64(a.grid_inert[i])+a.grid_inert_low[i]>=0
            for k in 1:2
                enzyme=reinterpret(Float32,a.enzyme_bits[i,k]);bound=a.bound[i,k]
                good&=isfinite(enzyme)&&enzyme>=0&&isfinite(bound)&&bound>=0&&isfinite(feed.fields[i,k])&&feed.fields[i,k]>=0
            end
            for k in 1:12
                value=Float64(a.grid_ledger[i,k])+a.grid_ledger_low[i,k]
                good&=isfinite(value)&&value>=0
            end
            !good && dh_fail!(e.status,64)
        end
    end
end
