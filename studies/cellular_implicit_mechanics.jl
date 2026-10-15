# Experimental opt-in solver. Not included by Interior's default module.
# Study callers may Base.include(Interior, this_file); no cache is persisted.
struct _ImplicitContactEdge
    i::Int
    j::Int
    normal::NTuple{3,Float64}
    stiffness::Float64
    tangential::Float64
end
@inline function _implicit_dot(a,b)
    total=0.
    @inbounds for i in eachindex(a,b);total+=a[i]*b[i];end
    total
end
function _implicit_edges(grid,ids,x,radii,bonds,index,repulsion,stiffness,pairs,h,exact_curvature)
    edges=_ImplicitContactEdge[];diagonal=fill(inv(h),size(x))
    function add(i,j,k,rest)
        delta=_cellular_delta(grid,x,i,j);distance=_cellular_length(delta)
        normal=distance>0 ? delta./distance : (ids[i]<ids[j] ? (1.,0.,0.) : (-1.,0.,0.))
        tangent=exact_curvature&&distance>0 ? k*(1-rest/distance) : 0.
        push!(edges,_ImplicitContactEdge(i,j,normal,k,tangent))
        for d in 1:3
            value=k*normal[d]^2+abs(tangent)*(1-normal[d]^2)
            diagonal[d,i]+=value;diagonal[d,j]+=value
        end
    end
    for (i,j) in pairs
        _cellular_length(_cellular_delta(grid,x,i,j))<radii[i]+radii[j] && add(i,j,repulsion,radii[i]+radii[j])
    end
    for e in eachindex(bonds.endpoints)
        a,b=bonds.endpoints[e];i=get(index,a,0);j=get(index,b,0)
        (i==0 || j==0) && continue
        _cellular_length(_cellular_delta(grid,x,i,j))>bonds.rest_lengths[e] && add(i,j,stiffness,bonds.rest_lengths[e])
    end
    edges,diagonal
end
function _implicit_hessian!(out,v,edges,h)
    @. out=v/h
    @inbounds for edge in edges
        i=edge.i;j=edge.j;n=edge.normal
        projection=0.
        for d in 1:3;projection+=n[d]*(v[d,i]-v[d,j]);end
        for d in 1:3
            value=edge.stiffness*n[d]*projection+edge.tangential*(v[d,i]-v[d,j]-n[d]*projection)
            out[d,i]+=value;out[d,j]-=value
        end
    end
    out
end
function _implicit_cg(gradient,diagonal,edges,h,maxiter)
    step=zeros(size(gradient));r=-gradient;z=r./diagonal;p=copy(z);hp=similar(p)
    rz=_implicit_dot(r,z);initial=_implicit_dot(r,r);iterations=0
    for iteration in 1:maxiter
        initial==0 && break
        _implicit_hessian!(hp,p,edges,h)
        denominator=_implicit_dot(p,hp)
        if !(denominator>0&&isfinite(denominator))
            iterations==0&&copyto!(step,z)
            break
        end
        alpha=rz/denominator
        @. step=step+alpha*p
        @. r=r-alpha*hp
        iterations=iteration
        _implicit_dot(r,r)<=.0025initial && break
        @. z=r/diagonal
        next_rz=_implicit_dot(r,z)
        beta=next_rz/rz;rz=next_rz
        @. p=z+beta*p
    end
    step,iterations
end
function _implicit_state(grid,ids,y,displacement,radii,bonds,index,repulsion,stiffness,h)
    x=y+displacement;nd=_cellular_ndim(grid)
    for i in axes(x,2),d in 1:3
        x[d,i]=d>nd ? 0. : mod(x[d,i],Float64(grid.extent[d]))
    end
    pairs=cellular_pairs(grid,x,radii)
    state=_cellular_energy_forces(grid,ids,x,radii,bonds,index,repulsion,stiffness,Val(true),pairs)
    gradient=displacement/h-state.forces
    objective=state.energy+_implicit_dot(displacement,displacement)/(2h)
    (;x,pairs,state,gradient,objective)
end

"""Study-only backward-Euler passive relaxation of the existing soft potential.

Minimizes E(mod(y+s))+||s||²/(2h), starting at s=0. The search matrix
uses the exact radial pair Hessian plus I/h, truncating CG on nonpositive
curvature. The positive diagonal preconditioner uses absolute transverse
curvature. `exact_curvature=false` retains a normal-block Gauss–Newton study
reference. Armijo tests use true potential and newly discovered contacts at
every trial.
Returns `converged=false` without modifying any input when residual/line search
fails. A caller must reduce h or reject the step, never silently accept it.
"""
function cellular_proximal_relax(grid,ids,y,radii,bonds::CellBonds,h;
    repulsion=10.,stiffness=2.,tolerance=1e-7,max_iterations=80,max_cg=80,exact_curvature=true,max_evaluations=10000)
    _cellular_validate_geometry(grid,y,radii);_cellular_validate_bonds(bonds)
    isfinite(h)&&h>0 && isfinite(tolerance)&&tolerance>0 || throw(ArgumentError("positive h and tolerance required"))
    all(k->isfinite(k)&&k>=0,(repulsion,stiffness)) || throw(ArgumentError("invalid stiffness"))
    length(ids)==length(radii)&&length(unique(ids))==length(ids) || throw(ArgumentError("unique cell IDs required"))
    max_iterations>0&&max_cg>0&&max_evaluations>0 || throw(ArgumentError("positive solver iteration limits required"))
    index=Dict(id=>i for (i,id) in enumerate(ids));s=zeros(size(y))
    current=_implicit_state(grid,ids,y,s,radii,bonds,index,repulsion,stiffness,h)
    initial_energy=current.state.energy
    scale=1+maximum(abs,current.gradient;init=0.)
    evaluations=1;cg_iterations=0;iterations=0;converged=false
    # Stay within a single local proximal chart, including across periodic seams.
    chart_limit=minimum(Float64(grid.extent[d]) for d in 1:_cellular_ndim(grid))/4
    for iteration in 0:max_iterations
        residual=maximum(abs,current.gradient;init=0.)
        if residual<=tolerance*scale
            converged=true;break
        end
        iteration==max_iterations && break
        edges,diagonal=_implicit_edges(grid,ids,current.x,radii,bonds,index,repulsion,stiffness,current.pairs,h,exact_curvature)
        direction,count=_implicit_cg(current.gradient,diagonal,edges,h,max_cg)
        cg_iterations+=count;slope=_implicit_dot(current.gradient,direction)
        slope<0 && isfinite(slope) || break
        alpha=1.;accepted=false
        for _ in 1:35
            evaluations>=max_evaluations&&break
            trial_s=s+alpha*direction
            if maximum(abs,trial_s;init=0.)<=chart_limit
                trial=_implicit_state(grid,ids,y,trial_s,radii,bonds,index,repulsion,stiffness,h)
                evaluations+=1
                armijo=trial.objective<=current.objective+1e-4alpha*slope
                # Near a stationary point, objective differences can round away
                # before the gradient residual meets tolerance. Permit only a
                # roundoff-sized change with a strictly improved residual;
                # retain the exact initial-objective passivity certificate.
                roundoff=8eps(max(1.,abs(current.objective)))
                refined=trial.objective<=current.objective+roundoff &&
                    _implicit_dot(trial.gradient,trial.gradient)<.99*_implicit_dot(current.gradient,current.gradient)
                if isfinite(trial.objective) && trial.objective<=initial_energy && (armijo||refined)
                    s=trial_s;current=trial;accepted=true;break
                end
            end
            alpha/=2
        end
        accepted || break
        iterations+=1
    end
    residual=maximum(abs,current.gradient;init=0.)
    # Starting from s=0 and monotone true-objective descent certifies passivity.
    passive_loss=initial_energy-current.state.energy
    (;positions=current.x,converged,residual,iterations,cg_iterations,evaluations,
      potential_before=initial_energy,potential_after=current.state.energy,
      relaxation_energy_loss=passive_loss,proximal_objective=current.objective,
      displacement=s)
end

function _study_funded_drift(grid,nd,x,velocity,funds,h,motor_cost,potential)
    velocity=copy(velocity);start_energy=potential(x);n=size(x,2)
    y=similar(x);weights=zeros(n);spend=zeros(n)
    for _ in 1:80
        @. y=x+h*velocity
        for i in 1:n,d in 1:nd;y[d,i]=mod(y[d,i],Float64(grid.extent[d]));end
        energy=potential(y)
        for i in 1:n;weights[i]=sum(abs2,view(velocity,:,i));end
        weight_sum=sum(weights);cost=max(0.,energy-start_energy)
        for i in 1:n;spend[i]=motor_cost*h*weights[i]+(weight_sum>0 ? cost*weights[i]/weight_sum : 0.);end
        if funds===nothing||all(i->spend[i]<=funds[i],eachindex(spend))
            return (;positions=y,weights,spend,drag=motor_cost*h*sum(weights),potential_work=energy-start_energy)
        end
        for i in 1:n
            spend[i]>funds[i]&&(velocity[:,i].*=.9sqrt(funds[i]/spend[i]))
        end
    end
    throw(ArgumentError("split motor funding failed; no state committed"))
end

_study_potential_pairs(relaxer,grid,x,radii)=nothing

"""Opt-in study integrator: funded motor drift followed by implicit relaxation.

Uses the existing soft potential, motor-work sharing rule and consent-based
bond refresh. It changes discretization: topology is refreshed at accepted
macrostep endpoints. This is a new numerical experiment, not an exact replay
replacement for `cellular_mechanics!`. Passive solver failures halve h and retry;
no partially accepted step is published. `symmetric_motor=true` uses two half motor drifts around the passive step,
with the same debit rule applied to each half. Motor displacement remains bounded by
0.2 minimum radius, independently of stiffness. Passive displacement is bounded
by the local periodic proximal chart, not a hard-sphere constraint.
"""
function cellular_implicit_mechanics!(grid,ids,positions,radii,bonds::CellBonds,dt::Real;
    adhesion=ones(length(ids)),motility=zeros(3,length(ids)),repulsion=10.,stiffness=2.,
    break_strain=.75,max_degree::Integer=6,adhesion_threshold=.5,bond_range=0.,
    energy=nothing,motor_cost=1.,tolerance=1e-7,max_iterations=80,max_cg=80,
    max_substeps=10000,max_step=Inf,exact_curvature=true,max_evaluations=10000,
    relaxer=cellular_proximal_relax,symmetric_motor=false)
    _cellular_validate_geometry(grid,positions,radii);_cellular_validate_bonds(bonds)
    n=length(ids);nd=_cellular_ndim(grid)
    length(radii)==length(adhesion)==n && size(motility)==(3,n) || throw(DimensionMismatch("cell action dimensions differ"))
    length(unique(ids))==n || throw(ArgumentError("unique cell IDs required"))
    all(a->isfinite(a)&&0<=a<=1,adhesion)&&all(isfinite,motility) || throw(ArgumentError("invalid adhesion or motility"))
    isfinite(dt)&&dt>=0&&all(k->isfinite(k)&&k>=0,(repulsion,stiffness,break_strain,bond_range))&&
        0<adhesion_threshold<=1&&max_degree>=0&&max_substeps>0 || throw(ArgumentError("invalid mechanics parameters"))
    isfinite(motor_cost)&&motor_cost>0&&max_step>0 || throw(ArgumentError("invalid motor cost or max step"))
    if energy!==nothing
        length(energy)==n&&all(e->isfinite(e)&&e>=0,energy) || throw(ArgumentError("invalid motor energy budget"))
    end
    max_evaluations>0 || throw(ArgumentError("positive objective work budget required"))
    # Own every mutable candidate until the entire requested interval succeeds.
    x=Float64.(positions);working_bonds=deepcopy(bonds);funds=energy===nothing ? nothing : Float64.(energy)
    for i in 1:n,d in 1:3;x[d,i]=d>nd ? 0. : mod(x[d,i],Float64(grid.extent[d]));end
    index=Dict(id=>i for (i,id) in enumerate(ids))
    potential(z)=_cellular_energy_forces(grid,ids,z,radii,working_bonds,index,repulsion,stiffness,
        Val(false),_study_potential_pairs(relaxer,grid,z,radii)).energy
    potential_before=potential(x)
    refresh=_cellular_refresh_bonds!(grid,ids,x,radii,working_bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness)
    formed=refresh.formed;broken=refresh.broken;topology_energy_loss=refresh.loss
    spent_motor_work=zeros(n);active_distance=zeros(n);motor_drag_loss=0.;motor_potential_work=0.;relaxation_energy_loss=0.
    motors=Float64.(motility);nd==2&&(motors[3,:].=0.)
    remaining=Float64(dt);substeps=0;nonlinear_iterations=0;cg_iterations=0;evaluations=0;retries=0;max_residual=0.
    while remaining>0&&n>0
        substeps<max_substeps || throw(ArgumentError("implicit mechanics exceeded macrostep limit"))
        velocity=copy(motors)
        if funds!==nothing
            for i in 1:n
                drag=motor_cost*remaining*sum(abs2,view(velocity,:,i))
                drag>funds[i]&&(velocity[:,i].*=sqrt(funds[i]/drag))
            end
        end
        speed=maximum(sqrt(sum(abs2,view(velocity,:,i))) for i in 1:n)
        h=min(remaining,max_step,speed>0 ? .2minimum(radii)/speed : Inf)
        # Absorb a rounding-sized remainder into the final actual step. Without
        # this, e.g. .1 - 4*.025 can leave ~7e-18: a proximal force residual
        # cannot be resolved by coordinate changes smaller than an ulp.
        # No physical interval is skipped and every energy/work gate still runs.
        remaining-h<=8eps(Float64(dt))&&(h=remaining)
        start_energy=potential(x);accepted=false
        for attempt in 1:50
            motor_h=symmetric_motor ? h/2 : h
            funded_velocity=copy(velocity);motor_x=similar(x);spend=zeros(n);weights=zeros(n);motor_energy=start_energy;funded=false
            for _ in 1:80
                @. motor_x=x+motor_h*funded_velocity
                for i in 1:n,d in 1:nd;motor_x[d,i]=mod(motor_x[d,i],Float64(grid.extent[d]));end
                motor_energy=potential(motor_x)
                for i in 1:n;weights[i]=sum(abs2,view(funded_velocity,:,i));end
                weight_sum=sum(weights);cost=max(0.,motor_energy-start_energy)
                for i in 1:n;spend[i]=motor_cost*motor_h*weights[i]+(weight_sum>0 ? cost*weights[i]/weight_sum : 0.);end
                if funds===nothing||all(i->spend[i]<=funds[i],eachindex(spend))
                    funded=true;break
                end
                for i in 1:n
                    spend[i]>funds[i]&&(funded_velocity[:,i].*=.9sqrt(funds[i]/spend[i]))
                end
            end
            funded || throw(ArgumentError("implicit motor funding failed"))
            result=relaxer(grid,ids,motor_x,radii,working_bonds,h;
                repulsion,stiffness,tolerance,max_iterations,max_cg,exact_curvature,
                max_evaluations=max_evaluations-evaluations)
            evaluations+=result.evaluations;nonlinear_iterations+=result.iterations;cg_iterations+=result.cg_iterations
            if result.converged&&result.potential_after<=motor_energy
                second=if symmetric_motor
                    remainder=funds===nothing ? nothing : funds-spend
                    _study_funded_drift(grid,nd,result.positions,velocity,remainder,h/2,motor_cost,potential)
                else
                    nothing
                end
                max_residual=max(max_residual,result.residual)
                motor_potential_work+=motor_energy-start_energy+(second===nothing ? 0. : second.potential_work)
                relaxation_energy_loss+=motor_energy-result.potential_after
                motor_drag_loss+=motor_cost*motor_h*sum(abs2,funded_velocity)+(second===nothing ? 0. : second.drag)
                for i in 1:n
                    spent_motor_work[i]+=spend[i]+(second===nothing ? 0. : second.spend[i])
                    active_distance[i]+=motor_h*sqrt(weights[i])+(second===nothing ? 0. : h/2*sqrt(second.weights[i]))
                    if funds!==nothing
                        funds[i]-=spend[i]
                        second!==nothing&&(funds[i]-=second.spend[i])
                    end
                end
                x=second===nothing ? result.positions : second.positions;accepted=true;break
            end
            evaluations<max_evaluations || throw(ArgumentError("implicit objective work budget exhausted; no state committed"))
            retries+=1;h/=2
        end
        accepted&&h>0&&remaining-h<remaining || throw(ArgumentError("implicit mechanical solve failed; no state committed"))
        remaining=max(0.,remaining-h);substeps+=1
        refresh=_cellular_refresh_bonds!(grid,ids,x,radii,working_bonds,index,adhesion,adhesion_threshold,break_strain,max_degree,bond_range,stiffness)
        formed+=refresh.formed;broken+=refresh.broken;topology_energy_loss+=refresh.loss
    end
    converted=eltype(positions).(x);potential_after=potential(converted)
    energy_residual=potential_after-potential_before-motor_potential_work+relaxation_energy_loss+topology_energy_loss
    positions.=converted;bonds.endpoints=working_bonds.endpoints;bonds.rest_lengths=working_bonds.rest_lengths
    energy!==nothing&&(energy.=funds)
    (;potential_before,potential_after,motor_potential_work,relaxation_energy_loss,topology_energy_loss,
      energy_residual,formed,broken,bonds=length(bonds.endpoints),substeps,spent_motor_work,motor_drag_loss,
      active_distance,nonlinear_iterations,cg_iterations,evaluations,retries,max_residual)
end
