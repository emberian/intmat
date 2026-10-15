# Preserved allocation-heavy mechanics oracle from the pinned repair-weave pilot v2.
# Constructor, state layout and validation are unchanged by the scalar-loop patch.
module BondNetworkReference
using LinearAlgebra
using Interior: World, BondNetwork
import Interior
function check_dt(w,p,dt)
    isfinite(dt) && dt>0 || throw(ArgumentError("dt"))
    isfinite(Float32(dt)) && Float32(dt)>0 || throw(ArgumentError("Float32 dt"))
    Interior.validate(w,p)
    Float32(dt)<=stability_limit(w,p) || throw(ArgumentError("stability"))
end
function bond_forces(p::BondNetwork)
    f=zeros(size(p.position))
    for e in eachindex(p.bond_ids)
        p.damage[e]==1 && continue
        i,j=p.bonds[:,e];delta=p.position[:,j]-p.position[:,i];l=norm(delta)
        l>1e-14 || continue
        force=(1-p.damage[e])*p.stiffness[e]*(l-p.rest_length[e])/l .* delta
        f[:,i].+=force;f[:,j].-=force
    end
    f
end
function mechanical_energy(p::BondNetwork)
    elastic=0.
    for e in eachindex(p.bond_ids)
        i,j=p.bonds[:,e];extension=norm(p.position[:,j]-p.position[:,i])-p.rest_length[e]
        elastic+=.5*(1-p.damage[e])*p.stiffness[e]*extension^2
    end
    kinetic=sum(.5*p.mass[i]*sum(abs2,p.velocity[:,i]) for i in eachindex(p.mass))
    (elastic=elastic,kinetic=kinetic,damage=p.damage_energy,damping=p.damping_energy,
     external_work=p.external_work,total=elastic+kinetic+p.damage_energy+p.damping_energy)
end
function break_bonds!(p::BondNetwork)
    for e in eachindex(p.bond_ids)
        i,j=p.bonds[:,e];extension=norm(p.position[:,j]-p.position[:,i])-p.rest_length[e]
        strain=extension/p.rest_length[e];p.peak_strain[e]=max(p.peak_strain[e],strain)
        if p.damage[e]<1 && strain>=p.threshold[e]
            p.damage_energy+=.5*(1-p.damage[e])*p.stiffness[e]*extension^2
            p.damage[e]=1.
        end
    end
    p
end
function stability_limit(w::World,p::BondNetwork)
    # Gershgorin bound on the current central-force tangent, including compressed
    # transverse stiffness; conservative safety margin for velocity Verlet.
    rows=zeros(length(p.mass))
    for e in eachindex(p.bond_ids)
        p.damage[e]==1 && continue
        i,j=p.bonds[:,e];l=norm(p.position[:,j]-p.position[:,i])
        tangent=p.stiffness[e]*(1-p.damage[e])*max(1.,abs(1-p.rest_length[e]/max(l,1e-14)))
        rows[i]+=tangent;rows[j]+=tangent
    end
    bound=maximum(2 .* rows ./ p.mass)
    bound==0 ? Inf : .35/sqrt(bound)
end
function mechanical_grips!(p,t)
    for i in eachindex(p.grip)
        s=p.grip[i];s==0 && continue
        p.position[:,i].=p.reference[:,i]
        p.position[1,i]+=s*p.stretch_rate*t/2
        p.position[2,i]+=s*p.shear_rate*t/2
        p.velocity[:,i].=0
    end
end
function advance!(w::World,p::BondNetwork,dt::Real)
    check_dt(w,p,dt);dt=Float64(dt)
    # Grip displacement is a kinematic external operation; book its exact change
    # in elastic energy separately from integration error and damage loss.
    before=mechanical_energy(p).elastic
    mechanical_grips!(p,p.elapsed+dt)
    p.external_work+=mechanical_energy(p).elastic-before
    break_bonds!(p)
    f=bond_forces(p); decay=exp(-p.damping*dt/2)
    for i in eachindex(p.mass)
        p.grip[i]!=0 && continue
        old=sum(abs2,p.velocity[:,i]);p.velocity[:,i].*=decay
        p.damping_energy+=.5*p.mass[i]*(old-sum(abs2,p.velocity[:,i]))
        p.velocity[:,i].+=dt/2 .* (f[:,i]+p.load[:,i])./p.mass[i]
        dx=dt .* p.velocity[:,i];p.position[:,i].+=dx
        p.external_work+=dot(p.load[:,i],dx)
    end
    break_bonds!(p);f=bond_forces(p)
    for i in eachindex(p.mass)
        p.grip[i]!=0 && continue
        p.velocity[:,i].+=dt/2 .* (f[:,i]+p.load[:,i])./p.mass[i]
        old=sum(abs2,p.velocity[:,i]);p.velocity[:,i].*=decay
        p.damping_energy+=.5*p.mass[i]*(old-sum(abs2,p.velocity[:,i]))
    end
    p.elapsed+=dt;w
end
end
