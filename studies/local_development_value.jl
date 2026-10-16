"""Local developmental accounting; these values are not physical energy stores."""
module LocalDevelopmentValue
@inline usable(q,energy,a,b)=energy+q.efficiency*(q.potentials[1]*a+q.potentials[2]*b)
@inline substrate_coefficient(q)=Float64(q.efficiency)*(Float64(q.potentials[1])+Float64(q.potentials[2]))/2
@inline replacement_coefficient(q)=substrate_coefficient(q)+Float64(q.growth_work)
@inline replacement_value(q,mass,energy,a,b)=usable(q,energy,a,b)+replacement_coefficient(q)*mass
# A relaxed instantaneous capacity: rates, energy-capacity overflow and native
# proportional fuel withdrawal can make the bound unattainable in finite time.
@inline function buildable_capacity(q,mass,energy,a,b)
    lambda=replacement_coefficient(q)
    fuelbound=lambda>0 ? usable(q,energy,a,b)/lambda : oftype(mass,Inf)
    mass+min(2a,2b,fuelbound)
end
end
