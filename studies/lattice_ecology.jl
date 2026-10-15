"""Fixed-voxel ecology: funded local transfers and inherited finite-locus rules.
A separate constitutive study; it neither imports nor migrates running arenas.
"""
module LatticeEcologyStudy
using Interior, Random, Statistics, SHA, Serialization
using Interior: KernelAbstractions
using KernelAbstractions: @kernel, @index
import Interior: advance!,validate,stability_limit,process_ports,record_process_fields,owned_material_mass
export LatticeEcology,lattice_world,lattice_parameters,lattice_statistics,
    lattice_checkpoint,restore_lattice,lattice_genotypes,lattice_damage!,ENV_CHANNELS
const MODEL_SOURCE=read(@__FILE__,String)
const MODEL_SHA=bytes2hex(sha256(MODEL_SOURCE))
const ENV_CHANNELS=(:fuel_a,:fuel_b,:precursor,:matrix,:spent,:reserve_a,:reserve_b,:occupancy_observation)
const HIDDEN=4
const INPUTS=14
const OUTPUTS=11
# No physical atomics: each local stage owns one voxel; each contact stage owns
# a disjoint periodic nearest-neighbor matching. Hashes are counter-based events.
@inline function mix32(x::UInt32)
    x=(x⊻(x>>16))*0x7feb352d;x=(x⊻(x>>15))*0x846ca68b;x⊻(x>>16)
end
@inline random01(seed,site,tick,stream)=Float32((mix32(seed⊻UInt32(site)*0x9e3779b9⊻UInt32(tick)*0x85ebca6b⊻UInt32(stream)*0xc2b2ae35)>>8)+UInt32(1))*5.9604645f-8
@inline squash(x)=x/(1f0+abs(x))
@inline gate(x)=inv(1f0+exp(-clamp(x,-20f0,20f0)))
@inline function heat_add!(heat,i,x)
    # Compensated Float32 owner-local accumulation; no Float64 device requirement.
    old=heat[i,1];next=old+x
    correction=abs(old)>=abs(x) ? (old-next)+x : (x-next)+old
    heat[i,1]=next;heat[i,2]+=correction
    nothing
end
@inline coords(site,dims)=(rem(site-1,dims[1])+1,rem(div(site-1,dims[1]),dims[2])+1,div(site-1,dims[1]*dims[2])+1)
@inline index(i,j,k,dims)=i+dims[1]*((j-1)+dims[2]*(k-1))
@inline function neighbor(site,dims,axis,direction)
    c=coords(site,dims)
    v=ntuple(d->d==axis ? mod(c[d]+direction-1,dims[d])+1 : c[d],Val(3))
    index(v...,dims)
end
function lattice_parameters(;feed=.04,fuel_diffusion=.12,precursor_diffusion=.08,
        matrix_barrier=8.,matrix_decay=.003,uptake=.7,efficiency=.65,
        growth=.9,growth_work=.35,construction=.06,construction_work=.4,
        erosion=.025,erosion_work=.12,basal_work=.025,neural_work=.005,
        move_rate=.8,move_work=.04,division_rate=1.5,division_work=.06,
        minimum_mass=.15,division_mass=.65,capacity=1.2,energy_capacity=2.,
        lifespan=24.,response_time=.7,mutation_probability=.04,
        taxis=2.,recurrence=1.)
    q=(;feed,fuel_diffusion,precursor_diffusion,matrix_barrier,matrix_decay,uptake,efficiency,
       growth,growth_work,construction,construction_work,erosion,erosion_work,basal_work,neural_work,
       move_rate,move_work,division_rate,division_work,minimum_mass,division_mass,capacity,
       energy_capacity,lifespan,response_time,mutation_probability,taxis,recurrence)
    all(x->isfinite(x)&&x>=0,values(q)) || throw(ArgumentError("finite nonnegative lattice parameters required"))
    0<efficiency<=1 && 0<=mutation_probability<=1 && 0<=recurrence<=1 || throw(ArgumentError("invalid probability/gain"))
    all(>(0),(growth_work,construction_work,erosion_work,move_work,division_work,minimum_mass,energy_capacity,lifespan,response_time)) || throw(ArgumentError("positive work/scales required"))
    2minimum_mass<=division_mass<=capacity || throw(ArgumentError("division must fund two minimum cells within capacity"))
    p=map(Float32,q)
    all(isfinite,values(p)) && all(>(0),(p.growth_work,p.construction_work,p.erosion_work,p.move_work,p.division_work,p.minimum_mass,p.energy_capacity,p.lifespan,p.response_time,p.efficiency)) || throw(ArgumentError("positive scales must remain finite and nonzero in Float32"))
    2p.minimum_mass<=p.division_mass<=p.capacity || throw(ArgumentError("Float32 division/capacity mismatch"))
    isfinite(p.energy_capacity*p.capacity)||throw(ArgumentError("body-scaled battery capacity overflows Float32"))
    p
end
lattice_world(grid;seed=11,backend=CPU())=World(grid;seed,backend,fields=(),scratch_fields=())
function lattice_genotypes(;loci=12,seed=911,strength=.065)
    1<=loci<=14 || throw(ArgumentError("shared lookup supports 1–14 binary loci"))
    isfinite(strength)&&strength>=0&&isfinite(Float32(strength))||throw(ArgumentError("finite nonnegative allele strength required"))
    rng=MersenneTwister(seed);G=1<<loci
    base=.16f0.*randn(rng,Float32,INPUTS,OUTPUTS)
    for h in 1:HIDDEN;base[h,h]+=.5f0;end
    # Disclosed viable founder priors, not evolved anatomy or learned behavior.
    for (a,bias) in enumerate((0.,1.6,1.3,-1.,-1.5,-.7,1.8));base[end,HIDDEN+a]=Float32(bias);end
    basis=Float32(strength).*randn(rng,Float32,loci,INPUTS,OUTPUTS)
    bank=Array{Float32}(undef,G,INPUTS,OUTPUTS)
    for output in 1:OUTPUTS,input in 1:INPUTS,code in 0:G-1
        value=base[input,output]
        for bit in 0:loci-1;value+=(((code>>bit)&1)==1 ? 1f0 : -1f0)*basis[bit+1,input,output];end
        bank[code+1,input,output]=value
    end
    all(isfinite,bank)||throw(ArgumentError("allele table overflows Float32"))
    bank
end
mutable struct LatticeEcology{A,I,J,U,W,Q}
    environment::A
    scratch::A
    cells::A
    hidden::A
    actions::A
    heat::A
    work::A
    genotype::I
    lineage::I
    generation::I
    events::J
    ids::U
    parents::U
    weights::W
    parameters::Q
    dims::NTuple{3,Int}
    spacing::NTuple{3,Float32}
    axes::Vector{Int}
    loci::Int
    seed::UInt32
    tick::Int
end
function LatticeEcology(w;parameters=lattice_parameters(),loci=12,occupancy=.1,
        founder_genotypes=8,initial_mass=.32,precursor=.6,fuel=.25,reserve=8.,niches=false,
        genotype_seed=911,initial_codes=nothing)
    0<=occupancy<=1 && 0<founder_genotypes<=(1<<loci) || throw(ArgumentError("invalid occupancy/founder genotype count"))
    all(x->isfinite(x)&&x>=0&&isfinite(Float32(x)),(precursor,fuel,reserve,initial_mass)) || throw(ArgumentError("invalid inventory"))
    isfinite(Float32(fuel)+(niches ? 2f0 : 1f0)*Float32(reserve))||throw(ArgumentError("combined fuel/reserve exceeds Float32 range"))
    q=lattice_parameters(;parameters...);q.minimum_mass<=Float32(initial_mass)<=q.capacity || throw(ArgumentError("founder mass outside capacity"))
    0<=w.seed<=typemax(UInt32)||throw(ArgumentError("counter-based lattice seed must fit UInt32"))
    dims=w.grid.size;all(n->n==1||iseven(n),dims) || throw(ArgumentError("matching requires even periodic sides (or singleton axes)"))
    V=prod(dims);V<typemax(Int32) || throw(ArgumentError("voxel indexing exceeds Int32"))
    rng=MersenneTwister(w.seed);host=zeros(Float32,V,8);host[:,1:2].=Float32(fuel);host[:,3].=Float32(precursor)
    host[:,6:7].=Float32(reserve)
    if niches
        isfinite(2Float32(reserve))||throw(ArgumentError("niche reserve overflows Float32"))
        for site in 1:V
            x=coords(site,dims)[1];left=x<=dims[1]÷2
            host[site,6]=left ? 2Float32(reserve) : 0;host[site,7]=left ? 0 : 2Float32(reserve)
        end
    end
    cells=zeros(Float32,V,3);g=zeros(Int32,V);lin=zeros(Int32,V);ids=zeros(UInt64,V)
    codes=initial_codes===nothing ? rand(rng,1:(1<<loci),founder_genotypes) : collect(initial_codes)
    !isempty(codes)&&all(x->x isa Integer&&1<=x<=(1<<loci),codes)||throw(ArgumentError("invalid founder code"))
    for site in 1:V
        if rand(rng)<occupancy
            g[site]=rand(rng,codes);lin[site]=site;ids[site]=site
            cells[site,1]=Float32(initial_mass);cells[site,2]=min(Float32(initial_mass),q.energy_capacity*Float32(initial_mass))
        end
    end
    z(T,ds...)=KernelAbstractions.zeros(w.backend,T,ds...)
    env=z(Float32,V,8);copyto!(env,host);c=z(Float32,V,3);copyto!(c,cells)
    gi=z(Int32,V);copyto!(gi,g);li=z(Int32,V);copyto!(li,lin);id=z(UInt64,V);copyto!(id,ids)
    weights=lattice_genotypes(;loci,seed=genotype_seed);bank=z(Float32,size(weights)...);copyto!(bank,weights)
    p=LatticeEcology(env,z(Float32,V,8),c,z(Float32,V,HIDDEN),z(Float32,V,3),z(Float32,V,2),z(Float32,V,5),
        gi,li,z(Int32,V),z(Int32,V,5),id,z(UInt64,V),bank,q,dims,w.grid.spacing,
        [d for d in 1:3 if dims[d]>1],loci,UInt32(w.seed),0)
    validate(w,p);p
end
function validate(w,p::LatticeEcology)
    w.seed==Int(p.seed)||throw(ArgumentError("lattice/world seed mismatch"))
    w.grid.size==p.dims&&w.grid.spacing==p.spacing || throw(ArgumentError("lattice geometry mismatch"))
    V=prod(p.dims)
    size(p.environment)==size(p.scratch)==(V,8)&&size(p.cells)==(V,3)&&size(p.hidden)==(V,HIDDEN)&&size(p.actions)==(V,3) || throw(DimensionMismatch("lattice arrays"))
    size(p.weights)==(1<<p.loci,INPUTS,OUTPUTS) || throw(DimensionMismatch("genotype lookup"))
    size(p.heat)==(V,2)&&size(p.work)==size(p.events)==(V,5)&&all(a->length(a)==V,(p.genotype,p.lineage,p.generation,p.ids,p.parents))||throw(DimensionMismatch("lattice identity/ledger arrays"))
    arrays=(p.environment,p.scratch,p.cells,p.hidden,p.actions,p.heat,p.work,p.genotype,p.lineage,p.generation,p.events,p.ids,p.parents,p.weights)
    all(a->KernelAbstractions.get_backend(a)==w.backend,arrays)||throw(ArgumentError("lattice backend mismatch"))
    for i in 1:length(arrays),j in 1:i-1;Base.mightalias(arrays[i],arrays[j])&&throw(ArgumentError("lattice owner arrays must not alias"));end
    lattice_parameters(;p.parameters...)
    0<=p.tick<typemax(UInt32)||throw(ArgumentError("lattice clock outside representable range"))
    nothing
end
function stability_limit(w,p::LatticeEcology)
    inverse=sum((inv(Float64(w.grid.spacing[d])^2) for d in p.axes);init=0.)
    min(1.,inverse>0 ? .9/(2inverse*max(p.parameters.fuel_diffusion,p.parameters.precursor_diffusion)) : 1.,
        .9/max(2p.parameters.move_rate,p.parameters.division_rate,eps(Float32)))
end
@kernel function diffuse_kernel!(out,env,cells,dims,spacing,q,dt)
    site=@index(Global,Linear)
    @inbounds begin
        for channel in 1:7
            value=env[site,channel];flux=0f0
            diffusion=channel<=2 ? q.fuel_diffusion : channel==3 ? q.precursor_diffusion : 0f0
            if diffusion>0
                for d in 1:3
                    if dims[d]>1
                        for direction in (-1,1)
                            other=neighbor(site,dims,d,direction)
                            face=inv(1f0+q.matrix_barrier*.5f0*(env[site,4]+env[other,4]))
                            weight=((dt*diffusion)/spacing[d])/spacing[d]
                            flux+=weight*face*(env[other,channel]-value)
                        end
                    end
                end
            end
            out[site,channel]=value+flux
        end
        # Local finite high-potential reservoirs feed dissolved fuel; no source
        # ledger is needed because both sides remain in the owned inventory.
        for channel in 1:2
            before=out[site,channel+5]
            remaining=max(0f0,before-min(before,dt*q.feed))
            amount=before-remaining # actual representable debit; sub-ULP requests stall
            out[site,channel]+=amount;out[site,channel+5]=remaining
        end
        decay=min(out[site,4],dt*q.matrix_decay*out[site,4])
        out[site,4]-=decay;out[site,3]+=decay
        out[site,8]=cells[site,1]>0 ? 1f0 : 0f0
    end
end
# Charge at least the representable positive requested work, never overdrawing
# the battery. Upward rounding of a debit dissipates an explicit work quantum.
@inline function charge_work(energy,cost)
    if !(0f0<cost<=energy);return (energy,0f0);end
    after=max(0f0,energy-cost);spent=energy-after
    if spent<cost&&after>0f0
        after=prevfloat(after);spent=energy-after
    end
    (after,spent)
end
@inline function funded_flux(request,energy,unit_work)
    amount=max(0f0,min(request,energy/unit_work));cost=amount*unit_work
    if cost>energy
        amount=max(0f0,prevfloat(amount));cost=amount*unit_work
    end
    # A nonrepresentable positive price cannot buy a free physical transfer.
    if !(amount>0f0&&0f0<cost<=energy);return (0f0,energy,0f0);end
    after,spent=charge_work(energy,cost)
    (amount,after,spent)
end
@inline function network_output(bank,code,inputs,output)
    value=0f0
    for col in 1:INPUTS;value=muladd(bank[code,col,output],inputs[col],value);end
    value
end
@kernel function metabolize_kernel!(env,old,cells,hidden,actions,genotype,lineage,generation,ids,parents,events,heat,work,bank,dims,q,dt)
    site=@index(Global,Linear)
    @inbounds begin
        for channel in 1:8;env[site,channel]=old[site,channel];end
        code=genotype[site]
        if code>0
            b=cells[site,1];energy=cells[site,2];age=cells[site,3]
            occupied=0f0;neighbors=0f0
            for d in 1:3
                if dims[d]>1
                    for direction in (-1,1);occupied+=old[neighbor(site,dims,d,direction),8];neighbors+=1f0;end
                end
            end
            density=neighbors>0 ? occupied/neighbors : 0f0
            # Explicit immutable tuple: a closure over b/energy/age (mutated
            # later below) boxes them as Any and destroys device inference.
            inputs=(q.recurrence*hidden[site,1],q.recurrence*hidden[site,2],
                q.recurrence*hidden[site,3],q.recurrence*hidden[site,4],
                squash(old[site,1]),squash(old[site,2]),squash(old[site,3]),squash(old[site,4]),
                squash(b),squash(energy),age/q.lifespan,density,squash(old[site,1]-old[site,2]),1f0)
            targets=ntuple(j->tanh(network_output(bank,code,inputs,j)),Val(HIDDEN))
            outputs=ntuple(j->gate(network_output(bank,code,inputs,HIDDEN+j)),Val(7))
            preference=outputs[1];feeding=outputs[2]
            for channel in 1:2
                fraction=channel==1 ? preference : 1f0-preference
                potential=channel==1 ? 4f0 : 3f0
                stock=env[site,channel];room=max(0f0,q.energy_capacity*b-energy)
                amount=min(stock,dt*q.uptake*b*feeding*fraction*stock/(.15f0+stock),room/(q.efficiency*potential))
                delivered=q.efficiency*potential*amount
                env[site,channel]-=amount;env[site,5]+=amount;energy+=delivered
                heat_add!(heat,site,potential*amount-delivered)
            end
            requested=dt*(q.basal_work+q.neural_work)*b
            request=min(energy,requested);energy,paid=charge_work(energy,request);heat_add!(heat,site,paid);work[site,1]+=paid
            funding=requested>0 ? min(1f0,paid/requested) : 1f0
            alpha=-expm1(-dt/q.response_time)*funding
            for j in 1:HIDDEN;hidden[site,j]+=alpha*(targets[j]-hidden[site,j]);end
            request=min(env[site,3],max(0f0,q.capacity-b-env[site,4]),dt*q.growth*b*outputs[3]*env[site,3]/(.2f0+env[site,3]))
            grown,energy,spent=funded_flux(request,energy,q.growth_work)
            env[site,3]-=grown;b+=grown;heat_add!(heat,site,spent);work[site,2]+=spent
            request=min(max(0f0,b-q.minimum_mass),dt*q.construction*b*outputs[4])
            built,energy,spent=funded_flux(request,energy,q.construction_work)
            b-=built;env[site,4]+=built;heat_add!(heat,site,spent);work[site,3]+=spent
            request=min(env[site,4],dt*q.erosion*b*outputs[5])
            removed,energy,spent=funded_flux(request,energy,q.erosion_work)
            env[site,4]-=removed;env[site,3]+=removed;heat_add!(heat,site,spent);work[site,3]+=spent
            excess=max(0f0,energy-q.energy_capacity*b);energy-=excess;heat_add!(heat,site,excess)
            age+=dt
            if energy<=0 || age>=q.lifespan
                env[site,3]+=b;heat_add!(heat,site,energy);b=0;energy=0;age=0
                genotype[site]=0;lineage[site]=0;generation[site]=0;ids[site]=0;parents[site]=0;events[site,2]+=1
                for j in 1:HIDDEN;hidden[site,j]=0;end
                for j in 1:3;actions[site,j]=0;end
            else
                actions[site,1]=preference;actions[site,2]=outputs[6];actions[site,3]=outputs[7]
            end
            cells[site,1]=b;cells[site,2]=energy;cells[site,3]=age
        else
            for j in 1:3;actions[site,j]=0;end
        end
    end
end
@kernel function matching_kernel!(cells,hidden,actions,env,genotype,lineage,generation,ids,parents,events,heat,work,dims,q,axis,parity,seed,tick,loci,dt)
    site=@index(Global,Linear)
    @inbounds begin
        coordinate=coords(site,dims)[axis]
        if rem(coordinate,2)==parity
            other=neighbor(site,dims,axis,1)
            a=genotype[site]>0;b=genotype[other]>0
            if a!=b
                donor=a ? site : other;target=a ? other : site
                mass=cells[donor,1];energy=cells[donor,2]
                half=mass*.5f0;gap=max(0f0,q.capacity-env[target,4])
                birth_possible=mass>=q.division_mass&&half>=q.minimum_mass&&half<=gap&&energy>q.division_work
                wants_birth=random01(seed,donor,tick,1)<dt*q.division_rate*actions[donor,3]
                if birth_possible&&wants_birth
                    energy,spent=charge_work(energy,q.division_work);heat_add!(heat,donor,spent);work[donor,5]+=spent
                    cells[donor,1]=mass-half;cells[target,1]=half
                    cells[donor,2]=energy*.5f0;cells[target,2]=energy-cells[donor,2];cells[target,3]=0
                    code=genotype[donor]
                    if random01(seed,target,tick,2)<q.mutation_probability
                        bit=min(loci-1,Int(floor(random01(seed,target,tick,3)*loci)))
                        code=Int32((Int(code)-1)⊻(1<<bit))+1;events[target,4]+=1
                    end
                    genotype[target]=code;lineage[target]=lineage[donor];generation[target]=generation[donor]+1
                    parents[target]=ids[donor];ids[target]=UInt64(tick)*UInt64(length(genotype))+UInt64(target)
                    for j in 1:HIDDEN;hidden[target,j]=0;end
                    for j in 1:3;actions[target,j]=0;end
                    events[target,1]+=1
                else
                    wants_birth && !birth_possible && (events[donor,5]+=1)
                    preference=actions[donor,1]
                    local_food=preference*env[donor,1]+(1f0-preference)*env[donor,2]
                    target_food=preference*env[target,1]+(1f0-preference)*env[target,2]
                    drive=2f0*gate(q.taxis*(target_food-local_food))
                    cost=q.move_work*mass
                    if mass<=gap&&0f0<cost<energy&&random01(seed,donor,tick,4)<dt*q.move_rate*actions[donor,2]*drive
                        energy,spent=charge_work(energy,cost);heat_add!(heat,donor,spent);work[donor,4]+=spent
                        for j in 1:3;cells[target,j]=cells[donor,j];cells[donor,j]=0;end
                        cells[target,2]=energy
                        for j in 1:HIDDEN;hidden[target,j]=hidden[donor,j];hidden[donor,j]=0;end
                        for j in 1:3;actions[target,j]=actions[donor,j];actions[donor,j]=0;end
                        genotype[target]=genotype[donor];genotype[donor]=0
                        lineage[target]=lineage[donor];lineage[donor]=0;generation[target]=generation[donor];generation[donor]=0
                        ids[target]=ids[donor];ids[donor]=0;parents[target]=parents[donor];parents[donor]=0;events[target,3]+=1
                    end
                end
            end
        end
    end
end
function advance!(w::World,p::LatticeEcology,dt::Real)
    Interior.check_dt(w,p,dt);p.tick<typemax(UInt32)-1||throw(ArgumentError("event clock exhausted"))
    dt=Float32(dt);backend=w.backend;V=prod(p.dims);q=p.parameters
    diffuse_kernel!(backend,128)(p.scratch,p.environment,p.cells,p.dims,p.spacing,q,dt;ndrange=V)
    metabolize_kernel!(backend,128)(p.environment,p.scratch,p.cells,p.hidden,p.actions,p.genotype,p.lineage,p.generation,p.ids,p.parents,p.events,p.heat,p.work,p.weights,p.dims,q,dt;ndrange=V)
    p.tick+=1
    if !isempty(p.axes)
        axis=p.axes[mod(div(p.tick-1,2),length(p.axes))+1];parity=mod(p.tick-1,2)
        matching_kernel!(backend,128)(p.cells,p.hidden,p.actions,p.environment,p.genotype,p.lineage,p.generation,p.ids,p.parents,p.events,p.heat,p.work,p.dims,q,axis,parity,p.seed,p.tick,p.loci,dt;ndrange=V)
    end
    w
end
@kernel function damage_kernel!(env,cells,hidden,actions,genotype,lineage,generation,ids,parents,events,heat,dims,center,radius)
    site=@index(Global,Linear)
    @inbounds begin
        coordinate=coords(site,dims);distance=0f0
        for d in 1:3
            if dims[d]>1
                delta=Float32(coordinate[d])-.5f0-center[d]
                delta-=floor(delta/Float32(dims[d])+.5f0)*Float32(dims[d])
                distance+=delta*delta
            end
        end
        if distance<radius*radius
            env[site,3]+=cells[site,1]+env[site,4];env[site,4]=0
            heat_add!(heat,site,cells[site,2])
            genotype[site]>0 && (events[site,2]+=1)
            for j in 1:3;cells[site,j]=0;actions[site,j]=0;end
            for j in 1:HIDDEN;hidden[site,j]=0;end
            genotype[site]=0;lineage[site]=0;generation[site]=0;ids[site]=0;parents[site]=0
        end
    end
end
"""Imposed injury: structural body/matrix becomes precursor; stored work becomes heat.
Fuel and finite reserve remain. This is a declared external perturbation, not a
spontaneously evolved event or a free reset of material/energy.
"""
function lattice_damage!(w,p;center=ntuple(d->Float32(p.dims[d])/2,3),radius=3.)
    validate(w,p);isfinite(radius)&&radius>0||throw(ArgumentError("positive finite injury radius required"))
    damage_kernel!(w.backend,128)(p.environment,p.cells,p.hidden,p.actions,p.genotype,p.lineage,p.generation,p.ids,p.parents,p.events,p.heat,p.dims,Float32.(center),Float32(radius);ndrange=prod(p.dims))
    w
end
process_ports(::LatticeEcology)=ProcessPorts(owned=(:environment,:cells,:hidden,:genotype,:lineage,:generation,:ids,:parents,:heat,:work,:events,:weights),
    drives="finite local fuel reservoirs with potentials 4 and 3; structural precursor has zero chemical potential",
    memory="immobile paid matrix and inherited binary genotype; recurrent state moves with cells and resets at birth",
    law="positivity-bounded local metabolism, symmetric finite-volume transport, paid growth/remodeling and disjoint nearest-neighbor birth/hopping")
record_process_fields(p::LatticeEcology)=(environment=p.environment,cells=p.cells,hidden=p.hidden,genotype=p.genotype,lineage=p.lineage,generation=p.generation,ids=p.ids,parents=p.parents,heat=p.heat,work=p.work,events=p.events,weights=p.weights,channels=ENV_CHANNELS,parameter_names=Tuple(keys(p.parameters)),parameters=Tuple(values(p.parameters)),model_source_sha=MODEL_SHA,tick=p.tick,loci=p.loci)
owned_material_mass(p::LatticeEcology)=sum(Float64,Array(p.environment)[:,1:7])+sum(Float64,Array(p.cells)[:,1])
function lattice_statistics(w,p::LatticeEcology)
    e=Array(p.environment);c=Array(p.cells);g=Array(p.genotype);h=Array(p.hidden);events=Array(p.events);work=Array(p.work)
    live=findall(>(0),g);n=length(live);heat=sum(Float64,Array(p.heat));popcounts=[count(==(v),g[live]) for v in unique(g[live])]
    diversity=n>0 ? exp(-sum(x->(x/n)*log(x/n),popcounts)) : 0.
    mass=sum(Float64,e[:,1:7])+sum(Float64,c[:,1]);potential=4(sum(Float64,e[:,1])+sum(Float64,e[:,6]))+3(sum(Float64,e[:,2])+sum(Float64,e[:,7]))+sum(Float64,c[:,2])
    arrays=(e,c,h,Array(p.actions),Array(p.heat),work)
    (time=w.time,tick=p.tick,population=n,births=sum(events[:,1]),deaths=sum(events[:,2]),hops=sum(events[:,3]),mutations=sum(events[:,4]),blocked_births=sum(events[:,5]),
     material=mass,potential=potential,heat=heat,energy_accounted=potential+heat,
     precursor=sum(Float64,e[:,3]),biomass=sum(Float64,c[:,1]),matrix=sum(Float64,e[:,4]),spent=sum(Float64,e[:,5]),
     fuel=[sum(Float64,e[:,i]) for i in 1:2],reserve=[sum(Float64,e[:,i]) for i in 6:7],
     work=[sum(Float64,work[:,j]) for j in 1:5],unique_genotypes=length(popcounts),effective_genotypes=diversity,
     founder_lineages=length(unique(Array(p.lineage)[live])),max_generation=n>0 ? maximum(Array(p.generation)[live]) : 0,
     minimum_physical=min(minimum(e[:,1:7]),minimum(c)),capacity_excess=maximum(c[:,1]+e[:,4])-p.parameters.capacity,
     hidden_maxabs=maximum(abs,h),finite=all(a->all(isfinite,a),arrays),
     resident_bytes=sum(sizeof,(p.environment,p.scratch,p.cells,p.hidden,p.actions,p.heat,p.work,p.genotype,p.lineage,p.generation,p.events,p.ids,p.parents,p.weights)),
     hypothetical_per_voxel_genome_bytes=prod(p.dims)*INPUTS*OUTPUTS*sizeof(Float32),shared_genotype_bytes=sizeof(p.weights))
end
function lattice_checkpoint(w,p,path)
    ispath(path)&&throw(ArgumentError("checkpoint exists"))
    serialize(path,(format="lattice-ecology-v1",source_sha=MODEL_SHA,source=MODEL_SOURCE,
        shared_sources=Interior.arena_sources(),julia_version=string(VERSION),grid=w.grid,seed=w.seed,time=w.time,state=map(Array,(p.environment,p.cells,p.hidden,p.actions,p.heat,p.work,p.genotype,p.lineage,p.generation,p.events,p.ids,p.parents,p.weights)),
        parameters=p.parameters,loci=p.loci,tick=p.tick))
end
function restore_lattice(path;backend=CPU())
    c=deserialize(path);c.format=="lattice-ecology-v1"&&c.source_sha==MODEL_SHA&&c.shared_sources==Interior.arena_sources()&&c.julia_version==string(VERSION)||throw(ArgumentError("exact lattice study source required"))
    w=lattice_world(c.grid;seed=c.seed,backend);p=LatticeEcology(w;parameters=c.parameters,loci=c.loci,occupancy=0,founder_genotypes=1,initial_mass=c.parameters.minimum_mass)
    for (a,b) in zip((p.environment,p.cells,p.hidden,p.actions,p.heat,p.work,p.genotype,p.lineage,p.generation,p.events,p.ids,p.parents,p.weights),c.state);copyto!(a,b);end
    p.tick=c.tick;w.time=c.time;w,p
end
end
