#!/usr/bin/env julia
"""Study-only episodic action learning around the unchanged native coupled step.

Physical episodes are independent fresh branches. Readout marks persist between
training episodes; this is an engineered training assay, not autonomous lifetime
learning. All actions, uptake, secretion, diffusion and costs use production laws.
"""
module NativeActionCredit
using Interior, KernelAbstractions, Random, Statistics, SHA, TOML
using Interior: JSON3
include(joinpath(@__DIR__,"../scripts/native_constructed_garden.jl"))

# Explicit Gaussian action perturbation, not a fluctuation inferred from recurrent
# state. Bias input is exactly one in dc_sense!. No forward dynamics are duplicated.
@kernel function install!(marks,learned,noise,h,rows)
    i,c,a=@index(Global,NTuple)
    marks[i,c,h+rows[a]]=learned[i,c,a]+(c==size(learned,2) ? noise[i,a] : 0f0)
end
@kernel function eligible!(trace,inputs,noise,invvariance)
    i,c,a=@index(Global,NTuple)
    trace[i,c,a]+=inputs[i,c]*noise[i,a]*invvariance
end
@kernel function reinforce!(learned,trace,advantage,rate,bound)
    i,c,a=@index(Global,NTuple)
    learned[i,c,a]=clamp(learned[i,c,a]+rate*advantage[i]*trace[i,c,a],-bound,bound)
end

function fixture(profile,context;seed=211)
    cfg=TOML.parsefile(profile)
    cfg["seed"]=seed;cfg["grid"]["shape"]=[8,8,8];cfg["grid"]["extent"]=[1.,1.,1.]
    cfg["cellular"]["founders"]=4;cfg["cellular"]["capacity"]=64
    # A nondividing assay isolates within-controller change from founder sorting.
    cfg["cellular"]["growth_rate"]=0.;cfg["cellular"]["mutation"]=0.
    cfg["cellular"]["adaptation_strength"]=.3
    for section in ("resources","feedstock"),key in ("source_rates","source_budgets","sink_rates","sink_budgets")
        cfg[section][key].=0.
    end
    cfg["resources"]["initial"]=[context[1],context[2],0.]
    cfg["native"]["dt"]=.02
    cfg["native"]["energy_type"]="Float64";cfg["native"]["inventory_type"]="Float64"
    w,p,_=native_constructed_garden(CPU(),cfg)
    (;grid=w.grid,seed=w.seed,checkpoint=Interior.device_habitat_checkpoint(w,p),config=cfg)
end
usable(p)=[Float64(p.cells.state.energy[i])+Float64(p.cells.parameters.efficiency)*sum(Float64(p.cells.state.metabolites[i,k])*Float64(p.cells.parameters.potentials[k]) for k in 1:2) for i in 1:4]
audit(p)=Interior.device_habitat_audit(Interior.device_habitat_snapshot(p))
function episode(f,learned;rows=(4,5,7),sigma=.15,seed=1,horizon=4.,dt=.02)
    w=World(f.grid,CPU(),NamedTuple(),NamedTuple(),0.,f.seed)
    p=Interior.device_habitat_restore(w,f.checkpoint);s=p.cells.state
    h=Int(p.cells.geometry.hidden);C=Int(p.cells.geometry.capacity)
    size(learned)==(C,2h+Interior.CELLULAR_SENSORS,length(rows)) || error("readout shape")
    trace=zero(learned);noise=zeros(Float32,C,length(rows));rng=MersenneTwister(seed)
    before=audit(p);u0=usable(p);ids=copy(s.ids);steps=round(Int,horizon/dt)
    isapprox(steps*dt,horizon;atol=1e-8)||error("horizon must divide dt")
    started=time()
    for t in 1:steps
        randn!(rng,noise);noise .*=Float32(sigma)
        install!(CPU(),64)(s.marks,learned,noise,h,rows;ndrange=size(learned))
        Interior.device_cellular_step!(w,p.cells,dt;anabolism=p.habitat,check_status=false)
        # Inputs remain the exact pre-action observations used by dc_controller!.
        sigma>0 && eligible!(CPU(),64)(trace,s.inputs,noise,Float32(inv(sigma^2));ndrange=size(trace))
    end
    Interior.device_cellular_status!(p.cells)
    all(s.active[1:4].==1) && s.ids==ids || error("death/division makes this fixed-ID fixture invalid")
    after=audit(p);reward=usable(p).-u0
    material=after.material-before.material-after.supplied+before.supplied+after.removed-before.removed
    potential=after.potential-before.potential+after.work-before.work+after.heat-before.heat-after.supplied_potential+before.supplied_potential+after.removed_potential-before.removed_potential
    after.finite && abs(material)<1e-5 && abs(potential)<1e-5 || error("physical audit failed")
    ledger(k)=[Float64(s.ledger[i,k])+Float64(s.ledger_low[i,k]) for i in 1:4]
    (;reward,trace,work=ledger(1),uptake_a=ledger(3),uptake_b=ledger(4),
        mass=collect(s.mass[1:4]),material_residual=material,potential_residual=potential,
        elapsed_seconds=time()-started)
end
compact(x)=NamedTuple{Tuple(filter(!=(:trace),collect(keys(x))))}(x)
function run(profile,out;episodes=48,horizon=4.,rate=.003,sigma=.15,seed=211,noise_seed=73000)
    episodes>0 && isfinite(horizon) && horizon>0 && isfinite(rate) && rate>=0 && isfinite(sigma) && sigma>0 || error("invalid training parameters")
    ispath(out)&&error("fresh output required");mkpath(out)
    driver=read(@__FILE__);write(joinpath(out,"executed-driver.jl"),driver)
    contexts=((.1,.5),(.5,.1),(.2,.4),(.4,.2))
    fixtures=[fixture(profile,c;seed) for c in contexts]
    f=first(fixtures);h=Int(f.checkpoint.cells.snapshot.geometry.hidden)
    C=64;rows=(4,5,7);initial=zeros(Float32,C,2h+Interior.CELLULAR_SENSORS,3)
    before=[compact(episode(x,initial;sigma=0.,horizon)) for x in fixtures]
    arms=Dict{String,Any}()
    for (name,sign) in (("credit",1f0),("no_credit",0f0),("reversed_credit",-1f0))
        learned=copy(initial);baseline=hcat((x.reward for x in before[1:2])...)
        history=Any[]
        for trial in 1:episodes
            context=mod1(trial,2)
            result=episode(fixtures[context],learned;rows,sigma,seed=noise_seed+trial,horizon)
            advantage=zeros(Float32,C);advantage[1:4].=result.reward.-baseline[:,context]
            reinforce!(CPU(),64)(learned,result.trace,advantage,Float32(rate)*sign,.3f0;ndrange=size(learned))
            baseline[:,context].=.8baseline[:,context].+.2result.reward
            push!(history,(;trial,context,advantage=collect(advantage[1:4]),compact(result)...))
            trial%8==0 && (println(name," trial=",trial," return=",mean(result.reward)," marks=",maximum(abs,learned));flush(stdout))
        end
        after=[compact(episode(x,learned;sigma=0.,horizon)) for x in fixtures]
        arms[name]=(;history,after,learned=collect(learned),mark_max=maximum(abs,learned),
            heldout_delta=[after[j].reward.-before[j].reward for j in 3:4])
        open(io->JSON3.write(io,arms[name]),joinpath(out,name*".json"),"w")
    end
    report=(;schema="native-action-credit-v1",episodes,horizon,rate,sigma,seed,noise_seed,rows,contexts,before,
        source_pin=f.checkpoint.source,cellular_source_pin=f.checkpoint.cells.source,driver_sha256=bytes2hex(sha256(driver)),
        factory_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"../scripts/native_constructed_garden.jl")))),config=f.config,
        scope="CPU 8^3 coupled production physics, four nondividing founders, independent physical episodes; only readout marks persist across episodes; no autonomous lifetime-learning or arena-speed claim",
        reward="terminal minus initial battery + efficiency * stored food potential; physical sources zero; no reward creates matter or energy",
        summary=Dict(k=>(;mark_max=v.mark_max,heldout_delta=v.heldout_delta,
            heldout_mean_delta=mean(Iterators.flatten(v.heldout_delta))) for (k,v) in arms))
    open(io->JSON3.pretty(io,report),joinpath(out,"report.json"),"w")
    report
end
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    length(ARGS)>=2||error("PROFILE FRESH_OUTPUT [EPISODES=48] [HORIZON=4]")
    NativeActionCredit.run(ARGS[1],ARGS[2];episodes=length(ARGS)>2 ? parse(Int,ARGS[3]) : 48,horizon=length(ARGS)>3 ? parse(Float64,ARGS[4]) : 4.)
end
