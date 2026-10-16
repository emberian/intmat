"""Pre-action usable-flow control variate; offline native-history gate.

This does not alter utility. Child gifts only undo native partition jumps,
using the same direct-parent convention as the current cohort estimator.
"""
module LocalFlowBaseline
using Interior: JSON3
using Statistics, Random

function adjusted(rows,events)
    gifts=Dict{Int,Vector{Tuple{Int,Float64}}}()
    for e in events
        push!(get!(gifts,Int(e.parent_id),Tuple{Int,Float64}[]),(round(Int,e.time/.02),Float64(e.child_usable_before_copy)))
    end
    [Float64(r.usable)+sum(g for (step,g) in get(gifts,Int(r.id),Tuple{Int,Float64}[]) if step<=Int(r.step);init=0.) for r in rows]
end
function features(rows,value;tau=1.)
    flow=zeros(length(rows));last=Dict{Int,Int}()
    for (i,r) in enumerate(rows)
        id=Int(r.id)
        if haskey(last,id)
            j=last[id];dt=(Int(r.step)-Int(rows[j].step))*.02
            dt>0||error("nonpositive resident clock")
            a=-expm1(-dt/tau);flow[i]=(1-a)*flow[j]+a*(value[i]-value[j])/dt
        end
        last[id]=i
    end
    flow
end
function run(input,out;horizon=4.)
    isfile(out)&&error("fresh output required")
    d=JSON3.read(read(input,String));rows=d.rows;y=adjusted(rows,d.events);flow=features(rows,y)
    lookup=Dict((Int(r.id),Int(r.step))=>i for (i,r) in enumerate(rows));lag=round(Int,horizon/.02)
    samples=Any[]
    for (i,r) in enumerate(rows)
        step=Int(r.step);id=Int(r.id)
        # Actual acting clock1, held-out20–28, labels through32.
        if 1000<=step<1400 && step%50==0 && haskey(lookup,(id,step+lag)) && haskey(lookup,(id,step-lag))
            j=lookup[(id,step+lag)];past=lookup[(id,step-lag)]
            target=y[j]-y[i]
            push!(samples,(;id,step,target,flow_baseline=horizon*flow[i],previous_return=y[i]-y[past],generation=r.generation))
        end
    end
    rms(v)=sqrt(mean(abs2,v))
    metrics=(;zero_rms=rms([s.target for s in samples]),previous_return_rms=rms([s.target-s.previous_return for s in samples]),
        flow_rms=rms([s.target-s.flow_baseline for s in samples]),samples=length(samples))
    # Gaussian score expectation: antithetic common draws remove odd numerical
    # fluctuation and demonstrate preservation of a known linear action effect.
    rng=MersenneTwister(211);noise=.3 .* randn(rng,10000);epsilon=vcat(noise,-noise)
    effect=1.2e-4;score=epsilon/.3^2
    r0=samples[1].target;b=samples[1].flow_baseline
    plain=(r0.+effect.*epsilon).*score;centered=(r0.+effect.*epsilon.-b).*score
    score_check=(;effect,expected_sample=effect*mean(epsilon.^2)/.3^2,plain_mean=mean(plain),centered_mean=mean(centered),
        modeled_plain_variance=metrics.zero_rms^2/.3^2+2effect^2,
        modeled_previous_variance=metrics.previous_return_rms^2/.3^2+2effect^2,
        modeled_flow_variance=metrics.flow_rms^2/.3^2+2effect^2,
        note="known linear-effect Gaussian fixture; variances use all248 observed offsets/residuals, not native noisy-policy samples")
    report=(;schema="local-flow-baseline-v1",horizon,flow_time=1.,metrics,score_check,samples,
        scope="offline control-variate gate using fixed-policy native histories; no native policy change or action-gradient claim",
        definition="pre-action exponentially averaged derivative of U plus cumulative direct-child partition gifts; baseline H*flow; reward stays existing birth-adjusted delta U")
    open(io->JSON3.write(io,report),out,"w");report
end
end
