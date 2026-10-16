"""Local input homeostasis with independently fixed temporal poles.

Standalone controller prototype. The toy work account is NOT native inventory.
No recurrent native weights or physical laws are changed.
"""
module HomeostaticTemporalMemory
using LinearAlgebra, Statistics
using Interior: JSON3

mutable struct State
    center::Float64
    variance::Float64
    fast::Float64
    slow::Float64
    battery::Float64
    work::Float64
end
State(;battery=1.)=State(0.,.25^2/2,0.,0.,battery,0.)
options(;dt=.2,homeostasis=20.)=(dt=dt,ah=-expm1(-dt/homeostasis),
    af=-expm1(-dt/.4),as=-expm1(-dt/2.),floor=.01^2,kappa=1e-4)
gain(s,q)=inv(sqrt(max(q.floor,s.variance)))
function observe!(s,u,q;adapt=true,yoked_payment=0.)
    # Causal drive uses the previously installed homeostat.
    v=tanh((u-s.center)*gain(s,q))
    s.fast+=q.af*(v-s.fast);s.slow+=q.as*(v-s.slow)
    delta=u-s.center
    center=s.center+q.ah*delta
    variance=(1-q.ah)*(s.variance+q.ah*delta^2)
    proposal=q.kappa*(abs(center-s.center)+abs(log(inv(sqrt(max(q.floor,variance)))/gain(s,q))))
    paid=0.
    if adapt && proposal<=s.battery
        s.center=center;s.variance=variance;s.battery-=proposal;s.work+=proposal;paid=proposal
    elseif !adapt
        paid=min(yoked_payment,s.battery);s.battery-=paid;s.work+=paid
    end
    (;features=(1.,s.fast,s.slow),paid,drive=v)
end

function replay(;adapt=true,yoked=nothing,battery=1.,homeostasis=20.)
    q=options(;homeostasis);s=State(;battery)
    nphase=800;n=3nphase;features=zeros(n,3);target=zeros(n)
    payments=zeros(n);centers=zeros(n);gains=zeros(n);inputs=zeros(n)
    for t in 1:n
        phase=cld(t,nphase);offset=phase==2 ? .6 : 0.;amplitude=phase==2 ? .025 : .25
        clock=(t-1)*q.dt
        u=offset+amplitude*sin(2pi*clock/8)
        r=observe!(s,u,q;adapt,yoked_payment=yoked===nothing ? 0. : yoked[t])
        features[t,:].=r.features;payments[t]=r.paid;centers[t]=s.center;gains[t]=gain(s,q);inputs[t]=u
        # Forecast relative phase one physical time ahead, NOT absolute stock.
        target[t]=sin(2pi*(clock+1.)/8)
    end
    # Equal-size linear readouts, separately calibrated on the first regime.
    train=401:800;X=features[train,:];w=(X'X+1e-6I)\(X'target[train])
    prediction=features*w
    metrics=[(;phase,late_mse=mean((prediction[(phase-1)*nphase+401:phase*nphase].-target[(phase-1)*nphase+401:phase*nphase]).^2),
        early_mse=mean((prediction[(phase-1)*nphase+1:(phase-1)*nphase+200].-target[(phase-1)*nphase+1:(phase-1)*nphase+200]).^2)) for phase in 1:3]
    (;features,target,prediction,payments,centers,gains,inputs,weights=w,metrics,work=s.work,battery=s.battery)
end
function run(out)
    isfile(out)&&error("fresh output required")
    a=replay();fixed=replay(adapt=false,yoked=a.payments);denied=replay(battery=0.)
    # Too-fast adaptation is an explicit anti-control: it removes task signal.
    fast=replay(homeostasis=.2)
    report=(;schema="homeostatic-temporal-memory-v1",scope="bounded sensory replay; toy work account, no native physical simulation or reward-trained recurrent weights",
        dt=.2,filter_times=(.4,2.),homeostasis_time=20.,regimes=((offset=0.,amplitude=.25),(offset=.6,amplitude=.025),(offset=0.,amplitude=.25)),
        adaptive=(metrics=a.metrics,weights=a.weights,work=a.work),fixed_yoked=(metrics=fixed.metrics,weights=fixed.weights,work=fixed.work),
        unpaid=(metrics=denied.metrics,weights=denied.weights,work=denied.work),fast_homeostasis=(metrics=fast.metrics,weights=fast.weights,work=fast.work),
        samples=[(;time=(t-1)*.2,input=a.inputs[t],target=a.target[t],adaptive=a.prediction[t],fixed=fixed.prediction[t],
            center=a.centers[t],gain=a.gains[t]) for t in 1:10:length(a.target)])
    open(io->JSON3.write(io,report),out,"w");report
end
end
