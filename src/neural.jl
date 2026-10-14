"""Persistent two-channel neural field with local recurrent tanh dynamics.

`mode=:oscillator` uses an oscillatory recurrent prior; `:random` samples fixed
weights reproducibly from `seed`. These are configurable dynamical rules, not
trained networks. Latent state belongs to this process and survives `step!` calls.
The first channel writes physical activity; nutrient, deposit and pressure feed
back into both channels. Material conversion remains the reaction process's job.

`recurrence` is `(w11,w12,w21,w22)` and `bias` contains two channel offsets.
Use alongside Deposition and transport in either two or three dimensions.
"""
struct NeuralField{A}
    latent1::A
    latent2::A
    scratch1::A
    scratch2::A
    recurrence::NTuple{4,Float32}
    bias::NTuple{2,Float32}
    gain::Float32
    response_time::Float32
    diffusion::Float32
    material_feedback::Float32
    nutrient_feedback::Float32
    pressure_feedback::Float32
end
function NeuralField(w::World; seed=w.seed, mode=:oscillator, recurrence=nothing,
                     bias=(0f0,0f0), gain=1f0, response_time=0.5f0,
                     diffusion=0.003f0, material_feedback=0.7f0,
                     nutrient_feedback=0.3f0, pressure_feedback=0.4f0,
                     initial_noise=0.04f0)
    mode=Symbol(mode)
    mode in (:oscillator,:random) || throw(ArgumentError("neural mode must be oscillator or random"))
    isfinite(initial_noise) && initial_noise>=0 || throw(ArgumentError("initial_noise must be finite and nonnegative"))
    rng=MersenneTwister(seed)
    weights=recurrence===nothing ? (mode==:oscillator ? (1.6f0,-1.2f0,1.2f0,1.6f0) : Tuple(1.2f0 .* randn(rng,Float32,4))) : Tuple(Float32.(recurrence))
    length(weights)==4 || throw(ArgumentError("recurrence requires four weights"))
    offsets=Tuple(Float32.(bias))
    length(offsets)==2 || throw(ArgumentError("bias requires two values"))
    arrays=ntuple(_->KernelAbstractions.zeros(w.backend,Float32,w.grid.size...),4)
    a=Array(w.fields.activity)
    copyto!(arrays[1],tanh.(a .+ Float32(initial_noise).*randn(rng,Float32,w.grid.size...)))
    copyto!(arrays[2],tanh.(Float32(initial_noise).*randn(rng,Float32,w.grid.size...)))
    p=NeuralField(arrays...,weights,offsets,Float32(gain),Float32(response_time),Float32(diffusion),Float32(material_feedback),Float32(nutrient_feedback),Float32(pressure_feedback))
    validate(p)
    p
end
function validate(p::NeuralField)
    all(isfinite,(p.recurrence...,p.bias...,p.gain,p.response_time,p.diffusion,p.material_feedback,p.nutrient_feedback,p.pressure_feedback)) || throw(ArgumentError("neural parameters must be finite"))
    p.response_time>0 && p.diffusion>=0 && p.gain>=0 || throw(ArgumentError("neural response_time must be positive; diffusion and gain nonnegative"))
end
function validate(w::World,p::NeuralField)
    validate(p)
    for a in (p.latent1,p.latent2,p.scratch1,p.scratch2)
        size(a)==w.grid.size || throw(ArgumentError("neural latent grid does not match world"))
        KernelAbstractions.get_backend(a)==w.backend || throw(ArgumentError("neural latent backend does not match world"))
    end
end
# Convex explicit mixing of neighboring states and a bounded tanh target.
function stability_limit(w,p::NeuralField)
    spatial=2f0*p.diffusion*sum(inv(w.grid.spacing[d]^2) for d in 1:(w.grid.size[3]==1 ? 2 : 3))
    0.9f0/(inv(p.response_time)+spatial)
end
@kernel function neural_kernel!(uout,vout,activity,u,v,n,d,pressure,h,dt,weights,bias,gain,tau,D,mfeedback,nfeedback,pfeedback)
    i,j,k=@index(Global,NTuple)
    x=u[i,j,k]; y=v[i,j,k]
    material=d[i,j,k]/(1f0+d[i,j,k])
    food=n[i,j,k]/(1f0+n[i,j,k])
    sense=nfeedback*food-mfeedback*material+pfeedback*tanh(pressure[i,j,k])
    tx=tanh(weights[1]*x+weights[2]*y+bias[1]+sense)
    ty=tanh(weights[3]*x+weights[4]*y+bias[2]-sense)
    un=x+dt*((tx-x)/tau+D*laplacian(u,i,j,k,h))
    vn=y+dt*((ty-y)/tau+D*laplacian(v,i,j,k,h))
    uout[i,j,k]=un; vout[i,j,k]=vn
    activity[i,j,k]=gain*un
end
function advance!(w::World,p::NeuralField,dt::Real)
    check_dt(w,p,dt)
    f=w.fields
    launch!(w,neural_kernel!,p.scratch1,p.scratch2,f.activity,p.latent1,p.latent2,f.nutrient,f.deposit,f.pressure,w.grid.spacing,Float32(dt),p.recurrence,p.bias,p.gain,p.response_time,p.diffusion,p.material_feedback,p.nutrient_feedback,p.pressure_feedback)
    copyto!(p.latent1,p.scratch1); copyto!(p.latent2,p.scratch2)
    w
end
