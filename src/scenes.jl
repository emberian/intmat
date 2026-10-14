"""Construct a reproducible art experiment. Parameters and rules are nondimensional.

The ellipsoid used for observation is a crop of a periodic computational domain,
not a simulated solid wall. Activity is seeded deliberately; no learned rule is claimed.
"""
function build_scene(config::AbstractDict; backend=CPU())
    sim = get(config,"simulation",Dict())
    n = Int(get(sim,"resolution",40)); dim=Int(get(sim,"dimension",3))
    dim in (2,3) || throw(ArgumentError("dimension must be 2 or 3"))
    seed=Int(get(sim,"seed",17)); rng=MersenneTwister(seed)
    g=Grid(n,n,dim==3 ? n : 1)
    w=World(g;seed,backend)
    activity=zeros(Float32,g.size); nutrient=similar(activity)
    # Five separated excitations and a weak helical seam; this is initial geometry.
    for I in CartesianIndices(activity)
        x,y,z=ntuple(d->(I[d]-.5f0)*g.spacing[d]-g.extent[d]/2,3)
        dim==2 && (z=0f0)
        value=0f0
        for l in 0:4
            angle=Float32(2pi*l/5)
            cx=.40f0*cos(angle); cy=.40f0*sin(angle); cz=dim==3 ? .3f0*sin(2angle) : 0f0
            value+=.95f0*exp(-((x-cx)^2+(y-cy)^2+(z-cz)^2)/.018f0)
        end
        activity[I]=value+.008f0*(rand(rng,Float32)-.5f0)
        nutrient[I]=.7f0+.3f0*exp(-2f0*(x*x+y*y+z*z))
    end
    copyto!(w.fields.activity,activity);copyto!(w.fields.nutrient,nutrient)
    processes=Any[]
    opts(name)= (; (Symbol(k)=>Float32(v) for (k,v) in get(config,name,Dict()) if k!="enabled")...)
    enabled(name,default=true)=get(get(config,name,Dict()),"enabled",default)
    enabled("wave") && push!(processes,WaveForcing(;opts("wave")...))
    if enabled("neural",false)
        neural=get(config,"neural",Dict())
        nopts=(; (Symbol(k)=>v for (k,v) in neural if k!="enabled")...)
        push!(processes,NeuralField(w;nopts...))
    end
    enabled("excitable") && push!(processes,Excitable(;opts("excitable")...))
    enabled("transport") && push!(processes,Diffusion(:nutrient,Float32(get(get(config,"transport",Dict()),"rate",0.002));material_feedback=Float32(get(get(config,"transport",Dict()),"material_feedback",0.0))))
    enabled("deposition") && push!(processes,Deposition(;opts("deposition")...))
    particles=nothing
    part=get(config,"particles",Dict())
    count=Int(get(part,"count",600))
    count>=0 || throw(ArgumentError("particle count must be nonnegative"))
    if count>0
        backend isa CPU || throw(ArgumentError("particle scenes currently require CPU; set particles.count=0 for GPU scenes"))
        rule=ParticleRule(; (Symbol(k)=>(k=="hidden" ? Int(v) : Float32(v)) for (k,v) in part if k!="count")...)
        positions=zeros(Float32,3,count)
        for i in 1:count
            # Seed mobile carriers inside an ellipsoid, in grid coordinates.
            p=randn(rng,Float32,3);p .*= rand(rng,Float32)^(1f0/3f0)*.7f0/sqrt(sum(abs2,p))
            for d in 1:3
                positions[d,i]=dim==2 && d==3 ? 0f0 : (p[d]+1f0)/g.spacing[d]
            end
        end
        hidden=.25f0.*randn(rng,Float32,size(rule.recurrence,1),count)
        particles=ParticleCarrier(positions;rule,hidden)
        push!(processes,particles)
    end
    (world=w,processes=Tuple(processes),particles=particles)
end

"""Optical interpretation of simulation fields; mask is a viewing crop only."""
function optical_fields(w::World)
    deposit=Array(w.fields.deposit);activity=Array(w.fields.activity)
    mask=similar(deposit)
    for I in CartesianIndices(mask)
        x,y,z=ntuple(d->2f0*(I[d]-.5f0)/size(mask,d)-1f0,3)
        size(mask,3)==1 && (z=0f0)
        r=sqrt((x/.87f0)^2+(y/.94f0)^2+(z/.80f0)^2)
        mask[I]=clamp((1f0-r)*30f0,0f0,1f0)
    end
    (density=.10f0 .+ 1.5f0.*deposit,activity=activity,structure=deposit,mask=mask)
end
