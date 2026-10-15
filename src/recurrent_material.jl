"""Eight-channel, fixed recurrent material. Four coupled oscillator pairs span
fast/slow clocks. Local state writes persistent solid from finite nutrient;
solid changes symmetric spatial communication and recurrent thresholds. No
training or reproductive evolution is claimed. GPU arrays remain resident."""
mutable struct RecurrentMaterial{A,R}
    state::A
    scratch::A
    resource_buffers::R
    coupling::Float32
    feedback::Float32
    writing::Float32
    erosion::Float32
    contrast::Float32
end
function RecurrentMaterial(w::World;seed=w.seed,coupling=.12,feedback=1.2,writing=.25,erosion=.06,contrast=.7)
    rng=MersenneTwister(seed)
    host=clamp.(.12f0.*randn(rng,Float32,w.grid.size...,8),-1f0,1f0)
    a=KernelAbstractions.zeros(w.backend,Float32,w.grid.size...,8)
    b=similar(a);copyto!(a,host);fill!(b,0)
    p=RecurrentMaterial(a,b,(w.fields.nutrient,w.scratch.nutrient,w.fields.deposit,w.scratch.deposit),Float32.((coupling,feedback,writing,erosion,contrast))...)
    validate(w,p)
    p
end
function validate(w::World,p::RecurrentMaterial)
    (w.fields.nutrient===p.resource_buffers[1] || w.fields.nutrient===p.resource_buffers[2]) &&
    (w.fields.deposit===p.resource_buffers[3] || w.fields.deposit===p.resource_buffers[4]) ||
        throw(ArgumentError("recurrent material resource buffers were replaced; fields are borrowed until the next stage"))
    p.state!==p.scratch || throw(ArgumentError("latent buffers must be distinct"))
    size(p.state)==size(p.scratch) || throw(DimensionMismatch("latent scratch shape"))
    for array in (p.state,p.scratch,p.resource_buffers...)
        KernelAbstractions.get_backend(array)==w.backend || throw(ArgumentError("recurrent material backend mismatch"))
    end
    for i in 1:4
        size(p.resource_buffers[i])==w.grid.size || throw(DimensionMismatch("resource buffer shape"))
        for j in 1:i-1
            p.resource_buffers[i]!==p.resource_buffers[j] || throw(ArgumentError("resource buffers must be distinct"))
        end
    end
    (w.scratch.nutrient===p.resource_buffers[1] || w.scratch.nutrient===p.resource_buffers[2]) &&
    (w.scratch.deposit===p.resource_buffers[3] || w.scratch.deposit===p.resource_buffers[4]) &&
    w.fields.nutrient!==w.scratch.nutrient && w.fields.deposit!==w.scratch.deposit ||
        throw(ArgumentError("resource source and destination must use opposite owned buffers"))
    size(p.state)==(w.grid.size...,8) || throw(DimensionMismatch("recurrent material state"))
    all(x->isfinite(x)&&x>=0,(p.coupling,p.feedback,p.writing,p.erosion,p.contrast)) || throw(ArgumentError("invalid recurrent material rates"))
end
function stability_limit(w,p::RecurrentMaterial)
    spatial=sum((inv(Float64(w.grid.spacing[d])^2) for d in 1:3 if w.grid.size[d]>1);init=0.0)
    min(.05f0,.9/(1+2p.coupling*spatial),.9/(p.writing+.05spatial+eps(Float32)))
end
@inline function recurrent_face(a,b)
    1f0/(1f0+4f0*(a+b))
end
@kernel function recurrent_material_kernel!(out,state,nout,dout,n,d,activity,h,dt,coupling,feedback,writing,erosion,contrast)
    i,j,k=@index(Global,NTuple)
    nx,ny,nz=size(n); pos=(i,j,k); dims=(nx,ny,nz)
    solid=d[i,j,k]; food=n[i,j,k]
    # Two constitutions meet without an imposed temporal forcing.
    constitution=tanh(8f0*(Float32(i)/Float32(nx)-.5f0))
    memory=solid/(1f0+solid)
    drive=0f0
    for channel in 1:8
        pair=(channel+1)÷2; first=2pair-1; second=2pair
        x=state[i,j,k,first]; y=state[i,j,k,second]
        previous=pair==1 ? 7 : first-2
        clock=pair==1 ? 1f0 : pair==2 ? .43f0 : pair==3 ? .17f0 : .065f0
        lateral=0f0
        value=state[i,j,k,channel]
        for axis in 1:3
            if dims[axis]>1
                for side in (-1,1)
                    q=phase_neighbor(pos,axis,side,dims)
                    lateral+=recurrent_face(solid,d[q...])*(state[q...,channel]-value)/h[axis]^2
                end
            end
        end
        odd=isodd(channel)
        # Rotational recurrent prior, coupled across timescales. Mineral changes
        # the threshold and cross-timescale coupling, not just its appearance.
        recurrent=odd ? 1.55f0*x-(1.9f0+contrast*constitution)*y :
            (1.9f0-contrast*constitution)*x+1.55f0*y
        cross=.28f0*state[i,j,k,previous]*(1f0+feedback*memory)
        bias=odd ? -.5f0*feedback*memory+.06f0*tanh(food) : .25f0*feedback*memory
        target=tanh(recurrent+cross+bias)
        updated=value+dt*(clock*(target-value)+coupling*clock*lateral)
        out[i,j,k,channel]=updated
        channel==1 && (drive=updated)
    end
    # Opposite halves of the fast cycle build and erase. Both exchanges are
    # funded locally, with equal/opposite resource entries; no source term.
    growth=min(food,dt*writing*max(0f0,drive-.15f0)*food)
    loss=min(solid,dt*erosion*max(0f0,-drive-.1f0)*solid)
    nf=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                q=phase_neighbor(pos,axis,side,dims)
                nf+=.025f0*recurrent_face(solid,d[q...])*(n[q...]-food)/h[axis]^2
            end
        end
    end
    nout[i,j,k]=food-growth+loss+dt*nf
    dout[i,j,k]=solid+growth-loss
    activity[i,j,k]=drive
end
function advance!(w::World,p::RecurrentMaterial,dt::Real)
    check_dt(w,p,dt)
    recurrent_material_kernel!(w.backend,256)(p.scratch,p.state,w.scratch.nutrient,w.scratch.deposit,
        w.fields.nutrient,w.fields.deposit,w.fields.activity,w.grid.spacing,Float32(dt),
        p.coupling,p.feedback,p.writing,p.erosion,p.contrast;ndrange=w.grid.size)
    p.state,p.scratch=p.scratch,p.state
    oldn=w.fields.nutrient; oldd=w.fields.deposit
    w.fields=merge(w.fields,(nutrient=w.scratch.nutrient,deposit=w.scratch.deposit))
    w.scratch=merge(w.scratch,(nutrient=oldn,deposit=oldd))
    w
end
function record_process_fields(p::RecurrentMaterial)
    (state=p.state,coupling=p.coupling,feedback=p.feedback,writing=p.writing,erosion=p.erosion,contrast=p.contrast)
end

process_ports(::RecurrentMaterial)=ProcessPorts(reads=(:nutrient,:deposit),writes=(:nutrient,:deposit,:activity),publishes=(:activity,),owned=(:state,),memory="eight shared-basis recurrent channels and persistent world deposit; resource arrays ping-pong, borrowed references expire at advance",law="four coupled recurrent oscillator pairs; material-dependent symmetric communication and conservative nutrient/solid conversion")

# The zero-copy resource swap invalidates fixed aliases held by other stages.
# These traits inspect references only; ownership/advancement stays with
# process_children, so borrowed media and reservoirs are never advanced twice.
composition_rebound_arrays(p::RecurrentMaterial)=p.resource_buffers
composition_borrowed_arrays(p::ActiveMedium)=(p.membrane,)
composition_borrowed_arrays(p::ResourcePool)=p.fields
composition_references(p::ResourcePool)=p.medium===nothing ? () : (p.medium,)
composition_references(p::NeuralEcology)=p.medium===nothing ? () : (p.medium,)
composition_references(p::MetabolicEcology)=(p.controller,p.resources,p.boundary)
composition_references(p::CellularEcology)=(p.resources,p.medium)
