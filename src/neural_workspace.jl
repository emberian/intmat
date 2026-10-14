"""Reusable, contiguous observation/controller storage; grid bins use CSR arrays.

No world or ecology state is owned here. Rebuild after population or grid changes;
otherwise frame preparation reuses every buffer and preserves synchronous sensing.
"""
mutable struct NeuralWorkspace
    dims::NTuple{3,Int}
    active::Int
    ids::Vector{Int}
    snapshot::Matrix{Float32}
    cells::Vector{CartesianIndex{3}}
    counts::Vector{Int}
    offsets::Vector{Int}
    cursors::Vector{Int}
    indices::Vector{Int}
    inputs::Matrix{Float32}
    outputs::Matrix{Float32}
end
function NeuralWorkspace(w,p)
    n=length(p.ids); h=size(p.hidden,1); dims=w.grid.size; nc=prod(dims)
    NeuralWorkspace(dims,n,zeros(Int,n),zeros(Float32,3,n),Vector{CartesianIndex{3}}(undef,n),zeros(Int,nc),zeros(Int,nc+1),zeros(Int,nc),zeros(Int,n),zeros(Float32,h+NEURAL_ECOLOGY_SENSORS,n),zeros(Float32,h+9,n))
end
@inline function _workspace_value(fields,channel,c)
    field=fields[channel]
    field===nothing ? (channel==6 ? 1f0 : 0f0) : field[c]
end
function workspace_sense!(w,p,s::NeuralWorkspace)
    n=length(p.ids); h=size(p.hidden,1); dims=w.grid.size
    if s.dims!=dims || size(s.snapshot,2)<n || size(s.inputs,1)!=h+NEURAL_ECOLOGY_SENSORS
        replacement=NeuralWorkspace(w,p)
        capacity=max(n,2size(s.snapshot,2))
        if capacity>n
            replacement.snapshot=zeros(Float32,3,capacity); replacement.cells=Vector{CartesianIndex{3}}(undef,capacity); replacement.indices=zeros(Int,capacity); replacement.ids=zeros(Int,capacity)
            replacement.inputs=zeros(Float32,h+NEURAL_ECOLOGY_SENSORS,capacity); replacement.outputs=zeros(Float32,h+9,capacity)
        end
        for field in fieldnames(NeuralWorkspace); setfield!(s,field,getfield(replacement,field)); end
    end
    s.active=n; copyto!(view(s.snapshot,:,1:n),p.positions); copyto!(view(s.ids,1:n),p.ids); fill!(s.counts,0)
    linear=LinearIndices(dims)
    @inbounds for i in 1:n
        c=_neural_cell(view(s.snapshot,:,i),dims); s.cells[i]=c; s.counts[linear[c]]+=1
    end
    s.offsets[1]=1
    @inbounds for c in eachindex(s.counts)
        s.offsets[c+1]=s.offsets[c]+s.counts[c]; s.cursors[c]=s.offsets[c]
    end
    @inbounds for i in 1:n
        c=linear[s.cells[i]]; s.indices[s.cursors[c]]=i; s.cursors[c]+=1
    end
    fields=(w.fields.nutrient,w.fields.deposit,ntuple(k->p.medium===nothing ? nothing : _neural_medium_field(p,(:signal,:phi,:interface,:permeability)[k]),4)...)
    _workspace_sensors!(p,s,fields)
    s
end
function _workspace_sensors!(p,s,fields)
    n=s.active; h=size(p.hidden,1); dims=s.dims; linear=LinearIndices(dims)
    @inbounds for i in 1:n
        c=s.cells[i]
        for k in 1:h; s.inputs[k,i]=p.hidden[k,i]; end
        for channel in 1:6; s.inputs[h+channel,i]=_neural_squash(_workspace_value(fields,channel,c)); end
        s.inputs[h+7,i]=_neural_squash(p.storage[i]); s.inputs[h+8,i]=_neural_squash(p.energy[i]); s.inputs[h+9,i]=p.age[i]/p.lifespan
        for channel in 1:6, d in 1:3
            gradient=0f0
            if dims[d]!=1
                lo=CartesianIndex(ntuple(k->k==d ? mod1(c[k]-1,dims[k]) : c[k],3))
                hi=CartesianIndex(ntuple(k->k==d ? mod1(c[k]+1,dims[k]) : c[k],3))
                gradient=_neural_squash((_workspace_value(fields,channel,hi)-_workspace_value(fields,channel,lo))/2)
            end
            s.inputs[h+9+(channel-1)*3+d,i]=gradient
        end
        r=p.neighbor_radius; dxsum=0f0; dysum=0f0; dzsum=0f0; count=0
        # First periodic occurrence of each axis reproduces the old Set traversal,
        # including small grids where the neighborhood wraps more than once.
        zrange=dims[3]==1 ? (0:0) : (-r:min(r,-r+dims[3]-1))
        for dz in zrange, dy in -r:min(r,-r+dims[2]-1), dx in -r:min(r,-r+dims[1]-1)
            q=CartesianIndex(mod1(c[1]+dx,dims[1]),mod1(c[2]+dy,dims[2]),mod1(c[3]+dz,dims[3])); bin=linear[q]
            for slot in s.offsets[bin]:(s.offsets[bin+1]-1)
                j=s.indices[slot]; j==i && continue
                x=Float32(mod(s.snapshot[1,j]-s.snapshot[1,i]+dims[1]/2,dims[1])-dims[1]/2)
                y=Float32(mod(s.snapshot[2,j]-s.snapshot[2,i]+dims[2]/2,dims[2])-dims[2]/2)
                z=dims[3]==1 ? 0f0 : Float32(mod(s.snapshot[3,j]-s.snapshot[3,i]+dims[3]/2,dims[3])-dims[3]/2)
                x*x+y*y+z*z>r*r && continue
                dxsum+=x; dysum+=y; dzsum+=z; count+=1
            end
        end
        s.inputs[h+28,i]=_neural_squash(count)
        denominator=Float32(max(count,1)*r)
        s.inputs[h+29,i]=dxsum/denominator; s.inputs[h+30,i]=dysum/denominator; s.inputs[h+31,i]=dzsum/denominator; s.inputs[h+32,i]=1f0
    end
    s
end
function workspace_prepare!(w,p,s::NeuralWorkspace)
    workspace_sense!(w,p,s)
    neural_controller!(view(s.outputs,:,1:s.active),p.weights,p.plastic_marks,view(s.inputs,:,1:s.active))
    s
end
"""Expose a frame for ecology execution; optional evaluation follows sensor edits."""
function neural_frame!(s::NeuralWorkspace,w,p;evaluate=true)
    workspace_sense!(w,p,s)
    evaluate && neural_controller!(view(s.outputs,:,1:s.active),p.weights,p.plastic_marks,view(s.inputs,:,1:s.active))
    (;ids=view(s.ids,1:s.active),cells=view(s.cells,1:s.active),bins=s,inputs=view(s.inputs,:,1:s.active),outputs=view(s.outputs,:,1:s.active))
end
function _neural_allocate!(field,s::NeuralWorkspace,requests)
    allocated=zeros(Float32,length(requests))
    @inbounds for bin in eachindex(s.counts)
        first=s.offsets[bin]; last=s.offsets[bin+1]-1
        total=0.; for slot in first:last; total+=Float64(requests[s.indices[slot]]); end
        available=Float64(field[bin]); total==0 && continue
        spent=0.
        for slot in first:last
            i=s.indices[slot]
            exact=total<=available ? Float64(requests[i]) : available*(Float64(requests[i])/total)
            amount=Float32(exact); Float64(amount)>exact && (amount=prevfloat(amount))
            allocated[i]=amount; spent+=Float64(amount)
        end
        field[bin]=Float32(max(0.,available-spent))
    end
    allocated
end
