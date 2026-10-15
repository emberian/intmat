# Cell slots are contiguous in Julia's column-major arrays. Grow directly into
# the final array instead of routing dense genotype storage through generic cat.
# New slots are intentionally uninitialized and must be filled before returning
# from the complete lifecycle operation; existing slots are copied exactly once.
function _cellular_append_slots(old::Array{T,N},count::Integer) where {T,N}
    count>=0 || throw(ArgumentError("nonnegative appended cell count required"))
    dimensions=ntuple(d->d==N ? size(old,d)+count : size(old,d),Val(N))
    result=Array{T,N}(undef,dimensions)
    copyto!(result,1,old,1,length(old))
    result
end

@inline function _cellular_copy_slot!(array::Array{T,N},destination,source) where {T,N}
    width=length(array)÷size(array,N)
    copyto!(array,(destination-1)*width+1,array,(source-1)*width+1,width)
    nothing
end

"""Reusable dense genome slots, addressed in logical cell order through `slots`."""
mutable struct CellularGenomeBank
    rows::Int
    columns::Int
    capacity::Int
    slots::Vector{Int}
    free_slots::Vector{Int}
    next_slot::Int
    weights::Vector{Float32}
    marks::Vector{Float32}
end

struct CellularGenomeView{Kind} <: AbstractArray{Float32,3}
    bank::CellularGenomeBank
end
Base.IndexStyle(::Type{<:CellularGenomeView})=IndexLinear()
Base.size(a::CellularGenomeView)=(a.bank.rows,a.bank.columns,length(a.bank.slots))
@inline _cellular_genome_data(a::CellularGenomeView{:weights})=a.bank.weights
@inline _cellular_genome_data(a::CellularGenomeView{:marks})=a.bank.marks
Base.dataids(a::CellularGenomeView)=Base.dataids(_cellular_genome_data(a))

function CellularGenomeBank(weights::Array{Float32,3},marks::Array{Float32,3})
    size(weights)==size(marks) || throw(DimensionMismatch("genome and mark storage"))
    rows,columns,n=size(weights);capacity=nextpow(2,max(n,1));block=rows*columns
    bank=CellularGenomeBank(rows,columns,capacity,collect(1:n),Int[],n+1,
        zeros(Float32,block*capacity),zeros(Float32,block*capacity))
    copyto!(bank.weights,1,weights,1,length(weights))
    copyto!(bank.marks,1,marks,1,length(marks))
    bank
end

Base.@propagate_inbounds function Base.getindex(a::CellularGenomeView,i::Int)
    @boundscheck checkbounds(a,i)
    bank=a.bank;block=bank.rows*bank.columns;cell,offset=divrem(i-1,block)
    @inbounds _cellular_genome_data(a)[(bank.slots[cell+1]-1)*block+offset+1]
end
Base.@propagate_inbounds function Base.getindex(a::CellularGenomeView,row::Int,column::Int,cell::Int)
    @boundscheck checkbounds(a,row,column,cell)
    bank=a.bank
    @inbounds _cellular_genome_data(a)[(bank.slots[cell]-1)*bank.rows*bank.columns+(column-1)*bank.rows+row]
end
Base.@propagate_inbounds function Base.setindex!(a::CellularGenomeView,value,i::Int)
    @boundscheck checkbounds(a,i)
    bank=a.bank;block=bank.rows*bank.columns;cell,offset=divrem(i-1,block)
    @inbounds _cellular_genome_data(a)[(bank.slots[cell+1]-1)*block+offset+1]=value
    a
end
Base.@propagate_inbounds function Base.setindex!(a::CellularGenomeView,value,row::Int,column::Int,cell::Int)
    @boundscheck checkbounds(a,row,column,cell)
    bank=a.bank
    @inbounds _cellular_genome_data(a)[(bank.slots[cell]-1)*bank.rows*bank.columns+(column-1)*bank.rows+row]=value
    a
end
Base.similar(::CellularGenomeView,::Type{T},dims::Dims) where T=Array{T}(undef,dims)
function Base.copy(a::CellularGenomeView)
    result=Array{Float32,3}(undef,size(a));bank=a.bank;block=bank.rows*bank.columns
    for (i,slot) in enumerate(bank.slots)
        copyto!(result,(i-1)*block+1,_cellular_genome_data(a),(slot-1)*block+1,block)
    end
    result
end
Base.Array(a::CellularGenomeView)=copy(a)
function Base.fill!(a::CellularGenomeView,value)
    bank=a.bank;block=bank.rows*bank.columns;data=_cellular_genome_data(a)
    for slot in bank.slots
        fill!(view(data,(slot-1)*block+1:slot*block),value)
    end
    a
end
function Base.all(::typeof(isfinite),a::CellularGenomeView)
    bank=a.bank;block=bank.rows*bank.columns;data=_cellular_genome_data(a)
    for slot in bank.slots
        all(isfinite,view(data,(slot-1)*block+1:slot*block)) || return false
    end
    true
end

function _cellular_reserve_genomes!(bank,required)
    required<=bank.capacity && return bank
    capacity=max(required,2bank.capacity);block=bank.rows*bank.columns
    for data in (bank.weights,bank.marks)
        old=length(data);resize!(data,capacity*block)
        fill!(view(data,old+1:length(data)),0f0)
    end
    bank.capacity=capacity
    bank
end
function _cellular_append_genomes!(bank,count)
    required=max(0,count-length(bank.free_slots))
    _cellular_reserve_genomes!(bank,bank.next_slot+required-1)
    for _ in 1:count
        slot=if isempty(bank.free_slots)
            value=bank.next_slot;bank.next_slot+=1;value
        else
            pop!(bank.free_slots)
        end
        push!(bank.slots,slot)
    end
    bank
end
function _cellular_select_genomes!(bank::CellularGenomeBank,indices)
    keep=falses(length(bank.slots));keep[indices].=true
    slots=bank.slots[indices]
    for (i,slot) in enumerate(bank.slots)
        keep[i] || push!(bank.free_slots,slot)
    end
    bank.slots=slots
    bank
end
@inline function _cellular_copy_slot!(array::CellularGenomeView,destination,source)
    bank=array.bank;block=bank.rows*bank.columns;data=_cellular_genome_data(array)
    copyto!(data,(bank.slots[destination]-1)*block+1,data,(bank.slots[source]-1)*block+1,block)
    nothing
end

function _cellular_mutate_genome!(rng,weights,slot,mutation,probability)
    for column in axes(weights,2),row in axes(weights,1)
        rand(rng)<probability && (weights[row,column,slot]+=Float32(mutation)*randn(rng,Float32))
    end
    nothing
end
function _cellular_mutate_genome!(rng,weights::CellularGenomeView,slot,mutation,probability)
    bank=weights.bank;start=(bank.slots[slot]-1)*bank.rows*bank.columns
    data=_cellular_genome_data(weights)
    for column in 1:bank.columns,row in 1:bank.rows
        index=start+(column-1)*bank.rows+row
        rand(rng)<probability && (data[index]+=Float32(mutation)*randn(rng,Float32))
    end
    nothing
end

# Resolve the physical slot once per cell, preserving the reference Float32
# accumulation order for every output neuron.
function neural_controller!(outputs::AbstractMatrix{Float32},weights::CellularGenomeView{:weights},
        marks::CellularGenomeView{:marks},inputs::AbstractMatrix{Float32})
    size(weights)==size(marks) || throw(DimensionMismatch("controller marks"))
    rows,columns,n=size(weights)
    size(inputs)==(columns,n) && size(outputs)==(rows,n) || throw(DimensionMismatch("controller batch"))
    wb=weights.bank;mb=marks.bank;block=rows*columns
    @inbounds for i in 1:n
        wo=(wb.slots[i]-1)*block;mo=(mb.slots[i]-1)*block
        for row in 1:rows
            value=0f0
            for column in 1:columns
                offset=(column-1)*rows+row
                value+=(wb.weights[wo+offset]+mb.marks[mo+offset])*inputs[column,i]
            end
            outputs[row,i]=tanh(value)
        end
    end
    outputs
end

function _cellular_enable_indexed_storage!(p)
    p.weights isa CellularGenomeView && return p
    bank=CellularGenomeBank(p.weights,p.plastic_marks)
    p.weights=CellularGenomeView{:weights}(bank)
    p.plastic_marks=CellularGenomeView{:marks}(bank)
    p
end
function _cellular_validate_storage(p)
    weights=p.weights;marks=p.plastic_marks
    if weights isa CellularGenomeView || marks isa CellularGenomeView
        weights isa CellularGenomeView{:weights} && marks isa CellularGenomeView{:marks} &&
            weights.bank===marks.bank || throw(ArgumentError("genome bank ownership mismatch"))
        bank=weights.bank;block=bank.rows*bank.columns
        !Base.mightalias(bank.weights,bank.marks) || throw(ArgumentError("genotype and developmental marks must own separate storage"))
        bank.rows>0 && bank.columns>0 && 1<=bank.next_slot<=bank.capacity+1 &&
            length(bank.weights)==length(bank.marks)==bank.capacity*block || throw(ArgumentError("invalid genome bank storage"))
        length(bank.slots)+length(bank.free_slots)==bank.next_slot-1 || throw(ArgumentError("invalid genome slot accounting"))
        seen=falses(bank.next_slot-1)
        for slots in (bank.slots,bank.free_slots),slot in slots
            1<=slot<bank.next_slot && !seen[slot] || throw(ArgumentError("invalid or aliased genome slot"))
            seen[slot]=true
        end
    else
        weights isa Array{Float32,3} && marks isa Array{Float32,3} || throw(ArgumentError("unsupported cellular genome storage"))
    end
    nothing
end
function _cellular_grow_genomes!(p,count)
    if p.weights isa CellularGenomeView
        _cellular_append_genomes!(p.weights.bank,count)
    else
        p.weights=_cellular_append_slots(p.weights,count)
        p.plastic_marks=_cellular_append_slots(p.plastic_marks,count)
    end
    nothing
end
function _cellular_select_storage!(p,indices)
    if p.weights isa CellularGenomeView
        _cellular_select_genomes!(p.weights.bank,indices)
    else
        p.weights=p.weights[:,:,indices]
        p.plastic_marks=p.plastic_marks[:,:,indices]
    end
    nothing
end

function _cellular_storage_statistics(p)
    if p.weights isa CellularGenomeView
        bank=p.weights.bank;block=bank.rows*bank.columns
        (;mode="indexed",active_cells=length(bank.slots),capacity=bank.capacity,
          reusable_slots=length(bank.free_slots),
          active_genome_bytes=2block*length(bank.slots)*sizeof(Float32),
          reserved_genome_bytes=(length(bank.weights)+length(bank.marks))*sizeof(Float32),
          slot_index_bytes=(length(bank.slots)+length(bank.free_slots))*sizeof(Int))
    else
        bytes=(length(p.weights)+length(p.plastic_marks))*sizeof(Float32)
        (;mode="dense",active_cells=length(p.ids),capacity=length(p.ids),reusable_slots=0,
          active_genome_bytes=bytes,reserved_genome_bytes=bytes,slot_index_bytes=0)
    end
end
