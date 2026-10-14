"""Evaluate independent recurrent controllers without copying each genome.

Input columns contain hidden state followed by sensors. This CPU implementation
retains the reference Float32 arithmetic; the kernel accepts backend-native
arrays for accelerator experiments without silently transferring state.
"""
function neural_controller!(outputs::AbstractMatrix, weights::AbstractArray{T,3}, marks::AbstractArray{T,3}, inputs::AbstractMatrix) where T
    size(weights)==size(marks) || throw(DimensionMismatch("controller marks"))
    rows,cols,n=size(weights)
    size(inputs)==(cols,n) && size(outputs)==(rows,n) || throw(DimensionMismatch("controller batch"))
    @inbounds for i in 1:n, row in 1:rows
        value=zero(T)
        for col in 1:cols
            value+=(weights[row,col,i]+marks[row,col,i])*inputs[col,i]
        end
        outputs[row,i]=tanh(value)
    end
    outputs
end
@kernel function neural_controller_kernel!(outputs,weights,marks,inputs,cols)
    row,i=@index(Global,NTuple)
    value=zero(eltype(outputs))
    for col in 1:cols
        value+=(weights[row,col,i]+marks[row,col,i])*inputs[col,i]
    end
    outputs[row,i]=tanh(value)
end
"""Backend-native batched controller. Caller owns transfers and synchronization."""
function neural_controller_backend!(outputs,weights,marks,inputs)
    rows,cols,n=size(weights)
    size(marks)==size(weights) && size(inputs)==(cols,n) && size(outputs)==(rows,n) || throw(DimensionMismatch("controller batch"))
    backend=KernelAbstractions.get_backend(outputs)
    all(a->KernelAbstractions.get_backend(a)==backend,(weights,marks,inputs)) || throw(ArgumentError("controller backend mismatch"))
    n==0 && return outputs
    neural_controller_kernel!(backend)(outputs,weights,marks,inputs,cols;ndrange=(rows,n))
    KernelAbstractions.synchronize(backend)
    outputs
end
"""Drop-in CPU controller outputs for the ecology's existing sensor vectors."""
function neural_controller_outputs(weights,marks,inputs::AbstractVector)
    rows,cols,n=size(weights)
    size(marks)==size(weights) && length(inputs)==n && all(x->length(x)==cols,inputs) || throw(DimensionMismatch("controller inputs"))
    outputs=[Vector{eltype(weights)}(undef,rows) for _ in 1:n]
    @inbounds for i in 1:n, row in 1:rows
        value=zero(eltype(weights))
        for col in 1:cols
            value+=(weights[row,col,i]+marks[row,col,i])*inputs[i][col]
        end
        outputs[i][row]=tanh(value)
    end
    outputs
end
