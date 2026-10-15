# O(1) nested-grid ownership descriptor, with no stored fine-grid index map.
struct NestedBlocks
    fine::Int
    coarse::Int
    ratio::Int
end
struct NestedBlock
    first::Int
    fine::Int
    ratio::Int
end
NestedBlocks(fine,coarse)=fine>0&&coarse>0&&fine%coarse==0 ?
    NestedBlocks(fine,coarse,fine÷coarse) : throw(ArgumentError("positive nested integer grids required"))
Base.length(b::NestedBlocks)=b.coarse^3
Base.eachindex(b::NestedBlocks)=Base.OneTo(length(b))
Base.eltype(::Type{NestedBlocks})=NestedBlock
function Base.getindex(b::NestedBlocks,i::Int)
    1<=i<=length(b)||throw(BoundsError(b,i))
    x=mod(i-1,b.coarse);y=mod((i-1)÷b.coarse,b.coarse);z=(i-1)÷b.coarse^2
    NestedBlock(1+b.ratio*(x+b.fine*(y+b.fine*z)),b.fine,b.ratio)
end
Base.iterate(b::NestedBlocks,i=1)=i>length(b) ? nothing : (b[i],i+1)
Base.length(b::NestedBlock)=b.ratio^3
Base.eltype(::Type{NestedBlock})=Int
Base.iterate(b::NestedBlock)=Base.iterate(b,(0,0,0))
function Base.iterate(b::NestedBlock,state::NTuple{3,Int})
    x,y,z=state;z>=b.ratio&&return nothing
    site=b.first+x+b.fine*(y+b.fine*z)
    x+=1
    if x==b.ratio;x=0;y+=1;end
    if y==b.ratio;y=0;z+=1;end
    site,(x,y,z)
end
