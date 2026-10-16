"""Small additive statistics for a local conditional-action response model.

This is a study component, not an installed living policy. Each source owns its
original statistics. Never accumulate a neighbor's already pooled statistics.
"""
module LocalActionEvidence
using LinearAlgebra
export Evidence, observe!, fit, merge_originals
mutable struct Evidence
    gram::Matrix{Float64}
    rhs::Vector{Float64}
    count::Int
end
Evidence()=Evidence(zeros(3,3),zeros(3),0)
function observe!(s::Evidence, innovation, context, target)
    length(context)==2 || throw(ArgumentError("two context coordinates required"))
    v=Float64[innovation,innovation*context[1],innovation*context[2]]
    all(isfinite,v)&&isfinite(target)||throw(ArgumentError("nonfinite observation"))
    s.gram .+= v*v';s.rhs .+= v*target;s.count+=1
    s
end
"""Merge disjoint originals only; caller is responsible for stable source identity."""
function merge_originals(sources)
    out=Evidence()
    for s in sources
        out.gram .+= s.gram;out.rhs .+= s.rhs;out.count+=s.count
    end
    out
end
"""Shared context slopes with a separate, regularized action intercept per source.

Schur elimination retains source identity; simply merging all Grams would impose
one common intercept. The recipient index selects its private intercept. This is
ridge response estimation, not an unbiased policy-gradient claim.
"""
function fit(sources,recipient;intercept_ridge=1e-6,slope_ridge=1e-6)
    1<=recipient<=length(sources)||throw(ArgumentError("invalid recipient"))
    intercept_ridge>0&&slope_ridge>0||throw(ArgumentError("positive ridge required"))
    S=Matrix{Float64}(I,2,2)*slope_ridge;t=zeros(2)
    for s in sources
        a=s.gram[1,1]+intercept_ridge;c=s.gram[2:3,1]
        S .+= s.gram[2:3,2:3]-(c*c')/a
        t .+= s.rhs[2:3]-c*(s.rhs[1]/a)
    end
    slope=Symmetric(S)\t;s=sources[recipient]
    intercept=(s.rhs[1]-dot(s.gram[1,2:3],slope))/(s.gram[1,1]+intercept_ridge)
    (;intercept,slope,information=S,count=sum(s.count for s in sources))
end
end
