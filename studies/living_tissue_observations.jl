module TissueObservations
export tissue_components,tissue_topology,TissueTurnoverObservation,observe_turnover!
mutable struct TissueTurnoverObservation{B,C,A}
    occupied::B
    reopenings::C
    reformations::C
    minimum::A
    maximum::A
    low::Float32
    high::Float32
end
function TissueTurnoverObservation(density;low=.10f0,high=.24f0)
    0<=low<high || throw(ArgumentError("ordered nonnegative density thresholds required"))
    TissueTurnoverObservation(density.>=low,zeros(Int,size(density)),zeros(Int,size(density)),
        Float32.(density),Float32.(density),Float32(low),Float32(high))
end
function observe_turnover!(p::TissueTurnoverObservation,density)
    size(density)==size(p.occupied) || throw(DimensionMismatch("turnover observation grid"))
    for i in eachindex(density)
        if p.occupied[i] && density[i]<p.low
            p.reopenings[i]+=1;p.occupied[i]=false
        elseif !p.occupied[i] && density[i]>p.high
            # Initially thin sites can form without having reopened. A
            # reformation event requires a previously observed opening.
            p.reopenings[i]>0 && (p.reformations[i]+=1)
            p.occupied[i]=true
        end
        p.minimum[i]=min(p.minimum[i],density[i]);p.maximum[i]=max(p.maximum[i],density[i])
    end
    (reopened_sites=count(>(0),p.reopenings),reformed_sites=count(>(0),p.reformations),
        reopening_events=sum(p.reopenings),reformation_events=sum(p.reformations))
end
function tissue_components(mask)
    nx,ny=size(mask);seen=falses(nx,ny);queue=CartesianIndex{2}[];sizes=Int[]
    for start in CartesianIndices(mask)
        (!mask[start] || seen[start]) && continue
        empty!(queue);push!(queue,start);seen[start]=true;head=1
        while head<=length(queue)
            i,j=Tuple(queue[head]);head+=1
            for neighbor in (CartesianIndex(mod1(i-1,nx),j),CartesianIndex(mod1(i+1,nx),j),
                CartesianIndex(i,mod1(j-1,ny)),CartesianIndex(i,mod1(j+1,ny)))
                if mask[neighbor] && !seen[neighbor];seen[neighbor]=true;push!(queue,neighbor);end
            end
        end
        push!(sizes,length(queue))
    end
    sort!(sizes;rev=true)
end
function tissue_topology(state;threshold=.15f0)
    size(state,3)==1 || throw(ArgumentError("2D topology observation required"))
    mask=state[:,:,1,2].+state[:,:,1,3].>threshold
    components=tissue_components(mask);voids=tissue_components(.!mask)
    (threshold=threshold,occupied_components=length(components),occupied_sizes=components,
        void_components=length(voids),void_sizes=voids)
end
end
