using Interior,Test,Random
module FlatSupportOracle
using Interior: DeviceCellular
include(joinpath(@__DIR__,"fixtures","device_support_flat_reference.jl"))
end
@testset "Carry support iterator preserves frozen flat order and bits" begin
    rng=MersenneTwister(21015)
    for counts in ((1,1,1),(1,2,4),(3,5,7),(26,26,26),(26,3,1),(1,26,2))
        C=5;A=maximum(counts);dims=Int32.(counts)
        indices=zeros(Int32,C,A,3);distance2=zeros(Float32,C,A,3)
        radius2=.18f0^2
        for i in 1:C,d in 1:3,k in 1:counts[d]
            indices[i,k,d]=Int32(mod(k+2i-1,counts[d])+1)
            distance2[i,k,d]=i==2 ? 1f0 : .04f0*rand(rng,Float32)
        end
        normalization=zeros(Float32,C);count=zeros(Int32,C);fallback=zeros(Int32,C)
        for i in 2:C
            for z in 1:counts[3],y in 1:counts[2],x in 1:counts[1]
                r2=(distance2[i,x,1]+distance2[i,y,2])+distance2[i,z,3]
                if r2<radius2
                    count[i]+=1;normalization[i]+=(1f0-r2/radius2)^2
                end
            end
            if count[i]==0
                fallback[i]=1;count[i]=1;normalization[i]=1f0
            end
        end
        d=(indices=indices,distance2=distance2,normalization=normalization,fallback=fallback,
            count=count,candidates=Int32.(counts),dims=dims,radius2=radius2)
        for i in 1:C
            old=collect(FlatSupportOracle.device_support_edges(d,i))
            new=collect(Interior.device_support_edges(d,i))
            @test isequal(new,old)
        end
    end
end
