#!/usr/bin/env julia
using Test
include("living_resource_fusion.jl")
const F=LivingResourceFusion
using Interior,KernelAbstractions
@testset "implicit grants preserve current endpoints and ordered tree" begin
    for T in (Float32,Float64),n in (1,7,127,128,129,16641),mode in (:zero,:source,:sink,:mixed)
        fields=[T(mod(i*17+k*3,97))/T(19) for i in 1:n,k in 1:3]
        masks=[Float32(mod(i+k,11))/11 for i in 1:n,k in 1:3]
        e=(fields=fields,source_masks=masks,sink_masks=reverse(masks;dims=1),quota=fill(T(.7),3,2),status=zeros(Int32,1))
        q=(source_rates=mode in (:source,:mixed) ? (.01f0,.03f0,.07f0) : (0f0,0f0,0f0),sink_rates=mode in (:sink,:mixed) ? (.02f0,.04f0,.08f0) : (0f0,0f0,0f0))
        a=deepcopy(e);b=deepcopy(e);wa=Interior.DeviceResourceWorkspace(CPU(),n;inventory_type=T);wb=Interior.DeviceResourceWorkspace(CPU(),n;inventory_type=T)
        Interior.dc_resource_tree_grants!(CPU(),128)(a,q,.02f0,wa.grants;ndrange=(n,3))
        F.fused_first!(CPU(),128)(b,q,.02f0,first(wb.levels).high,first(wb.levels).low;ndrange=(size(first(wb.levels).high,1),3))
        ih=wa.grants;il=nothing
        for (i,l) in enumerate(wa.levels)
            Interior.dc_resource_tree_reduce!(CPU(),128)(l.high,l.low,ih,il,a.status,size(ih,1);ndrange=size(l.high));ih=l.high;il=l.low
            if i>1
                prior=wb.levels[i-1];other=wb.levels[i]
                Interior.dc_resource_tree_reduce!(CPU(),128)(other.high,other.low,prior.high,prior.low,b.status,size(prior.high,1);ndrange=size(other.high))
            end
        end
        @test isequal(a.fields,b.fields)
        @test a.status==b.status==[0]
        for (la,lb) in zip(wa.levels,wb.levels)
            @test isequal(la.high,lb.high)
            @test isequal(la.low,lb.low)
        end
    end
end
