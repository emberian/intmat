using Test
include("../studies/material_memory_transfer.jl")
const R=MaterialMemoryTransfer
row(id,step,u;generation=0)=(;id,step,generation,observation=fill(u,3),native_fast=fill(u,3),native_slow=fill(u,3))
@testset "causal normalization and resident identity" begin
    first=[row(1,0,.6),row(1,10,.65),row(2,13,-.5;generation=1),row(1,20,.7),row(2,20,-.4;generation=1)]
    a=R.encode(first);b=R.encode(vcat(first,[row(1,30,-.9),row(2,30,.9;generation=1)]))
    @test a.adaptive==b.adaptive[1:5,:]
    @test a.adaptive[1,:]==zeros(6)
    @test a.adaptive[3,:]==zeros(6) # newborn does not borrow parent's history
    @test a.births==1
    @test a.raw[3,:]==fill(-.5,3) # absolute scarcity signal remains available
    @test_throws ErrorException R.encode([row(1,10,.2),row(1,10,.3)])
    @test_throws ErrorException R.encode([row(1,10,.2),row(1,0,.3)])
end
@testset "future targets respect clocks and IDs" begin
    rows=[row(id,step,.1id;generation=id-1) for step in 0:10:1590 for id in 1:2 if id==1 || step>=900]
    p=R.indices(rows)
    @test !isempty(p.train) && !isempty(p.test)
    @test all(rows[j].id==rows[i].id && rows[j].step-rows[i].step==100 for (i,j) in vcat(p.train,p.test))
    @test all(100<=rows[i].step<700 && rows[j].step<800 for (i,j) in p.train)
    @test all(1000<=rows[i].step<1500 && rows[j].step<1600 for (i,j) in p.test)
    @test all(rows[i].id==1 for (i,_) in p.train)
    @test any(rows[i].id==2 for (i,_) in p.test)
end
