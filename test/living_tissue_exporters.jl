module PackingExporterChecks
using Test,Interior,TOML
using Interior:JSON3
include(joinpath(@__DIR__,"..","studies","living_tissue_packing_run.jl"))

@testset "shared exporter composes both material laws" begin
    mktempdir() do dir
        common=Dict("id"=>"exporter-composition","title"=>"Tiny recorded material",
            "description"=>"A composition check, not an art specimen.",
            "dimensions"=>[4,4,1],"extent"=>[4.,4.,1.],"seed"=>921,
            "noise_model"=>"bandlimited","dt"=>.01,"steps"=>4,"capture_every"=>2,
            "reserve"=>2.,"base"=>Dict{String,Any}(),"phase"=>Dict{String,Any}())
        original=joinpath(dir,"turnover.toml")
        open(io->TOML.print(io,common),original,"w")
        packing=deepcopy(common);packing["base"]=Dict("capacity"=>.9,"growth"=>0.,"biomass_diffusion"=>0.)
        packing["phase"]=Dict("capacity"=>.9,"growth"=>1.2)
        packed=joinpath(dir,"packing.toml");open(io->TOML.print(io,packing),packed,"w")
        a=joinpath(dir,"original");b=joinpath(dir,"packing")
        run_turnover_backend(original,a,CPU(),"cpu")
        run_packing_backend(packed,b,CPU(),"cpu")
        for (root,law) in ((a,"living-tissue-turnover-v1"),(b,"living-tissue-packing-v1"))
            report=JSON3.read(read(joinpath(root,"run.json"),String))
            @test report.law==law
            @test report.frame_steps==[0,2,4]
            @test report.step_budget==20000
            @test report.history[end].finite
            for name in ("precursor","biomass","matrix","spent","excitation","inhibition","memory","reserve")
                @test filesize(joinpath(root,"frames","0000004",name*".f32"))==4*4*4
            end
        end
        report=JSON3.read(read(joinpath(b,"run.json"),String))
        @test report.history[end].capacity_excess<=16eps(Float32)*report.history[end].capacity
        @test haskey(report.source_sha256,"studies/living_tissue_packing.jl")
        wa,pa=restore_turnover!(joinpath(a,"final.jls"))
        wb,pb=restore_packing!(joinpath(b,"final.jls"))
        @test wa.time==wb.time
        @test packing_statistics(wb,pb).material≈turnover_statistics(wa,pa).material atol=1e-4
        extended=deepcopy(packing);extended["step_budget"]=64000
        extended_recipe=joinpath(dir,"extended.toml");open(io->TOML.print(io,extended),extended_recipe,"w")
        run_packing_backend(extended_recipe,joinpath(dir,"extended"),CPU(),"cpu")
        @test JSON3.read(read(joinpath(dir,"extended","run.json"),String)).step_budget==64000
        extended["step_budget"]=100001
        open(io->TOML.print(io,extended),extended_recipe,"w")
        @test_throws ErrorException run_packing_backend(extended_recipe,joinpath(dir,"bad-budget"),CPU(),"cpu")
    end
end
end # module PackingExporterChecks
