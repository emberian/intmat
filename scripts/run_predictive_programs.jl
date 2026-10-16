#!/usr/bin/env julia
include("../studies/native_predictive_programs.jl")
include("cell_program_observatory.jl")
const P=NativePredictivePrograms
length(ARGS)==2||error("PROFILE FRESH_OUTPUT")
profile,out=abspath.(ARGS);ispath(out)&&error("fresh output required");mkpath(out)
reports=Dict{String,Any}()
for (name,predictive,forced) in (("predictive",true,-1),("persistence",false,-1),("balanced",false,3))
    f=P.N.fixture(profile);w,p=f.w,f.p;h=P.sidecar(p;predictive,forced)
    origin=P.N.C.A.audit(p)
    metadata=(;law="native-predictive-programs-v1",arm=name,options=h.program_options,recipe=f.cfg,
        source=P.Interior.DEVICE_CELLULAR_SOURCE_PIN,habitat_source=P.Interior.DEVICE_HABITAT_SOURCE_PIN,
        programs=[(;id=0,name="Harvest",description="No construction; inspect local state every .2"),(;id=1,name="Build A",description=".4 build at investment.3/allocation.9, then .8 harvest"),(;id=2,name="Build B",description=".4 build at investment.3/allocation.1, then .8 harvest"),(;id=3,name="Balanced build",description=".4 build at investment.3/allocation.5, then .8 harvest")],dt=.02,extent=f.cfg["grid"]["extent"],
        phases=Dict("0"=>"Harvest / reconsider","1"=>"Build dwell","2"=>"Harvest tail"),
        scope="authored finite program scaffold; learned local prediction; chemistry/capture heuristic is not a counterfactual oracle")
    recorder=CellProgramObservatory.Recorder(metadata)
    for step in 1:400
        P.L.advance!(w,p,h,.02;dt=.02)
        step%10==1 && CellProgramObservatory.capture!(recorder,P.observation(w,p,h))
    end
    frame=P.N.frame(w,p,h,origin)
    reports[name]=(;final=P.observation(w,p,h),ledger=frame.ledger,events=P.N.birth_events(h),funded=frame.funded_secretion,
        work=sum(h.model.work),updates=sum(h.model.updates),weights=Array(h.model.weights))
    CellProgramObservatory.write_observatory(recorder,joinpath(out,name*"-observatory.json"))
    open(io->P.JSON3.write(io,reports[name]),joinpath(out,name*"-report.json"),"w")
    println((;name,population=length(frame.cells),work=sum(h.model.work),funded=frame.funded_secretion));flush(stdout)
end
open(io->P.JSON3.write(io,reports),joinpath(out,"report.json"),"w")
