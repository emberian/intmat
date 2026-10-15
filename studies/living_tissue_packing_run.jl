# The same native field exporter and budgets serve both constitutive laws.
include("living_tissue_turnover_run.jl")
include("living_tissue_packing.jl")
using .LivingTissuePacking
function run_packing(recipe,output;backend_kind="cpu")
    backend=BASE.tissue_backend(backend_kind)
    Base.invokelatest(run_packing_backend,recipe,output,backend,backend_kind)
end
function run_packing_backend(recipe,output,backend,backend_kind)
    run_turnover_backend(recipe,output,backend,backend_kind;
        make_base=packing_base_parameters,make_phase=packing_parameters,
        make_process=PackingTissue,measure=packing_statistics,
        save_checkpoint=packing_checkpoint,
        extra_sources=("living_tissue_packing.jl","living_tissue_packing_run.jl"),
        law="living-tissue-packing-v1")
end
if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (2,3) || error("recipe fresh-output [cpu|amd|metal]")
    run_packing(ARGS[1],ARGS[2];backend_kind=length(ARGS)==3 ? ARGS[3] : "cpu")
end
