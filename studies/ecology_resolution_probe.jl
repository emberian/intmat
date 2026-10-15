# Narrow inventory/source normalization probe, not a high-resolution ecology run.
using Interior, Statistics
using Interior: JSON3
const ROOT=normpath(joinpath(@__DIR__,".."))
results=Any[]
for n in (32,64)
 scale=(32/n)^3
 w=World(Grid(n,n,n;extent=(24.,24.,24.));seed=117)
 p=ResourcePool(w;initial=(1.2scale,.8scale,0),source_rates=(.035scale,.02scale,0),sink_rates=(0,0,.01scale),source_budgets=(1e6,1e6,0),sink_budgets=(0,0,1e6),diffusion=(.03,.02,.01))
 initial=resource_mass(p)
 for tick=1:5;step!(w,p,.1);end
 stock=p.fields[1][cld(n,2),cld(n,2),cld(n,2)]
 # Current metabolic uptake uses stock mass in a saturation expression.
 current_request=.14stock/(.2+stock)
 normalized_concentration_request=.14(stock/scale)/(.2+stock/scale)
 push!(results,(resolution=n,extent=(24.,24.,24.),cell_mass_scale=scale,initial_mass=initial,final_mass=resource_mass(p),supplied=copy(p.supplied),removed=copy(p.removed),current_local_request=current_request,concentration_normalized_request=normalized_concentration_request,local_cell_stock=stock,stock_after_equivalent_request=max(0,stock-normalized_concentration_request)))
end
analytical128=(resolution=128,inventory_inflation_if_unscaled=64,cell_mass_scale=1/64,uniform_stock=1.2/64,current_local_request=.14*(1.2/64)/(.2+1.2/64),concentration_normalized_request=.14*1.2/(.2+1.2),single_cell_grant_cap=1.2/64,physical_kernel_cells_for_coarse_cell=64)
out=joinpath(ROOT,"outputs","refinement-20261015","ecology");mkpath(out)
open(io->JSON3.pretty(io,(purpose="Resource-only resolution normalization probe; no claim of coupled ecology equivalence",steps=5,dt=.1,actual_tests=results,analytical128=analytical128,limitations=["The128³ case is analytical, not allocated or simulated.","Per-cell stocks and source/sink rates need cell-volume scaling to preserve total budgets.","Current uptake saturation and sensors consume per-cell mass; these need concentration normalization.","Pointwise uptake/deposit/film at finer grids changes physical interaction support; use fixed-physical-radius exchange kernels.","Neighborhood radii and daughter offsets are in cell coordinates; preserve physical scale explicitly.","Diffusion already uses physical grid spacing; phase-interface units need checking separately."])),joinpath(out,"resolution-probe.json"),"w")
println(results);println(analytical128)
