# Small causal fixture: shared inventory and mineral-to-sound feedback.
using Interior, SHA
using Interior: JSON3
function fixture()
 n=24;dt=.004
 w=World(Grid(n,n,n));p1=ThermalWriter(w;path=:stationary,width=.18,heat=30,conversion=4,stop=.8);p2=ThermalWriter(w;path=:stationary,width=.16,heat=24,conversion=6,stop=.8)
 for _ in 1:200;step!(w,(p1,p2),dt);end
 d1=-expm1.(-p1.conversion.*p1.dose);d2=-expm1.(-p2.conversion.*p2.dose);expected=d1.+d2.-d1.*d2
 err=maximum(abs.(w.fields.deposit.-expected));@assert err<2e-6
 suppression=sum(Float64,d1.*d2);@assert suppression>1
 # Same pulse, same medium; remove just the mineral-to-speed connection.
 function sound(feedback)
  q=World(Grid(n,n,n));q.fields.deposit.=w.fields.deposit
  acoustic=AcousticField(q;speed=.4,material_feedback=feedback,source_amplitude=0,activity_gain=0)
  seed_acoustic!(q,acoustic;amplitude=1,radius=.16,center=(0,0,0),dt=dt)
  for _ in 1:180;step!(q,acoustic,dt);end
  Array(acoustic.pressure)
 end
 coupled=sound(14);ablated=sound(0);difference=maximum(abs.(coupled.-ablated));@assert difference>.01
 result=(schema="warm-thread-causal-fixture",shape=(n,n,n),shared_inventory_composition_max_error=err,independent_inventory_excess=suppression,material_sound_ablation_max_difference=difference,material_relative_error=(sum(Float64,w.fields.deposit)+sum(Float64,w.fields.nutrient)-n^3)/n^3,interpretation="At coarse, deliberately broader writer resolution, shared nutrient prevents both histories independently claiming the same material; fixed deposited material measurably changes propagation of an identical pulse. Not a full-resolution convergence result.")
 open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
end
fixture()
