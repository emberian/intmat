using Interior, TOML, JSON3, Printf
base=TOML.parsefile(joinpath(@__DIR__,"choir.toml"))
mkpath("outputs/scout")
variants=[("fast",.15,.18,.5,.4), ("brake",.12,.12,1.5,.5), ("lace",.18,.16,2.5,.65), ("slow",.2,.10,3.,.8), ("front",.08,.12,.25,.18), ("veil",.09,.09,.8,.22)]
summary=[]
for (name,threshold,recovery,feedback,deposition) in variants
 c=deepcopy(base);c["simulation"]["resolution"]=32;c["simulation"]["dimension"]=2;c["particles"]["count"]=0
 c["excitable"]["threshold"]=threshold;c["excitable"]["recovery_rate"]=recovery;c["excitable"]["material_feedback"]=feedback;c["deposition"]["rate"]=deposition
 c["wave"]["activity_gain"]=.04
 scene=build_scene(c);w=scene.world;dir="outputs/scout/$name";mkpath(dir)
 open("$dir/config.toml","w") do io;TOML.print(io,c);end
 for i in 1:1000
  step!(w,scene.processes,.02f0)
  if i%100==0
   for f in (:activity,:deposit);open("$dir/$(i)-$f.f32","w") do io;write(io,Array(getproperty(w.fields,f)));end;end
   push!(summary,(name=name,step=i,time=w.time,deposit_mean=sum(w.fields.deposit)/length(w.fields.deposit),activity_min=minimum(w.fields.activity),activity_max=maximum(w.fields.activity)))
  end
 end
 println(name," done")
end
open("outputs/scout/summary-all.json","w") do io;JSON3.pretty(io,summary);end
