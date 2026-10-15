# Amber archive: field-screened dielectric growth coupled to persistent thermal writing.
# The electrode roots and susceptibility envelope are designed; every grown cell
# is selected using a converged harmonic field and the evolving local material.
include(joinpath(@__DIR__,"branch_solver_prototype.jl"))
include(joinpath(@__DIR__,"branch_metal_multigrid.jl"))
using Metal, SHA, TOML, Statistics
using Interior.KernelAbstractions: @kernel, @index

@kernel function archive_heat!(out,temp,dose,food,solid,h,dt,position,power,width,D,cooling,threshold,conversion)
 i,j,k=@index(Global,NTuple)
 # Reuse Interior's thermal-writing law; the source follows the selected
 # electrical discharge site instead of an externally prescribed toolpath.
 T,exposure=Interior.thermal_exposure(temp,i,j,k,h,(2f0,2f0,2f0),dt,position,width,power,D,cooling,threshold)
 dose[i,j,k]+=exposure
 f=food[i,j,k];flux=f*(-expm1(-conversion*exposure))
 food[i,j,k]=f-flux;solid[i,j,k]+=flux;out[i,j,k]=T
end
@kernel function archive_gather!(v,t,d,f,u,temp,dose,food,indices)
 i=@index(Global,Linear);j=indices[i]
 v[i]=u[j];t[i]=temp[j];d[i]=dose[j];f[i]=food[j]
end
@kernel function archive_deposit!(food,solid,index,amount)
 i=@index(Global,Linear)
 if i==1;food[index]-=amount;solid[index]+=amount;end
end
function root_mask(n)
 m=falses(n,n,n)
 # Three angled electrode stubs, intentionally unequal and at different depths.
 for (x,y,z,dx) in ((-.34,.04,-.35,.10),(.28,.14,.12,-.07),(-.05,-.20,.43,.06))
  for s in range(0,1;length=ceil(Int,.14n))
   xyz=(x+dx*s,y+.035s,z+.10s)
   I=ntuple(d->clamp(round(Int,(xyz[d]+1)*n/2+.5),2,n-1),3);m[I...]=true
  end
 end
 # Fill each stub with face-connected steps (rounding a diagonal may skip faces).
 for k in 2:n-1,j in 2:n-1,i in 2:n-1
  m[i,j,k] || continue
  for (di,dj,dk) in ((1,1,0),(1,0,1),(0,1,1),(-1,1,0),(-1,0,1),(0,-1,1))
   m[i+di,j+dj,k+dk] || continue
   di!=0 ? (m[i+di,j,k]=true) : (m[i,j+dj,k]=true)
  end
 end
 m
end
function frontier_data(p,level,temp,dose,food,upload)
 sites=findall(p.frontier);indices=upload(Int32[LinearIndices(p.potential)[I] for I in sites]);a=ntuple(_->similar(level.u,Float32,length(sites)),4)
 archive_gather!(Interior.KernelAbstractions.get_backend(level.u),256)(a...,level.u,temp,dose,food,indices;ndrange=length(sites))
 sites,map(Array,a)
end
function choose_archive(p,sites,data,c)
 values,temps,doses,foods=data;weights=zeros(Float64,length(sites));n=size(p.potential,1)
 for (q,I) in enumerate(sites)
  foods[q]>=p.deposit_per_cell || continue
  x,y,z=ntuple(d->2(I[d]-.5)/n-1,3)
  # Restrict the active specimen to a chamfered prism without changing the
  # harmonic box boundary. A soft susceptibility margin prevents box grounding.
  envelope=clamp((.78-abs(x)-.12abs(z))*18,0,1)*clamp((.56-abs(y))*18,0,1)*clamp((.80-abs(z))*18,0,1)
  contacts=count(J->p.occupied[J],Interior.electrical_neighbors(I))
  memory=exp(clamp(c["memory_gain"]*doses[q],0,3))
  weights[q]=contacts*max(0.,1-Float64(values[q]))^p.eta*envelope*memory*foods[q]
 end
 total=sum(weights);total>0 || error("No viable frontier")
 target=Interior.electrical_uniform(p.seed,p.events+1)*total;acc=0.;chosen=lastindex(sites)
 for q in eachindex(sites);acc+=weights[q];if acc>target;chosen=q;break;end;end
 sites[chosen],(temperature=temps[chosen],dose=doses[chosen],nutrient=foods[chosen],weight=weights[chosen])
end
function archive_geometry(p,birth)
 sites=findall(p.occupied);ids=Dict(I=>q for (q,I) in enumerate(sites));edges=Tuple{Int,Int}[]
 for (q,I) in enumerate(sites),axis in 1:3
  J=I+CartesianIndex(ntuple(d->d==axis ? 1 : 0,3));haskey(ids,J)&&push!(edges,(q,ids[J]))
 end
 (positions=Tuple.(sites),edges=edges,birth=[get(birth,I,0) for I in sites],shape=size(p.potential))
end
function archive_run(recipe,out)
 c=TOML.parsefile(recipe);n=branch_grid_size(c["resolution"]);events=c["events"];capture=c["capture_every"];dt=Float32(c["dt"])
 # Conservative upper bound covers host/device multigrid, four material fields,
 # thermal scratch, one capture, and solver temporaries (not Blender).
 estimate=Int(64n^3+512*1024^2);estimate<c["memory_budget_mb"]*1024^2 || error("Estimated working set exceeds recipe memory budget")
 ispath(out)&&error("Refusing existing run: $out");mkpath(out);cp(recipe,joinpath(out,"recipe.toml"))
 source=Dict{String,String}()
 for f in ("src/electrical.jl","src/writing.jl","studies/branch_solver_prototype.jl","studies/branch_metal_multigrid.jl","studies/electrothermal_archive.jl")
  bytes=read(joinpath(dirname(@__DIR__),f));source[f]=bytes2hex(sha256(bytes));write(joinpath(out,basename(f)),bytes)
 end
 backend=MetalBackend();upload=Metal.MtlArray;mask=root_mask(n);p=process(n,mask);p.seed=UInt64(c["seed"]);p.eta=c["eta"];p.deposit_per_cell=.25f0;p.tolerance=c["tolerance"]
 levels=device_levels(hierarchy(p.potential,p.occupied),upload)
 temp=Metal.zeros(Float32,n,n,n);scratch=similar(temp);dose=Metal.zeros(Float32,n,n,n);food=Metal.ones(Float32,n,n,n);solid=Metal.zeros(Float32,n,n,n)
 h=ntuple(_->2f0/n,3);D=Float32(c["diffusion"]);cooling=Float32(c["cooling"])
 dt*(cooling+6D/h[1]^2)<.9 || error("Thermal explicit step violates positivity bound")
 for I in findall(mask);idx=LinearIndices(mask)[I];archive_deposit!(backend,1)(food,solid,idx,p.deposit_per_cell;ndrange=1);end
 birth=Dict(I=>0 for I in findall(mask));history=Any[];frames=Any[]
 initialmass=Float64(n)^3
 function capture_state(event;residual=NaN)
  directory=joinpath(out,"frames",lpad(event,5,'0'));mkpath(directory)
  fields=(temperature=Array(temp),dose=Array(dose),solid=Array(solid))
  for (name,a) in pairs(fields)
   path=joinpath(directory,string(name)*".f32")
   if get(c,"compress_captures",false)
    open(pipeline(`gzip -1`,stdout=path*".gz"),"w") do io;write(io,a);end
   else
    write(path,a)
   end
  end
  open(io->JSON3.pretty(io,archive_geometry(p,birth)),joinpath(directory,"geometry.json"),"w")
  open(io->JSON3.pretty(io,(shape=(n,n,n),extent=(2,2,2),time=event*dt,event=event,order="F",dtype="float32",compression=get(c,"compress_captures",false) ? "gzip lossless" : "none",fields=collect(keys(fields)))),joinpath(directory,"volume.json"),"w")
  hostfood=Array(food);mass=sum(Float64,hostfood)+sum(Float64,fields.solid)
  row=(event=event,time=event*dt,occupied=count(p.occupied),temperature_max=maximum(fields.temperature),dose_max=maximum(fields.dose),solid_max=maximum(fields.solid),nutrient_min=minimum(hostfood),material_relative_error=(mass-initialmass)/initialmass,residual=residual,finite=all(isfinite,fields.temperature)&&all(isfinite,fields.dose)&&all(isfinite,hostfood))
  row.finite&&row.nutrient_min>=0&&abs(row.material_relative_error)<1e-6 || error("Material invariant failed")
  push!(history,row);push!(frames,(step=event,time=Float64(event)*dt,directory=relpath(directory,out)));println(JSON3.write(row));flush(stdout)
 end
 solved=device_solve!(levels;tol=.8p.tolerance);capture_state(0;residual=solved.residual)
 elapsed=@elapsed for event in 1:events
  sites,data=frontier_data(p,first(levels),temp,dose,food,upload);I,localstate=choose_archive(p,sites,data,c)
  Interior.electrical_mark!(p,I);p.events+=1;birth[I]=event;device_mark!(levels,I)
  archive_deposit!(backend,1)(food,solid,LinearIndices(mask)[I],p.deposit_per_cell;ndrange=1)
  position=ntuple(d->Float32((I[d]-.5)*2/n-1),3)
  archive_heat!(backend,256)(scratch,temp,dose,food,solid,h,dt,position,Float32(c["heat"]),Float32(c["width"]),D,cooling,Float32(c["threshold"]),Float32(c["conversion"]);ndrange=(n,n,n));temp,scratch=scratch,temp
  solve_seconds=@elapsed solved=device_solve!(levels;tol=.8p.tolerance)
  event%8==0 && (println(JSON3.write((progress_event=event,solve_seconds=solve_seconds,vcycles=solved.cycles)));flush(stdout))
  if event%capture==0 || event==events
   p.potential.=Array(first(levels).u);r=Interior.potential_residual(p);r<=p.tolerance || error("Harmonic residual exceeded tolerance")
   capture_state(event;residual=r)
  end
 end
 write(joinpath(out,"final-nutrient.f32"),Array(food))
 report=(id=c["id"],title=c["title"],description=c["description"],family="interfaces",shape=(n,n,n),extent=(2,2,2),steps=events,time=events*Float64(dt),dt=dt,elapsed_seconds=elapsed,estimated_working_bytes=estimate,source_sha256=source,params=c,history=history,frames=frames,simulation="Quasistatic screened dielectric growth with converged 3D harmonic potential; discharge-centered heat diffusion/cooling and conservative thermal nutrient-to-mineral conversion reuse Interior ThermalWriter exposure law; accumulated dose enhances subsequent growth susceptibility while depleted nutrient inhibits it.",initialization="Three unequal designed conductive electrode stubs and a soft chamfered-prism susceptibility envelope. No pre-grown branches.",optics="Chosen amber chamber, metallic conductor thickness, cold emission for recently created sites, mineral isosurfaces and temperature volume are artistic optical mappings; no physical luminescence or Maxwell optical solve.")
 open(io->JSON3.pretty(io,report),joinpath(out,"run.json"),"w");println("READY ",out)
end
function archive_verify(output)
 function fixture(upload;gain=1.5)
  n=33;dt=.01f0;h=ntuple(_->2f0/n,3);mask=root_mask(n);p=process(n,mask);p.eta=4.8;p.seed=UInt64(7331);p.deposit_per_cell=.25f0
  levels=device_levels(hierarchy(p.potential,p.occupied),upload);backend=Interior.KernelAbstractions.get_backend(first(levels).u)
  temp=upload(zeros(Float32,n,n,n));scratch=similar(temp);dose=upload(zeros(Float32,n,n,n));food=upload(ones(Float32,n,n,n));solid=upload(zeros(Float32,n,n,n))
  for I in findall(mask);archive_deposit!(backend,1)(food,solid,LinearIndices(mask)[I],.25f0;ndrange=1);end
  device_solve!(levels;tol=3e-7,maxcycles=500);trajectory=Tuple{Int,Int,Int}[]
  for event=1:180
   sites,data=frontier_data(p,first(levels),temp,dose,food,upload);I,_=choose_archive(p,sites,data,Dict("memory_gain"=>gain));push!(trajectory,Tuple(I))
   Interior.electrical_mark!(p,I);p.events+=1;device_mark!(levels,I);archive_deposit!(backend,1)(food,solid,LinearIndices(mask)[I],.25f0;ndrange=1)
   position=ntuple(d->Float32((I[d]-.5)*2/n-1),3)
   archive_heat!(backend,256)(scratch,temp,dose,food,solid,h,dt,position,45f0,.07f0,.0032f0,.8f0,.03f0,.4f0;ndrange=(n,n,n));temp,scratch=scratch,temp
   device_solve!(levels;tol=3e-7,maxcycles=500)
  end
  p.potential.=Array(first(levels).u)
  (temperature=Array(temp),dose=Array(dose),nutrient=Array(food),solid=Array(solid),trajectory=trajectory,residual=Interior.potential_residual(p))
 end
 cpu=fixture(identity);gpu=fixture(Metal.MtlArray);ablated=fixture(identity;gain=0.)
 cpu.trajectory==gpu.trajectory || error("CPU/Metal growth choices differ")
 errors=Dict(string(key)=>maximum(abs.(getfield(cpu,key).-getfield(gpu,key))) for key in (:temperature,:dose,:nutrient,:solid))
 maximum(values(errors))<2e-5 || error("CPU/Metal field parity failed")
 balance=(sum(Float64,gpu.nutrient)+sum(Float64,gpu.solid)-33^3)/33^3
 abs(balance)<1e-6&&minimum(gpu.nutrient)>=0 || error("Material conservation failed")
 firstdifference=findfirst(i->cpu.trajectory[i]!=ablated.trajectory[i],eachindex(cpu.trajectory));firstdifference!==nothing || error("Dose feedback has no measured growth effect")
 result=(shape=(33,33,33),events=180,matched_cpu_metal_events=180,cpu_metal_maximum_errors=errors,material_relative_error=balance,min_nutrient=minimum(gpu.nutrient),harmonic_residual=gpu.residual,first_dose_feedback_ablation_difference=firstdifference,thermal_mineral_gain=sum(Float64,gpu.solid)-.25*(count(root_mask(33))+180),interpretation="Matched coarse fixture uses a spatially scaled heat width and diffusion. Demonstrates conservative thermal writing, portable coupled execution, and a causal dose-to-growth effect; not a universal morphology or physical electrical calibration claim.")
 open(io->JSON3.pretty(io,result),output,"w");println(JSON3.write(result))
end
if abspath(PROGRAM_FILE)==@__FILE__
 ARGS[1]=="--verify" ? archive_verify(ARGS[2]) : archive_run(ARGS[1],ARGS[2])
end
