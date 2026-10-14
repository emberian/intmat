# Bounded scaling measurements, not a launch of long high-resolution studies.
using Interior, Metal
using Interior: JSON3
Metal.functional() || error("Metal unavailable")
const RESULTS = Any[]
function setup(kind,n,backend)
    extent=kind=="phase" ? (32f0,32f0,32f0) : (2f0,2f0,2f0)
    w=World(Grid(n,n,n;extent);backend)
    if kind=="phase"
        p=PhaseSeparation(w;noise=0,bulk=1,kappa=.7,mobility=1)
        a=Array(p.phi)
        for k in 1:n,j in 1:n,i in 1:n
            x=(i-.5f0)/n; y=(j-.5f0)/n; z=(k-.5f0)/n
            a[i,j,k]=.03f0*(sin(2f0*Float32(pi)*2x)+sin(2f0*Float32(pi)*3y)+sin(2f0*Float32(pi)*2z))
        end
        copyto!(p.phi,a)
        return w,(p,),p
    end
    seed_activity!(w)
    processes=(WaveForcing(),Excitable(),Diffusion(:nutrient,.01;material_feedback=2),Deposition())
    w,processes,nothing
end
function run_trial(kind,n,backend,name;steps=20)
    w,processes,p=setup(kind,n,backend)
    limit=minimum(stability_limit(w,q) for q in processes)
    dt=Float32(.8limit)
    # Compile and settle all kernels, reductions and copies before measurement.
    for _ in 1:2; step!(w,processes,dt); end
    elapsed=@elapsed for _ in 1:steps; step!(w,processes,dt); end
    buffers=16+(kind=="phase" ? 4 : 0)
    seconds_per_step=elapsed/steps
    target_time=kind=="phase" ? 30.0 : 10.0
    total_steps=ceil(Int,target_time/dt)
    bytes_per_scalar=n^3*sizeof(Float32)
    result=(kind=kind,backend=name,resolution=n,extent=collect(w.grid.extent),stable_limit=limit,dt=dt,
        measured_steps=steps,seconds=elapsed,seconds_per_step=seconds_per_step,
        resident_scalar_buffers=buffers,resident_field_bytes=buffers*bytes_per_scalar,
        capture_fields=kind=="phase" ? 6 : 3,capture_bytes=(kind=="phase" ? 6 : 3)*bytes_per_scalar,
        estimated_target_time=target_time,estimated_steps=total_steps,estimated_seconds=seconds_per_step*total_steps,
        estimate_basis="same measured step cost and current-state timestep bound; future nonlinear restrictions and storage overhead excluded")
    push!(RESULTS,result)
    println(JSON3.write(result)); flush(stdout)
    nothing
end
# Same analytic initial field at each physical resolution and on both backends.
a,pa,phasea=setup("phase",12,CPU()); b,pb,phaseb=setup("phase",12,MetalBackend())
dt=Float32(.8min(stability_limit(a,phasea),stability_limit(b,phaseb)))
for _ in 1:10; step!(a,pa,dt); step!(b,pb,dt); end
parity=(phi=maximum(abs.(Array(phasea.phi).-Array(phaseb.phi))),interface=maximum(abs.(Array(phasea.interface).-Array(phaseb.interface))),
    activity=maximum(abs.(Array(a.fields.activity).-Array(b.fields.activity))),nutrient=maximum(abs.(Array(a.fields.nutrient).-Array(b.fields.nutrient))))
maximum(values(parity))<1e-5 || error("Phase GPU discrepancy: $parity")
println("phase CPU/Metal parity: ",parity)
for kind in ("dense","phase"), n in (32,64,128)
    n<=64 && run_trial(kind,n,CPU(),"CPU")
    run_trial(kind,n,MetalBackend(),"Metal")
    GC.gc()
end
out=joinpath(@__DIR__,"..","outputs","bigruns"); mkpath(out)
open(joinpath(out,"scaling.json"),"w") do io
    JSON3.pretty(io,(julia=string(VERSION),device="Apple M2 Max 38 GPU cores",phase_parity=parity,measurements=RESULTS,
        note="Buffers count explicit World and process Float32 arrays, excluding Julia/driver/runtime overhead; time predictions are estimates, not completed runs"))
end
