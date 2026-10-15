# Bounded CPU reference comparison of the writer's linear thermal Green function.
using Interior
using Interior: JSON3
function simpson(f,a,b,n=256)
    iseven(n) || error("even quadrature count required")
    dx=(b-a)/n
    dx/3*(f(a)+f(b)+sum((isodd(i) ? 4 : 2)*f(a+i*dx) for i in 1:n-1))
end
function gaussian_temperature(point,t;heat=30.,width=.1,D=.0003,cooling=4.,position=s->(0.,0.,0.))
    t==0 && return 0.
    simpson(0.,t) do age
        variance=width^2+2D*age; center=position(t-age)
        heat*(width^2/variance)^1.5*exp(-cooling*age-sum((point[d]-center[d])^2 for d in 1:3)/(2variance))
    end
end
function probe_green()
    rows=Any[]
    for n in (32,64,96)
        w=World(Grid(n,n,n));p=ThermalWriter(w;path=:stationary,width=.1,heat=30,diffusion=.0003,cooling=4,threshold=.05,conversion=18)
        for _ in 1:200;step!(w,p,.001);end
        probes=Any[]
        for offset in (0,round(Int,.1*n/2),round(Int,.2*n/2))
            cell=(n÷2+offset,n÷2,n÷2); point=ntuple(d->(cell[d]-.5)*w.grid.spacing[d]-1,3)
            temperature=gaussian_temperature(point,w.time)
            dose=simpson(t->max(0.,gaussian_temperature(point,t)-.05),0.,w.time,128)
            push!(probes,(;point,temperature_green=temperature,temperature_grid=p.temperature[cell...],dose_green=dose,dose_grid=p.dose[cell...]))
        end
        push!(rows,(;shape=w.grid.size,time=w.time,dt=.001,probes))
    end
    result=(;schema="bounded-thermal-green-comparison",source="stationary Gaussian; infinite-space Green kernel; periodic-image tails negligible here",rows)
    println(JSON3.write(result));result
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    result=probe_green()
    !isempty(ARGS) && open(io->JSON3.pretty(io,result),ARGS[1],"w")
end
