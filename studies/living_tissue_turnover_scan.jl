# Actual one-site process trajectories to screen turnover clocks before spatial
# experiments. These reductions cannot establish a spatial pattern or a limit cycle.
using Interior,Statistics,SHA
using Interior:JSON3
include("living_tissue_turnover.jl")
using .LivingTissueTurnover

function turnover_scan(output;profile="screen")
    ispath(output)&&error("fresh destination required")
    results=Any[]
    gains=profile=="reopen" ? (2.8f0,3.1f0) : (2.5f0,2.8f0,3.1f0)
    withdrawals=profile=="reopen" ? (1.2f0,1.8f0,2.4f0) : (.65f0,.85f0,1.05f0)
    constructions=profile=="reopen" ? (.2f0,) : (.2f0,.3f0)
    growths=profile=="reopen" ? (1f0,1.4f0) : (.5f0,)
    reserve=profile=="reopen" ? 12f0 : 6f0
    for gain in gains,withdrawal in withdrawals,construction in constructions,growth in growths
        w=tissue_world(Grid(1,1,1;extent=(1,1,1)))
        state=zeros(Float32,1,1,1,7);state[1,1,1,:].=(.6,.26,0,0,.03,0,0)
        b=turnover_base_parameters(neural_gain=gain,inhibition_gain=2.3,inhibition_time=4.5,growth=growth,capacity=profile=="reopen" ? 1.2 : .9)
        q=turnover_parameters(withdrawal=withdrawal,construction=construction,unbuilding=profile=="reopen" ? 1.2 : .9)
        p=TurnoverTissue(w;initial=state,reserve,base_parameters=b,parameters=q)
        initial=turnover_statistics(w,p);xs=Float32[];density=Float32[];matrix=Float32[];tail=Any[]
        occupied=true;opens=0;reforms=0
        for k in 1:10000
            step!(w,p,.04)
            if k%25==0
                x=p.tissue.state[1,1,1,5];solid=p.tissue.state[1,1,1,2]+p.tissue.state[1,1,1,3]
                push!(xs,x);push!(density,solid);push!(matrix,p.tissue.state[1,1,1,3])
                if occupied && solid<.1;opens+=1;occupied=false
                elseif !occupied && solid>.24;reforms+=1;occupied=true;end
                k>8000 && push!(tail,(time=w.time,state=vec(copy(p.tissue.state)),reserve=only(p.reserve)))
            end
        end
        final=turnover_statistics(w,p)
        push!(results,(gain,withdrawal,construction,growth,opens,reforms,min_density=minimum(density),
            max_density=maximum(density),tail_x_std=std(xs[321:end]),tail_density_std=std(density[321:end]),
            tail_matrix_std=std(matrix[321:end]),crossings=count(i->xs[i]*xs[i+1]<0,1:length(xs)-1),
            final,material_residual=final.material-initial.material,tail))
        println(gain," ",withdrawal," ",construction," growth ",growth," opens/reforms ",opens,"/",reforms,
            " densitymin ",minimum(density)," tailstd ",std(xs[321:end]));flush(stdout)
    end
    open(io->JSON3.pretty(io,(model_sha256=LivingTissueTurnover.TURNOVER_SHA,
        base_sha256=LivingTissueTurnover.LivingTissueStudy.MODEL_SHA,dt=.04,steps=10000,profile,reserve,results,
        interpretation="One-site screens of the actual composed finite-budget law. Diffusion is inactive. These do not certify a limit cycle, spatial instability, biological organization or adaptation.")),output,"w")
end
if abspath(PROGRAM_FILE)==@__FILE__;length(ARGS) in (1,2)||error("fresh-output.json [screen|reopen]");turnover_scan(ARGS[1];profile=length(ARGS)==2 ? ARGS[2] : "screen");end
