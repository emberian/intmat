# Seeded finite population used by full-loop probes, not a trained genome.
function device_cellular_population(side,backend;mechanical_substeps=16,bin_width=parse(Float32,get(ENV,"DEVICE_BIN_WIDTH","1")),capacity=parse(Int,get(ENV,"DEVICE_CAPACITY",string(2side^3))))
    n=side^3;extent=1.7side;dims=max(4,2side)
    grid=Grid(dims,dims,dims;extent=(extent,extent,extent))
    wcpu=World(grid,CPU(),NamedTuple(),NamedTuple(),0.,17)
    r=ResourcePool(wcpu;potentials=(4,4,0),initial=(.4,.4,0),diffusion=(.02,.02,.01))
    x=zeros(3,n)
    for (i,xyz) in enumerate(Iterators.product(0:side-1,0:side-1,0:side-1));x[:,i].=1.7.*(collect(xyz).+.5);end
    c=CellularEcology(wcpu,x,r;capacity=capacity,biomass=2,energy=1,mutation=.01,developmental_partition=.02)
    # A declared prior enables adhesion, junctions and division; other recurrent weights remain seeded random.
    h=size(c.hidden,1);c.weights[h+8,end,:].=2;c.weights[h+9,end,:].=2;c.weights[h+10,end,:].=2
    for i in 16:16:n;c.biomass[i]=.1;end
    # Explicit empty World works with both the pinned and current package:
    # this process owns every chemical/material buffer itself.
    w=World(wcpu.grid,backend,NamedTuple(),NamedTuple(),0.,wcpu.seed)
    w,Interior.DeviceCellular(w,c;mechanical_substeps=mechanical_substeps,bin_width=bin_width)
end
