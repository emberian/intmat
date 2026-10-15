@testset "physical conservative cell transport" begin
    # Qualified names permit standalone inclusion before the public module is wired.
    stencil=Interior.CellStencil
    sample=Interior.sample_concentrations
    gradients=Interior.sample_concentration_gradients
    allocate=Interior.allocate_cell_resources!
    deposit=Interior.deposit_cell_resource!
    rebuild=Interior.rebuild_cell_stencil!
    positions=[0.01 0.04 1.1; 0.02 0.06 1.2; 0.0 0.0 0.0]
    support_areas=Float64[]
    for n in (16,32,64)
        grid=Grid(n,n;extent=(2,2,9))
        world=World(grid)
        pool=ResourcePool(world;initial=(0,0,0))
        st=stencil(grid,positions,.3)
        @test st.cell_volume==4/n^2 # The inactive z extent is irrelevant.
        for s in 1:3; fill!(pool.fields[s],Float32(s*st.cell_volume)); end
        @test sample(pool,st)≈repeat(reshape([1.,2.,3.],3,1),1,3) rtol=1e-14
        @test all(iszero,gradients(pool,st))
        for cell in 1:3
            edges=st.offsets[cell]:st.offsets[cell+1]-1
            @test sum(st.weights[edges])≈1 atol=2e-15
            @test length(unique(st.sites[edges]))==length(edges)
            @test length(edges)>1
        end
        push!(support_areas,(st.offsets[2]-st.offsets[1])*st.cell_volume)
        @test abs(support_areas[end]-pi*.3^2)<.08
        before=[sum(Float64,a) for a in pool.fields]
        # Strong overlapping demand exercises synchronous competition.
        requests=[4.0 8.0 6.0; 3.0 1.0 9.0; 0.0 0.0 0.0]
        permutation=[3,1,2]
        other=deepcopy(pool)
        shuffled=stencil(grid,positions[:,permutation],.3)
        grants=allocate(pool,st,requests)
        permuted=allocate(other,shuffled,requests[:,permutation])
        @test grants[:,permutation]≈permuted rtol=2e-14 atol=1e-15
        @test all(s->pool.fields[s]≈other.fields[s],1:3)
        @test all(grants.>=0)
        @test all(grants.<=requests)
        @test all(a->minimum(a)>=0,pool.fields)
        for s in 1:3
            actual=before[s]-sum(Float64,pool.fields[s])
            @test sum(grants[s,:])≈actual atol=2e-14
            @test pool.withdrawn[s]≈actual atol=2e-14
        end
        before=sum(Float64,pool.fields[3])
        actual=deposit(pool,:byproduct,st,1,.035)
        @test 0<actual<=.035
        @test actual≈sum(Float64,pool.fields[3])-before atol=1e-14
        @test pool.returned[3]≈actual atol=1e-15
        @test count(>(Float32(3st.cell_volume)),pool.fields[3])==st.offsets[2]-1
    end
    @test maximum(support_areas)/minimum(support_areas)<1.3

    grid=Grid(16,16;extent=(2,2,1)); world=World(grid)
    pool=ResourcePool(world;initial=(0,0,0))
    # Same location and support: demand ratios are preserved under starvation.
    st=stencil(grid,zeros(3,2),.3)
    fill!(pool.fields[1],.01f0)
    grants=allocate(pool,st,[1.0 3.0; 0 0; 0 0])
    @test grants[1,2]≈3grants[1,1] rtol=1e-14
    @test count(iszero,pool.fields[1])==st.offsets[2]-1
    # Periodic images give the same support, including across corners.
    periodic=stencil(grid,[0.0 2.0 -2.0; 0.0 -2.0 4.0; 0.0 0.0 0.0],.3)
    @test sample(pool,periodic)[:,1]≈sample(pool,periodic)[:,2]
    @test sample(pool,periodic)[:,1]≈sample(pool,periodic)[:,3]
    broad=stencil(grid,zeros(3,1),10)
    @test length(broad.sites)==length(pool.fields[1])
    @test length(unique(broad.sites))==length(broad.sites)
    tiny=stencil(grid,zeros(3,1),1e-20)
    @test length(tiny.sites)==1
    @test tiny.weights==[1.0]
    @test_throws ArgumentError stencil(grid,zeros(3,1),0)
    @test_throws ArgumentError stencil(grid,zeros(3,1),Inf)
    @test_throws DimensionMismatch stencil(grid,zeros(2,1),1)
    @test_throws ArgumentError stencil(grid,fill(NaN,3,1),1)
    @test_throws DimensionMismatch allocate(pool,tiny,zeros(2,1))
    @test_throws ArgumentError allocate(pool,tiny,fill(-1.,3,1))
    @test_throws ArgumentError deposit(pool,1,tiny,1,-1)
    @test_throws BoundsError deposit(pool,1,tiny,2,1)
    @test_throws DimensionMismatch sample(ResourcePool(World(Grid(8,8))),tiny)

    # Sub-ulp requests/deposits do not manufacture ledger mass.
    fill!(pool.fields[1],1f0)
    @test allocate(pool,tiny,[1e-12;0;0;;])[1,1]==0
    @test deposit(pool,1,tiny,1,1e-12)==0
    @test deposit(pool,1,tiny,1,.125)==.125
    # Finite huge competing demands must not overflow their proportional sum.
    huge=stencil(grid,zeros(3,2),1e-20)
    pool.fields[1][huge.sites[1]]=1
    @test allocate(pool,huge,[1e308 1e308;0 0;0 0])[1,:]==[.5,.5]
    # Rebuilding removes stale supports while retaining reusable array objects.
    sites=st.sites; weights=st.weights
    rebuild(st,reshape([1.,1.,0.],3,1))
    @test st.sites===sites && st.weights===weights
    @test size(sample(pool,st))==(3,1)
    rebuild(st,zeros(3,0))
    @test isempty(st.sites) && isempty(st.site_heads)
    @test size(sample(pool,st))==(3,0)
    @test size(allocate(pool,st,zeros(3,0)))==(3,0)

    # Gradient units and singleton axes, with a sub-grid support at a known center.
    st=stencil(grid,reshape([.5625,.5625,0.],3,1),1e-4)
    for index in CartesianIndices(grid.size)
        x=(index[1]-.5)*2/16
        pool.fields[1][index]=Float32((2+sin(pi*x))*st.cell_volume)
    end
    g=gradients(pool,st)
    expected=pi*cos(pi*.5625)*sin(pi/8)/(pi/8)
    @test g[1,1,1]≈expected atol=2e-6
    @test all(iszero,g[:,2:3,:])
    grid3=Grid(8,12,16;extent=(2,3,4))
    st3=stencil(grid3,reshape([.5,.75,1.],3,1),.6)
    pool3=ResourcePool(World(grid3);initial=(0,0,0))
    @test st3.cell_volume==.25^3
    fill!(pool3.fields[1],Float32(7st3.cell_volume))
    @test sample(pool3,st3)[1,1]≈7 atol=1e-14
    @test all(iszero,gradients(pool3,st3))
    @test length(unique(st3.sites))==length(st3.sites)
end
