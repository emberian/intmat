include(joinpath(@__DIR__,"../studies/voltage_between_pieces.jl"))
using Test
function check()
 w=World(Grid(4,4,4));p=Interior.BondNetwork(w;size=(5,4,3),seed=7);n=length(p.node_ids);g=[(.3+.7sin(e)^2)/p.rest_length[e] for e in eachindex(p.bond_ids)];isolated=LinearIndices((5,4,3))[3,2,2]
 for e in eachindex(g);isolated in p.bonds[:,e]&&(g[e]=0);end
 v=zeros(n);buffers=ntuple(_->zeros(n),6);pinned=falses(n);solved=voltage!(v,p,g,buffers...,pinned;tol=1e-11)
 @test pinned[isolated];@test solved.floating_components==1
 free=findall(.!pinned);index=Dict(k=>i for (i,k) in enumerate(free));A=zeros(length(free),length(free));b=zeros(length(free))
 for e in eachindex(g)
  i,j=p.bonds[:,e]
  for (a,c) in ((i,j),(j,i))
   haskey(index,a)||continue;r=index[a];A[r,r]+=g[e]
   if haskey(index,c);A[r,index[c]]-=g[e];else;b[r]+=g[e]*(p.grip[c]==1 ? 1. : 0.);end
  end
 end
 reference=A\b;err=maximum(abs.(reference.-v[free]));@test err<1e-8;@test solved.relative_power_defect<1e-8;@test minimum(v)>=-1e-8&&maximum(v)<=1+1e-8
 # No free unknowns: every node is a floating singleton or an electrode.
 boundary=Dict{String,Any}()
 for mode in ("fully_isolated","all_electrode")
  bp=Interior.BondNetwork(w;size=(2,2,2),seed=11);bn=length(bp.node_ids)
  bg=mode=="fully_isolated" ? zeros(length(bp.bond_ids)) : ones(length(bp.bond_ids))
  mode=="fully_isolated"&&(bp.grip.=0)
  bv=fill(.37,bn);bb=ntuple(_->zeros(bn),6);pins=falses(bn)
  bs=voltage!(bv,bp,bg,bb...,pins;tol=1e-11)
  @test all(pins);@test bs.residual==0.;@test bs.iterations==0
  @test bs.relative_power_defect<1e-12
  if mode=="fully_isolated"
   @test bs.floating_components==bn;@test all(iszero,bv)
   @test bs.joule_power==bs.electrode_power==0.
  else
   @test bs.floating_components==0
   @test bv==Float64.((bp.grip.+1)./2)
   @test bs.joule_power>0
  end
  boundary[mode]=bs
 end
 result=(direct_linear_system_max_error=err,solver=solved,isolated_gauge_checked=true,maximum_principle=true,empty_free_set_fixtures=boundary)
 open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
end
check()
