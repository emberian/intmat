# Device block means and detached biological observations; no simulation updates.
using Interior.KernelAbstractions: @kernel,@index
using Serialization
@kernel function habitat_coarse!(dst,cells,habitat,raw,g,ncoarse)
    site=@index(Global,Linear)
    pos=Interior.dc_coords(Int32(site),(ncoarse,ncoarse,ncoarse));block=div(g.dims[1],ncoarse)
    # A uniform accumulator type is required when wide food inventories mix
    # with Float32 catalysts/bound fields; tuple widening in the loop is not
    # GPU-compilable. Only the final display observations round to Float32.
    sums=ntuple(_->0.0,Val(9))
    for z in 1:block,y in 1:block,x in 1:block
        fine=ntuple(d->(pos[d]-Int32(1))*block+(d==1 ? x : d==2 ? y : z),3)
        j=Interior.dc_index(fine...,g.dims)
        sums=let old=sums
            ntuple(k->old[k]+(k<=2 ? cells.fields[j,k]/g.volume :
            k<=4 ? raw.fields[j,k-2]/g.volume :
            k<=6 ? reinterpret(Float32,habitat.enzyme_bits[j,k-4])/g.volume :
            k<=8 ? habitat.bound[j,k-6]/g.volume : cells.permeability[j]),Val(9))
        end
    end
    for k in 1:9;dst[site,k]=sums[k]/Float32(block^3);end
end
function resident_bytes(p)
    groups=(p.cells.state,p.cells.environment,p.habitat.arrays,p.habitat.feedstock.arrays)
    sum(sizeof(eltype(a))*length(a) for group in groups for a in values(group))
end
function pilot_frame!(path,w,p,buffer,coarse)
    c=p.cells;h=p.habitat
    Interior.device_cellular_status!(c)
    habitat_coarse!(c.backend,128)(buffer,c.environment,h.arrays,h.feedstock.arrays,c.geometry,Int32(coarse);ndrange=coarse^3)
    fields=Array(buffer);active=findall(==(1),Array(c.state.active));ids=Array(c.state.ids)[active]
    lineages=Array(c.state.lineages)[active];selected=findall(i->ids[i]<=16 || lineages[i]<=4,eachindex(ids))
    snapshot=(time=w.time,completed_steps=c.attempted_steps,extent=c.geometry.extent,
        field_names=(:free_a,:free_b,:raw_a,:raw_b,:enzyme_a,:enzyme_b,:bound_a,:bound_b,:permeability),
        field_units="concentration; permeability dimensionless; each voxel is a block mean",fields=reshape(fields,coarse,coarse,coarse,9),
        ids=ids,parents=Array(c.state.parents)[active],lineages=lineages,generations=Array(c.state.generations)[active],
        position=Array(c.state.position)[active,:],radii=Array(c.state.radii)[active],mass=Array(c.state.mass)[active],
        bond_ids=Array(c.state.bond_ids)[active,:],trace_ids=ids[selected],
        hidden=Array(c.state.hidden)[active[selected],:],actions=Array(c.state.actions)[active[selected],:],
        energy=Array(c.state.energy)[active[selected]],metabolites=Array(c.state.metabolites)[active[selected],:])
    serialize(path,snapshot);filesize(path)
end
