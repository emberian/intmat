# Private scheduling prototype derived from the integrated DeviceCellular queue.
# No physiology/mechanics kernel is copied. Ordinary production APIs are untouched.
function device_graph_enqueue_prototype!(w,p,delta,execution,anabolism,step_index)
    g=p.geometry;q=p.parameters;s=p.state;e=p.environment
    C=Int(g.capacity);V=prod(Int.(g.dims))
    device_anabolic_begin!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_resource_quota!,e,q,delta;ndrange=3)
    device_cellular_launch!(p,dc_resource_exchange!,e,q,g,delta;ndrange=(V,3))
    device_cellular_launch!(p,dc_resource_finish!,e;ndrange=3)
    device_cellular_launch!(p,dc_diffuse!,e,q,g,delta;ndrange=(V,3))
    device_cellular_launch!(p,dc_commit_resources!,e;ndrange=(V,3))
    device_cellular_launch!(p,dc_clear_grid!,e;ndrange=V)
    device_cellular_launch!(p,dc_clear_heads!,e;ndrange=length(e.heads))
    device_cellular_launch!(p,dc_geometry!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_support!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_sense!,s,e,q,g;ndrange=C)
    if execution===:staged
        device_cellular_launch!(p,dc_controller!,s,e,g;ndrange=(C,Int(g.hidden)+CELLULAR_ACTIONS))
    else
        controller=p.ownership.stages.controller
        device_cellular_launch!(p,dc_controller!,controller.state,controller.environment,controller.geometry;ndrange=(C,Int(g.hidden)+CELLULAR_ACTIONS))
    end
    device_cellular_launch!(p,dc_regulate!,s,e,q,g,delta;ndrange=C)
    device_cellular_launch!(p,dc_allocate_sites!,e;ndrange=(V,2))
    device_cellular_launch!(p,dc_gather!,s,e,g;ndrange=C)
    # Keep the junction stage ABI small: whole-state tuples trigger a pathological
    # LLVM EarlyCSE alias walk in the pinned AMD compiler for some 3D specializations.
    junction=(active=s.active,ids=s.ids,bonds=s.bonds,bond_ids=s.bond_ids,mass=s.mass,
        metabolites=s.metabolites,metabolites_next=s.metabolites_next,actions=s.actions,
        exchange_scale=s.exchange_scale,ledger=s.ledger,ledger_low=s.ledger_low)
    junction_status=(status=e.status,);junction_q=(junction_rate=q.junction_rate,);junction_g=(degree=g.degree,)
    device_cellular_launch!(p,dc_exchange_scale!,junction,junction_status,junction_q,junction_g,delta;ndrange=C)
    device_cellular_launch!(p,dc_exchange_commit!,junction,junction_status,junction_q,junction_g,delta;ndrange=C)
    device_cellular_launch!(p,dc_commit_metabolites!,s,e;ndrange=C)
    device_anabolic_prepare!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_metabolism!,s,e,q,g,device_anabolic_argument(anabolism),delta;ndrange=C)
    device_cellular_launch!(p,dc_return_waste!,s,e,g;ndrange=C)
    device_cellular_launch!(p,dc_commit_waste!,e;ndrange=V)
    if execution===:staged
        device_cellular_mechanics!(p,delta,Val(false))
    else
        device_cellular_mechanics!(p,delta,Val(true))
    end
    device_cellular_launch!(p,dc_propose_bonds!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_commit_bonds!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_lifecycle!,s,e,q,g,step_index;ndrange=1)
    device_cellular_launch!(p,dc_divide!,s,e,q,g,step_index;ndrange=C)
    device_anabolic_inherit!(anabolism,w,p,step_index)
    device_cellular_launch!(p,dc_clear_occupancy!,e;ndrange=V)
    device_cellular_launch!(p,dc_support!,s,e,q,g;ndrange=C)
    device_cellular_launch!(p,dc_porosity_scatter!,s,e,g;ndrange=C)
    device_cellular_launch!(p,dc_porosity_commit!,e,q;ndrange=V)
    device_anabolic_finish!(anabolism,w,p,delta)
    device_cellular_launch!(p,dc_validate_complete!,s,e,g;ndrange=max(C,V))
    device_cellular_launch!(p,dc_complete_step!,e;ndrange=1)
    p
end
