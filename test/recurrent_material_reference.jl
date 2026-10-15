# Preserved pre-hoist Julia law for stencil arithmetic/performance comparisons.
@kernel function recurrent_material_reference_kernel!(out,state,nout,dout,n,d,activity,h,dt,coupling,feedback,writing,erosion,contrast)
    i,j,k=@index(Global,NTuple)
    nx,ny,nz=size(n); pos=(i,j,k); dims=(nx,ny,nz)
    solid=d[i,j,k]; food=n[i,j,k]
    # Two constitutions meet without an imposed temporal forcing.
    constitution=tanh(8f0*(Float32(i)/Float32(nx)-.5f0))
    memory=solid/(1f0+solid)
    drive=0f0
    for channel in 1:8
        pair=(channel+1)÷2; first=2pair-1; second=2pair
        x=state[i,j,k,first]; y=state[i,j,k,second]
        previous=pair==1 ? 7 : first-2
        clock=pair==1 ? 1f0 : pair==2 ? .43f0 : pair==3 ? .17f0 : .065f0
        lateral=0f0
        value=state[i,j,k,channel]
        for axis in 1:3
            if dims[axis]>1
                for side in (-1,1)
                    q=phase_neighbor(pos,axis,side,dims)
                    lateral+=recurrent_face(solid,d[q...])*(state[q...,channel]-value)/h[axis]^2
                end
            end
        end
        odd=isodd(channel)
        # Rotational recurrent prior, coupled across timescales. Mineral changes
        # the threshold and cross-timescale coupling, not just its appearance.
        recurrent=odd ? 1.55f0*x-(1.9f0+contrast*constitution)*y :
            (1.9f0-contrast*constitution)*x+1.55f0*y
        cross=.28f0*state[i,j,k,previous]*(1f0+feedback*memory)
        bias=odd ? -.5f0*feedback*memory+.06f0*tanh(food) : .25f0*feedback*memory
        target=tanh(recurrent+cross+bias)
        updated=value+dt*(clock*(target-value)+coupling*clock*lateral)
        out[i,j,k,channel]=updated
        channel==1 && (drive=updated)
    end
    # Opposite halves of the fast cycle build and erase. Both exchanges are
    # funded locally, with equal/opposite resource entries; no source term.
    growth=min(food,dt*writing*max(0f0,drive-.15f0)*food)
    loss=min(solid,dt*erosion*max(0f0,-drive-.1f0)*solid)
    nf=0f0
    for axis in 1:3
        if dims[axis]>1
            for side in (-1,1)
                q=phase_neighbor(pos,axis,side,dims)
                nf+=.025f0*recurrent_face(solid,d[q...])*(n[q...]-food)/h[axis]^2
            end
        end
    end
    nout[i,j,k]=food-growth+loss+dt*nf
    dout[i,j,k]=solid+growth-loss
    activity[i,j,k]=drive
end
