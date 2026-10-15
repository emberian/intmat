# Scalar site laws for CPU/AMDGPU habitat kernels. Inputs are validated finite,
# nonnegative Float32 inventories/rates, positive volume, and B <= min(M,C).
# Float64 intermediates are intentional; Metal does not support this first law.
@inline function dh_sorption_target(B::Float32,M::Union{Float32,Float64},C::Float32,
        kon::Float32,koff::Float32,volume::Float32,dt::Float32)
    b,m,c,on,off,v,t=Float64(B),Float64(M),Float64(C),Float64(kon),Float64(koff),Float64(volume),Float64(dt)
    c==0 && return 0f0
    on==0 && return Float32(b*exp(-off*t))
    a=on/v
    d=hypot(a*(m-c),sqrt(off*(2a*(m+c)+off)))
    d==0 && return Float32(c-(c-b)/(1+a*(c-b)*t))
    small=2a*m*c/(a*(m+c)+off+d)
    offset=b-small;decay=exp(-d*t)
    Float32(clamp(small+offset*decay/(1-a*offset*(-expm1(-d*t))/d),0.,min(m,c)))
end

# Actual endpoint changes pay for each transfer. Unreceivable representable
# donor debit belongs to inert material; the caller accounts its potential heat
# exactly once. Ordinary exchange never forces an unrepresentable desorption.
@inline function dh_sorption_exchange(free::Float32,bound::Float32,target::Float32)
    if target>bound
        request=min(Float64(free),Float64(target)-Float64(bound))
        nextfree,taken=dh_take_down(free,request)
        nextbound,gained=dh_add_down(bound,taken)
        return (;free=nextfree,bound=nextbound,inert=taken-gained)
    elseif target<bound
        request=Float64(bound)-Float64(target)
        # Extreme inventory ratios can erase a positive target in subtraction.
        Float64(bound)-request<Float64(target) && (request=prevfloat(request))
        nextfree,received=dh_add_down(free,request)
        nextbound,lost=dh_take_up(bound,received)
        return (;free=nextfree,bound=nextbound,inert=lost-received)
    end
    (;free,bound,inert=0.)
end

# Scaffold loss is compulsory: release what the soluble Float32 endpoint can
# receive, then retain any remaining excess as inert mass, never excess binding.
@inline function dh_sorption_collapse(free::Float32,bound::Float32,capacity::Float32)
    exchanged=dh_sorption_exchange(free,bound,min(bound,capacity))
    nextbound=min(exchanged.bound,capacity)
    (;free=exchanged.free,bound=nextbound,
      inert=exchanged.inert+(Float64(exchanged.bound)-Float64(nextbound)))
end
