"""Shared trainable tanh RNN, batched private state and exact truncated BPTT.

No objective, physical actuator, reward, identity policy or environment derivative
is embedded here. Callers supply observations, continuation masks and head VJPs.
Native buffers stay on the KA backend. Optimizer updates are external training,
not a claim that a physical resident paid for private parameter plasticity.
"""
module SharedRecurrentRuntime
using KernelAbstractions,LinearAlgebra,Random
export initialize,to_backend,allocate,step!,forward!,backward!,optimizer,update!,check_status!

function initialize(rng,D::Integer,H::Integer,O::Integer;T=Float32,recurrent_gain=.95)
    D>0 && H>0 && O>0 && T in (Float32,Float64) || throw(ArgumentError("positive dimensions and Float32/64 required"))
    gain=T(recurrent_gain)
    isfinite(gain) && gain>=0 || throw(ArgumentError("finite nonnegative representable recurrent gain required"))
    (Wx=randn(rng,T,H,D)/sqrt(T(D)),Wh=gain*Matrix(qr(randn(rng,T,H,H)).Q),
        b=zeros(T,H),Wy=randn(rng,T,O,H)/sqrt(T(H)),by=zeros(T,O))
end
function to_backend(backend,theta)
    map(theta) do a
        b=KernelAbstractions.zeros(backend,eltype(a),size(a)...);copyto!(b,a);b
    end
end
function dimensions(theta)
    H,D=size(theta.Wx);O=size(theta.Wy,1)
    size(theta.Wh)==(H,H) && size(theta.b)==(H,) && size(theta.Wy)==(O,H) && size(theta.by)==(O,) || throw(DimensionMismatch("recurrent parameter layout"))
    all(a->eltype(a)===eltype(theta.Wx)&&KernelAbstractions.get_backend(a)==KernelAbstractions.get_backend(theta.Wx),values(theta)) || throw(ArgumentError("parameter dtype/backend differ"))
    D,H,O
end
function compatible_arrays(theta,arrays)
    T=eltype(theta.Wx);backend=KernelAbstractions.get_backend(theta.Wx)
    all(a->eltype(a)===T&&KernelAbstractions.get_backend(a)==backend,arrays) || throw(ArgumentError("recurrent buffer dtype/backend differ"))
    nothing
end
function validate_gradient(theta,g)
    keys(g)==keys(theta) && map(size,g)==map(size,theta) || throw(DimensionMismatch("parameter gradient layout"))
    compatible_arrays(theta,values(g))
end
function allocate(theta,B::Integer,Tsteps::Integer)
    D,H,O=dimensions(theta);B>0 && Tsteps>0 || throw(ArgumentError("positive batch and horizon required"))
    new(d...)=similar(theta.Wx,eltype(theta.Wx),d...)
    P=sum(length,values(theta));blocks=cld(B*Tsteps,128)
    (h=new(H,B,Tsteps),y=new(O,B,Tsteps),da=new(H,B,Tsteps),dh0=new(H,B),g=map(a->similar(a),theta),
        partials=new(P,blocks),dims=(D,H,O,Int(B),Int(Tsteps)))
end
@inline function parameter(a,p,D,H,O)
    if p<=H*D;a.Wx[p]
    elseif p<=H*D+H*H;a.Wh[p-H*D]
    elseif p<=H*D+H*H+H;a.b[p-H*D-H*H]
    elseif p<=H*D+H*H+H+O*H;a.Wy[p-H*D-H*H-H]
    else;a.by[p-H*D-H*H-H-O*H]
    end
end
@inline function setparameter!(a,p,v,D,H,O)
    if p<=H*D;a.Wx[p]=v
    elseif p<=H*D+H*H;a.Wh[p-H*D]=v
    elseif p<=H*D+H*H+H;a.b[p-H*D-H*H]=v
    elseif p<=H*D+H*H+H+O*H;a.Wy[p-H*D-H*H-H]=v
    else;a.by[p-H*D-H*H-H-O*H]=v
    end
end
@kernel function recurrent_step!(hn,theta,x,hprev,continuation)
    r,b=@index(Global,NTuple)
    value=theta.b[r]
    for k in 1:size(x,1);value+=theta.Wx[r,k]*x[k,b];end
    if continuation[b]!=0
        for k in 1:size(hprev,1);value+=theta.Wh[r,k]*hprev[k,b];end
    end
    hn[r,b]=tanh(value)
end
@kernel function head_step!(y,theta,h)
    o,b=@index(Global,NTuple);value=theta.by[o]
    for k in 1:size(h,1);value+=theta.Wy[o,k]*h[k,b];end
    y[o,b]=value
end
function step!(hn,y,theta,x,hprev,continuation)
    D,H,O=dimensions(theta);B=size(x,2)
    size(x)==(D,B) && size(hprev)==size(hn)==(H,B) && size(y)==(O,B) && size(continuation)==(B,) || throw(DimensionMismatch("one-step recurrent buffers"))
    compatible_arrays(theta,(hn,y,x,hprev))
    KernelAbstractions.get_backend(continuation)==KernelAbstractions.get_backend(theta.Wx) && eltype(continuation)<:Real || throw(ArgumentError("continuation mask backend/type differ"))
    Base.mightalias(hn,hprev) && throw(ArgumentError("next hidden buffer must not alias previous hidden"))
    backend=KernelAbstractions.get_backend(theta.Wx)
    recurrent_step!(backend,(32,1))(hn,theta,x,hprev,continuation;ndrange=(H,B))
    head_step!(backend,(32,1))(y,theta,hn;ndrange=(O,B))
    y
end
@kernel function recurrent_block_step!(h,theta,x,h0,continuation,t)
    r,b=@index(Global,NTuple);value=theta.b[r]
    for k in 1:size(x,1);value+=theta.Wx[r,k]*x[k,b,t];end
    if continuation[b,t]!=0
        for k in 1:size(h,1);value+=theta.Wh[r,k]*(t==1 ? h0[k,b] : h[k,b,t-1]);end
    end
    h[r,b,t]=tanh(value)
end
@kernel function head_block!(y,theta,h)
    o,b,t=@index(Global,NTuple);value=theta.by[o]
    for k in 1:size(h,1);value+=theta.Wy[o,k]*h[k,b,t];end
    y[o,b,t]=value
end
function validate(c,theta,x,h0,continuation)
    D,H,O=dimensions(theta);_,_,_,B,Tsteps=c.dims
    c.dims[1:3]==(D,H,O) && size(x)==(D,B,Tsteps) && size(h0)==(H,B) && size(continuation)==(B,Tsteps) || throw(DimensionMismatch("recurrent block inputs"))
    size(c.h)==size(c.da)==(H,B,Tsteps) && size(c.y)==(O,B,Tsteps) && size(c.dh0)==(H,B) || throw(DimensionMismatch("recurrent cache"))
    size(c.partials)==(sum(length,values(theta)),cld(B*Tsteps,128)) || throw(DimensionMismatch("gradient partial workspace"))
    validate_gradient(theta,c.g)
    compatible_arrays(theta,(x,h0,c.h,c.y,c.da,c.dh0,c.partials))
    KernelAbstractions.get_backend(continuation)==KernelAbstractions.get_backend(theta.Wx) && eltype(continuation)<:Real || throw(ArgumentError("continuation mask backend/type differ"))
    Base.mightalias(c.h,h0) && throw(ArgumentError("block initial hidden buffer must not alias trajectory"))
    nothing
end
function forward!(c,theta,x,h0,continuation)
    validate(c,theta,x,h0,continuation);D,H,O,B,Tsteps=c.dims;backend=KernelAbstractions.get_backend(theta.Wx)
    for t in 1:Tsteps
        recurrent_block_step!(backend,(32,1))(c.h,theta,x,h0,continuation,t;ndrange=(H,B))
    end
    head_block!(backend,(32,1,1))(c.y,theta,c.h;ndrange=(O,B,Tsteps))
    c.y
end
@kernel function hidden_vjp!(da,theta,h,dy,continuation,t)
    r,b=@index(Global,NTuple);value=zero(eltype(da))
    for o in 1:size(dy,1);value+=theta.Wy[o,r]*dy[o,b,t];end
    if t<size(h,3) && continuation[b,t+1]!=0
        for k in 1:size(h,1);value+=theta.Wh[k,r]*da[k,b,t+1];end
    end
    da[r,b,t]=value*(one(value)-h[r,b,t]^2)
end
@kernel function initial_vjp!(dh0,theta,da,continuation)
    r,b=@index(Global,NTuple);value=zero(eltype(dh0))
    if continuation[b,1]!=0
        for k in 1:size(da,1);value+=theta.Wh[k,r]*da[k,b,1];end
    end
    dh0[r,b]=value
end
# Each partial visits a contiguous chunk of 128 (batch,time) samples in fixed
# batch-fastest order. Final reduction visits partials in increasing order.
# No floating atomic, host download, or gradient coupling across independent h.
@kernel function gradient_partials!(partials,da,dy,h,x,h0,continuation,D,H,O)
    p,block=@index(Global,NTuple);value=zero(eltype(partials));B=size(h,2);N=B*size(h,3)
    for sample in (block-1)*128+1:min(block*128,N)
        b=rem(sample-1,B)+1;t=div(sample-1,B)+1
        if p<=H*D
            row=rem(p-1,H)+1;k=div(p-1,H)+1;value+=da[row,b,t]*x[k,b,t]
        elseif p<=H*D+H*H
            k0=p-H*D;row=rem(k0-1,H)+1;k=div(k0-1,H)+1
            if continuation[b,t]!=0;value+=da[row,b,t]*(t==1 ? h0[k,b] : h[k,b,t-1]);end
        elseif p<=H*D+H*H+H
            value+=da[p-H*D-H*H,b,t]
        elseif p<=H*D+H*H+H+O*H
            k0=p-H*D-H*H-H;row=rem(k0-1,O)+1;k=div(k0-1,O)+1;value+=dy[row,b,t]*h[k,b,t]
        else
            value+=dy[p-H*D-H*H-H-O*H,b,t]
        end
    end
    partials[p,block]=value
end
@kernel function gradient_finish!(g,partials,D,H,O)
    p=@index(Global,Linear);value=zero(eltype(partials))
    for block in 1:size(partials,2);value+=partials[p,block];end
    setparameter!(g,p,value,D,H,O)
end
function backward!(c,theta,x,h0,continuation,dy)
    validate(c,theta,x,h0,continuation);size(dy)==size(c.y)||throw(DimensionMismatch("head VJP"));compatible_arrays(theta,(dy,))
    D,H,O,B,Tsteps=c.dims;backend=KernelAbstractions.get_backend(theta.Wx)
    for t in Tsteps:-1:1
        hidden_vjp!(backend,(32,1))(c.da,theta,c.h,dy,continuation,t;ndrange=(H,B))
    end
    initial_vjp!(backend,(32,1))(c.dh0,theta,c.da,continuation;ndrange=(H,B))
    gradient_partials!(backend,(32,1))(c.partials,c.da,dy,c.h,x,h0,continuation,D,H,O;ndrange=size(c.partials))
    gradient_finish!(backend,128)(c.g,c.partials,D,H,O;ndrange=size(c.partials,1))
    c.g
end
mutable struct Optimizer{A,N,S,T}
    m::A
    v::A
    norm::N
    status::S
    step::Int
    kind::Symbol
    rate::T
    beta1::T
    beta2::T
    epsilon::T
    clip::T
end
function optimizer(theta;kind=:adam,rate=.001,beta1=.9,beta2=.999,epsilon=1e-8,clip=1.)
    kind in (:adam,:sgd)||throw(ArgumentError("optimizer kind must be :adam or :sgd"))
    dimensions(theta);T=eltype(theta.Wx);backend=KernelAbstractions.get_backend(theta.Wx)
    r,b1,b2,eps0,cl=T.((rate,beta1,beta2,epsilon,clip))
    all(isfinite,(r,b1,b2,eps0,cl)) && r>0 && 0<=b1<1 && 0<=b2<1 && eps0>0 && cl>0 || throw(ArgumentError("optimizer options must remain finite and valid in parameter precision"))
    moments()=map(a->fill!(similar(a),zero(T)),theta)
    Optimizer(moments(),moments(),KernelAbstractions.zeros(backend,T,2),KernelAbstractions.zeros(backend,Int32,1),0,kind,r,b1,b2,eps0,cl)
end
# Scaled sum of squares avoids overflowing g² before clipping a finite VJP.
@kernel function gradient_norm!(norm,status,g,D,H,O,P,clip)
    if status[1]==0
        scale=zero(eltype(norm));ssq=one(eltype(norm));good=true
        for p in 1:P
            value=abs(parameter(g,p,D,H,O));good&=isfinite(value)
            if isfinite(value) && value!=0
                if scale<value;ssq=one(ssq)+ssq*(scale/value)^2;scale=value
                else;ssq+=(value/scale)^2;end
            end
        end
        if !good
            status[1]=Int32(1)
        else
            root=sqrt(ssq);norm[1]=scale*root
            norm[2]=scale>0 && scale>clip/root ? (clip/scale)/root : one(scale)
        end
    end
end
@kernel function optimizer_step!(theta,g,m,v,norm,status,D,H,O,rate,beta1,beta2,epsilon,clip,correction1,correction2,adam)
    p=@index(Global,Linear)
    if status[1]==0
        grad=parameter(g,p,D,H,O)*norm[2]
        if adam
            nextm=beta1*parameter(m,p,D,H,O)+(one(beta1)-beta1)*grad
            nextv=beta2*parameter(v,p,D,H,O)+(one(beta2)-beta2)*grad^2
            change=rate*(nextm/correction1)/(sqrt(nextv/correction2)+epsilon)
            setparameter!(m,p,nextm,D,H,O);setparameter!(v,p,nextv,D,H,O)
        else
            change=rate*grad
        end
        setparameter!(theta,p,parameter(theta,p,D,H,O)-change,D,H,O)
    end
end
function update!(theta,g,opt::Optimizer)
    D,H,O=dimensions(theta);P=sum(length,values(theta));backend=KernelAbstractions.get_backend(theta.Wx)
    validate_gradient(theta,g);validate_gradient(theta,opt.m);validate_gradient(theta,opt.v)
    size(opt.norm)==(2,) && size(opt.status)==(1,) || throw(DimensionMismatch("optimizer scalar buffers"))
    compatible_arrays(theta,(opt.norm,))
    KernelAbstractions.get_backend(opt.status)==backend && eltype(opt.status)===Int32 || throw(ArgumentError("optimizer status backend/type differ"))
    opt.step+=1
    gradient_norm!(backend,1)(opt.norm,opt.status,g,D,H,O,P,opt.clip;ndrange=1)
    optimizer_step!(backend,128)(theta,g,opt.m,opt.v,opt.norm,opt.status,D,H,O,opt.rate,opt.beta1,opt.beta2,opt.epsilon,opt.clip,
        one(opt.beta1)-opt.beta1^opt.step,one(opt.beta2)-opt.beta2^opt.step,opt.kind===:adam;ndrange=P)
    theta
end
function check_status!(opt::Optimizer)
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(opt.status))
    only(Array(opt.status))==0||throw(ArgumentError("recurrent optimizer halted: nonfinite gradient"))
    opt
end
end
