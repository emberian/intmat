# Independent positive-edge enumeration: source power, Joule ledger and flux divergence.
function lens_edge_checks(v,sigma,h,voltage)
 n=size(v,1);flux=zeros(size(v));power=0.;input=0.
 @inbounds for k in 1:n,j in 1:n,i in 1:n,axis in 1:3
  x=i+(axis==1);y=j+(axis==2);z=k+(axis==3)
  max(x,y,z)>n && continue
  a=sigma[i,j,k];b=sigma[x,y,z];g=2a*b/(a+b);drop=v[i,j,k]-v[x,y,z];f=g*drop/h^2
  flux[i,j,k]+=f;flux[x,y,z]-=f;power+=g*drop^2*h
  axis==1 && i==1 && (input+=g*drop*voltage*h)
 end
 residual=maximum(abs,@view flux[2:n-1,:,:])/maximum(sigma)
 (independent_residual=residual,independent_edge_power=power,independent_input_power=input,independent_power_defect=(power-input)/max(input,eps()))
end
