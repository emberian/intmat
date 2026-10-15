using Serialization, JSON3, SHA
root=ARGS[1];out=ARGS[2];ispath(out)&&error("fresh output required")
frames=Dict()
for arm in ("immobile_full","mobile_full","immobile_catalytic","mobile_catalytic")
 p=joinpath(root,arm,"frame-03000.jls");x=deserialize(p);n=x.shape[3];k=div(n,2)
 f=Dict{String,Any}("sha256"=>bytes2hex(sha256(read(p))),"time"=>x.time,"shape"=>collect(x.shape),"extent"=>collect(x.extent),"z_plane_pair"=>[k,k+1])
 for field in (:free,:raw,:bound,:enzyme)
  a=getproperty(x,field);s=(Float64.(a[:,:,k,1]).+a[:,:,k,2].+a[:,:,k+1,1].+a[:,:,k+1,2])./(2Float64(x.volume))
  f[string(field)]=[collect(s[:,j]) for j in axes(s,2)]
 end
 frames[arm]=f
end
open(io->JSON3.pretty(io,frames),out,"w")
