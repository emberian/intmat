using Interior, JSON3
function mode_check(n=32;steps=160)
 w=World(Grid(n,n,n));dt=.006;p=AcousticField(w;speed=.44,material_feedback=0,damping=0,source_amplitude=0,activity_gain=0)
 h=2/n;lambda=-4/h^2*(sin(pi*2/n)^2+sin(pi*3/n)^2+sin(pi/n)^2);theta=acos(1+.5dt^2*.44^2*lambda)
 mode=[Float32(cos(2pi*(2*(i-.5)+3*(j-.5)+(k-.5))/n)) for i=1:n,j=1:n,k=1:n]
 p.pressure.=mode;p.previous.=cos(theta).*mode
 for k=1:steps;Interior.advance!(w,p,dt);w.time+=dt;end
 err=maximum(abs.(p.pressure.-Float32(cos(steps*theta)).*mode));err<1e-4||error("Fourier discrete dispersion error $err")
 (n=n,steps=steps,max_error=err,discrete_phase=theta,cfl=dt/Interior.stability_limit(w,p))
end
println(JSON3.write(mode_check()))
