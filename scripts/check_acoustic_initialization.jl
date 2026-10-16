using Interior, Test
using Interior: JSON3
kind=get(ENV,"WAVE_BACKEND","cpu")
if kind=="metal";using Metal;backend=Metal.MetalBackend();else;backend=Interior.CPU();end
function checks(backend)
 n=16;dt=.006f0;w=World(Grid(n,n,n));p=AcousticField(w;speed=.44,material_feedback=0,damping=0,source_amplitude=0,activity_gain=0)
 mode=[Float32(cos(2pi*(2*(i-.5)+3*(j-.5)+(k-.5))/n)) for i=1:n,j=1:n,k=1:n];saved=copy(mode)
 Interior.seed_acoustic_field!(w,p,mode;dt);@test mode==saved
 lambda=-4/(2/n)^2*(sin(2pi/n)^2+sin(3pi/n)^2+sin(pi/n)^2)
 previous=mode.*Float32(1+.5Float64(dt)^2*Float64(p.speed)^2*lambda)
 e=maximum(abs.(p.previous.-previous));@test e<2e-6
 @test_throws DimensionMismatch Interior.seed_acoustic_field!(w,p,zeros(3,4,5);dt)
 @test_throws ArgumentError Interior.seed_acoustic_field!(w,p,mode;dt=-1)
 @test_throws ArgumentError Interior.seed_acoustic_field!(w,p,mode;dt=1)
 @test_throws ArgumentError Interior.seed_acoustic_field!(w,p,fill(Inf,n,n,n);dt)
 a=World(Grid(n,n,n));ap=AcousticField(a);seed_acoustic!(a,ap;amplitude=.8,radius=.13,center=(.1,.2,-.1),dt)
 expected=zeros(Float32,n,n,n)
 for k=1:n,j=1:n,i=1:n
  distance=sum((((i,j,k)[d]-.5f0)*a.grid.spacing[d]-a.grid.extent[d]/2-(.1,.2,-.1)[d])^2 for d=1:3)
  expected[i,j,k]=.8*exp(-distance/(2*.13^2))
 end
 @test Array(ap.pressure)==expected
 # Compare the general field route with the existing Gaussian convenience route.
 b=World(Grid(n,n,n));bp=AcousticField(b);Interior.seed_acoustic_field!(b,bp,Array(ap.pressure);dt)
 @test ap.pressure==bp.pressure;@test ap.previous==bp.previous
 g=World(Grid(n,n,n);backend);gp=AcousticField(g;speed=.44,material_feedback=0,damping=0,source_amplitude=0,activity_gain=0)
 deviceinput=similar(gp.pressure);copyto!(deviceinput,mode);Interior.seed_acoustic_field!(g,gp,deviceinput;dt)
 @test Array(gp.pressure)==p.pressure;@test Array(gp.previous)==p.previous
 for _=1:20;Interior.advance!(w,p,dt);Interior.advance!(g,gp,dt);end
 deviceerror=maximum(abs.(Array(gp.pressure).-p.pressure));@test deviceerror<1e-5
 (fourier_previous_error=e,device_host_initialization_and_20step_error=deviceerror,input_unchanged=true,gaussian_route_exact=true,validation_cases=4,backend=kind)
end
result=checks(backend);open(io->JSON3.pretty(io,result),ARGS[1],"w");println(JSON3.write(result))
