# Rebuild the old smooth Memory rings as a finite, resonantly rewritten interior.
# Geometry is explicitly seeded; pressure and mineral then alter one another.
include("interfacial_choir.jl")
function reliquary_initial(n;aperture=.85f0,frequency=1f0)
    solid=zeros(Float32,n,n,n); food=similar(solid); host=similar(solid)
    latent=zeros(Float32,n,n,n,8)
    for k=1:n,j=1:n,i=1:n
        x,y,z=ntuple(d->2f0*((i,j,k)[d]-.5f0)/n-1,3)
        support=0f0
        for layer=0:2
            # Unequal eccentric shells leave actual open sightlines into depth.
            xx=x-.055f0*layer; yy=y+.025f0*layer; zz=z+.035f0*layer
            r=sqrt(xx^2+(yy/.86f0)^2+(zz/1.08f0)^2)
            theta=atan(zz,xx); phi=atan(yy,sqrt(xx^2+zz^2))
            radius=.83f0-.205f0*layer+.028f0*sin(7theta+layer)*cos(5phi)+.017f0*cos(11phi-3theta)
            shell=exp(-((r-radius)/.026f0)^2)
            # Ragged aperture plus narrow branching fissures are initial cuts.
            opening=clamp((yy+.40f0-.10f0*sin(3theta+layer)-.09f0*cos(7theta-2layer))*18f0,0f0,1f0)
            veins=abs(sin(7theta+2f0*sin(5phi)+layer)*cos(8phi-3theta))
            fissure=clamp((veins-.055f0)*14f0,0f0,1f0)
            support=max(support,shell*opening*fissure)
        end
        phase=9atan(z,x)+12y+3sin(7x+4z)
        fraction=.49f0+.13f0*sin(phase)+.07f0*cos(11z-8y)
        solid[i,j,k]=support*fraction;food[i,j,k]=support*(1-fraction);host[i,j,k]=support
        for pair=1:4
            a=phase+(pair-1)*.8f0+pair*.7f0*sin(6z+3x)
            latent[i,j,k,2pair-1]=.7f0*sin(a);latent[i,j,k,2pair]=.7f0*cos(a)
        end
    end
    solid,food,latent,host
end
if abspath(PROGRAM_FILE)==(@__FILE__)
    run(ARGS[1],ARGS[2];initializer=reliquary_initial,
        source_files=("src/recurrent_material.jl","src/acoustics.jl","studies/interfacial_choir.jl","studies/resonant_reliquary.jl"),
        initialization_text="Three eccentric, corrugated mineral shells with deliberately cut ragged apertures and branching fissures; finite nutrient confined initially to those shells. Eight phase-textured recurrent channels. Geometry is seeded, not claimed to emerge from noise; subsequent pressure/mineral feedback is simulated.")
end
