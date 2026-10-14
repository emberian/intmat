@testset "persistent neural field" begin
    for nz in (1,5)
        w=World(Grid(6,6,nz)); v=World(w.grid)
        p=NeuralField(w;mode=:random,seed=31)
        q=NeuralField(v;mode=:random,seed=31)
        @test p.recurrence==q.recurrence
        @test p.latent1==q.latent1
        dt=min(0.01f0,stability_limit(w,p)/2)
        initial=copy(p.latent1)
        for _ in 1:10
            step!(w,(p,Deposition()),dt); step!(v,(q,Deposition()),dt)
        end
        @test w.fields.activity==v.fields.activity
        @test p.latent1!=initial
        @test w.fields.activity≈p.gain.*p.latent1
        @test material_mass(w)≈prod(w.grid.size)
        @test all(isfinite,p.latent1)
        @test_throws ArgumentError step!(w,p,2stability_limit(w,p))
    end
    # Opposing neural biases change physical reaction while conserving material.
    grow=World(Grid(4,4)); rest=World(grow.grid)
    fill!(grow.fields.deposit,0.3f0); fill!(grow.fields.nutrient,0.7f0)
    fill!(rest.fields.deposit,0.3f0); fill!(rest.fields.nutrient,0.7f0)
    p=NeuralField(grow;recurrence=(0,0,0,0),bias=(2,0),diffusion=0,nutrient_feedback=0,initial_noise=0)
    q=NeuralField(rest;recurrence=(0,0,0,0),bias=(-2,0),diffusion=0,nutrient_feedback=0,initial_noise=0)
    for _ in 1:50
        step!(grow,(p,Deposition()),0.02f0)
        step!(rest,(q,Deposition()),0.02f0)
    end
    @test sum(grow.fields.deposit)>0.3f0*16
    @test sum(rest.fields.deposit)<0.3f0*16
    @test material_mass(grow)≈material_mass(rest)
    @test_throws ArgumentError NeuralField(grow;response_time=0)
    @test_throws ArgumentError NeuralField(grow;mode=:trained)
    wrong=NeuralField(World(Grid(3,3)))
    oldactivity=copy(grow.fields.activity)
    @test_throws ArgumentError step!(grow,(WaveForcing(),wrong),0.01f0)
    @test grow.fields.activity==oldactivity
    noisy=NeuralField(grow;initial_noise=100)
    @test maximum(abs,noisy.latent1)<=1
    @test maximum(abs,noisy.latent2)<=1
    # Coupling is an input to neural dynamics, not a display-only feature.
    sensed=World(Grid(4,4)); control=World(sensed.grid)
    fill!(sensed.fields.deposit,1f0); fill!(sensed.fields.pressure,1f0)
    a=NeuralField(sensed;initial_noise=0); b=NeuralField(control;initial_noise=0)
    advance!(sensed,a,0.01f0); advance!(control,b,0.01f0)
    @test sensed.fields.activity!=control.fields.activity
end
