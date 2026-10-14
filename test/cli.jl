@testset "Unified application interface" begin
    root=normpath(joinpath(@__DIR__,".."))
    cli=joinpath(root,"bin","interior.jl")
    for command in (String[],["arena","--help"],["run","--help"],["bigrun","--help"],["study","--help"])
        text=read(`$(Base.julia_cmd()) --project=$root $cli $command`,String)
        @test !isempty(text)
    end
    for folder in ("arena","bigruns"), path in readdir(joinpath(root,"examples",folder);join=true)
        endswith(path,".toml") || continue
        @test islink(path)
        @test isfile(realpath(path))
    end
end

@testset "Material command preserves existing output" begin
    root=normpath(joinpath(@__DIR__,".."))
    cli=joinpath(root,"bin","interior.jl")
    mktempdir() do out
        sentinel=joinpath(out,"frame-0000-deposit.f32")
        write(sentinel,"existing recording")
        # Deliberately invalid recipe: ownership rejection must precede even parse,
        # and therefore precede backend loading and expensive scene allocation.
        missing_profile=joinpath(out,"missing.toml")
        error_output=IOBuffer()
        command=`$(Base.julia_cmd()) --project=$root $cli run $missing_profile $out cpu --steps 1 --no-render`
        process=run(pipeline(ignorestatus(command);stdout=devnull,stderr=error_output))
        @test !success(process)
        @test occursin("output must be a fresh or empty directory",String(take!(error_output)))
        @test read(sentinel,String)=="existing recording"
        @test readdir(out)==[basename(sentinel)]
    end
end
