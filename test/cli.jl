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
