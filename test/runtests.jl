using Test

@time @testset "BuoyData.jl" begin
    @time @testset "Test NDBC.jl" begin
        include("ndbc_script.jl")
    end
end
