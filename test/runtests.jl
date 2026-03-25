using Test, SafeTestsets

@time @testset verbose=true "BuoyData.jl" begin
    @time @safetestset "Test NDBC.jl" begin
        include("test_ndbc.jl")
    end
end
