import LinearAlgebra
using StaticArrays
using StaticExpmv
using Test

@testset "StaticExpmv.jl" begin
    A3 = LinearAlgebra.I+SMatrix{3,3}(randn(3,3)/10)
    v3 = SVector{3}(randn(3))
    @test expmv(A3,v3) ≈ exp(A3)*v3
    A7 = 10LinearAlgebra.I+SMatrix{7,7}(randn(7,7))
    v7 = SVector{7}(randn(7))
    @test expmv(A7,v7) ≈ exp(A7)*v7
    A11 = SDiagonal(SVector{11}(ifelse(iseven(k),1,-1) for k=1:11))+SMatrix{11,11}(randn(11,11))/100
    v11 = SVector{11}(randn(11))
    @test expmv(A11,v11) ≈ exp(A11)*v11
end
