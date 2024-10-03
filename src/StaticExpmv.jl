#= Based on
Al-Mohy, A. H., & Higham, N. J. (2011).
Computing the action of the matrix exponential, 
with an application to exponential integrators.
SIAM journal on scientific computing, 33(2), 488-511.
=#
module StaticExpmv
    import ExpmvCore: tolerance,theta,M_MAX,P_MAX
    import Base: @propagate_inbounds
    import LinearAlgebra: opnorm,norm,I,tr
    using StaticArrays

    @propagate_inbounds function calculate_s(α::T,m::I)::I where {T <: Number,I <: Integer}
        return ceil(I,α/theta(T,m))
    end
    @propagate_inbounds function parameter_search(nA::Number,m::I)::I where {I <: Integer}
        return m*calculate_s(nA,m)
    end
    @propagate_inbounds function parameters(A::StaticMatrix{N,N,T},t::Number)::Tuple{Int,Int} where {N,T}
        1 ≤ N ≤ 50 || throw(DomainError(N,"leading dimension of A must be ≤ 50; larger matrices require Higham's 1-norm estimation algorithm"))
        nA = opnorm(A,1)
        ntA = abs(t)*nA
        iszero(ntA) && return (0,1)
        @inbounds if nA ≤ 4theta(T,M_MAX)*P_MAX*(P_MAX+3)/(M_MAX*1)
            mo = argmin(Base.Fix1(parameter_search,nA),1:M_MAX)
            s = calculate_s(nA,mo)
            return (mo,s)
        else
            Aᵐ = A*A
            pη = abs(t)*√(opnorm(Aᵐ,1))
            (Cmo::Int,mo::Int) = (typemax(Int),1)
            for p ∈ 2:P_MAX
                Aᵐ *= A
                η = abs(t)*opnorm(Aᵐ,1)^(1/(p+1))
                α = max(pη,η)
                pη = η
                (Cmp::Int,mp::Int) = findmin(Base.Fix1(parameter_search,α),p*(p-1)-1:M_MAX)
                if (Cmp,mp) < (Cmo,mo)
                    (Cmo,mo) = (Cmp,mp)
                end
            end
            s = max(Cmo÷mo,1)
            return (mo,s)
        end
    end

    function expmv end
    @propagate_inbounds function expmv(t::Number,A::StaticMatrix{N,N,T},v::SVector{N}) where {N,T}
        T === StaticArrays.arithmetic_closure(T) || throw(ArgumentError("element type of A must equal its own arithmetic closure; see StaticArrays.arithmetic_closure"))
        μ = tr(A)/N
        A -= μ*I
        mo, s = parameters(A,t)
        F = v
        A *= t/s
        η = exp(μ*t/s)
        ϵ = tolerance(T)
        for _ ∈ 1:s
            c₁ = norm(v,Inf)
            for j ∈ 1:mo
                v = 1/j*(A*v)
                F += v
                c₂ = norm(v,Inf)
                c₁+c₂ ≤ ϵ*norm(F,Inf) && break
                c₁ = c₂
            end
            F *= η
            v = F
            all(isfinite,v) || break
        end
        return F
    end
    @propagate_inbounds function expmv(A::StaticMatrix{N,N,T},v::SVector{N}) where {N,T}
        return expmv(1,A,v)
    end
    export expmv
end
