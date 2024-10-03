module ExpmvCore
    using Roots
    import Base: @propagate_inbounds
    tolerance(::Type) = ldexp(1.0,-53)
    tolerance(::Type{Float32}) = ldexp(one(Float32),-24)
    @inline function T(M::Integer,x::Number)
        y = one(x)
        for m ∈ M:-1:1
            y = 1+x/m*y
        end
        y
    end
    T(M::Integer) = Base.Fix1(T,M)
    @inline function h(M::Integer,x::Number)
        log(max(zero(x),exp(-x)*T(M,x)))
    end
    h(M::Integer) = Base.Fix1(h,M)
    h̃(M::Integer,x::Number) = ifelse(isodd(M),-1,1)*h(M,-x)
    function θf((M,ϵ)::Tuple{<:Integer,<:Number},x::Number)
        h̃(M+1,x)/x-ϵ
    end
    # this could probably be improved significantly but
    # is not a priority as we calculate thetas only once
    function calc_thetas(m_max,::Type{T}) where {T <: AbstractFloat}
        thetas = Vector{BigFloat}(undef,m_max)
        (l,u) = BigFloat.((√(eps(T)),one(T)))
        ϵ = BigFloat(tolerance(T))
        for m=1:m_max
            thetas[m] = find_zero(Base.Fix1(θf,(m,ϵ)),(l,u);xatol=big(0.0),xrtol=big(0.0),atol=big(0.0),rtol=big(1e-17),verbose=false)
            l = thetas[m]
            u = l+0.3
        end
        v = T.(thetas)
        pushfirst!(v,eps(T))
    end
    const NUM_THETA = 99
    const P_MAX = 8
    const M_MAX = 55
    const THETAS32 = Tuple(calc_thetas(M_MAX,Float32))
    const THETAS64 = Tuple(calc_thetas(M_MAX,Float64))
    @propagate_inbounds theta(::Type{Float64},m::Integer) = THETAS64[m]
    @propagate_inbounds theta(::Type{Float32},m::Integer) = THETAS32[m]
    @propagate_inbounds theta(::Type{Complex{T}},m::Integer) where {T} = theta(T,m)
    @propagate_inbounds theta(::Type{T},::Integer) where {T} = throw(DomainError(T,"type must be either Float32 or Float64"))
    @propagate_inbounds theta(x::Number,m::Integer) = theta(typeof(x),m)

    export NUM_THETA,P_MAX,M_MAX,theta,tolerance
end