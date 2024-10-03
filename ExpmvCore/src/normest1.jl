_sgn(x::Number) = ifelse(iszero(x),one(x),sign(x))
function normest1pow(A::StaticMatrix{N,N,T},P::Integer) where {N,T}
	maxit = 20
	P > 0 || throw(DomainError(P,"power must be a positive integer"))
	w = T(1/√N)
	wp =  w
	wm = -w
	if isodd(N)
		Nm = N÷2
		Np = N-Nm
		δw = w/N
		wp -= δw
		wm -= δw
		n = √(Np*wp^2+Nm*wm^2)
		wp /= n
		wm /= n
	end
	X = SMatrix{N,2,T}(ifelse(isone(n),w,ifelse(isodd(m),wp,wm)) for m=1:N,n=1:2)
	est = zero(real(T))
	oldest = zero(real(T))
	indhist = zero(MVector{N,Bool})
	S = zero(SVector{N,T})
	for k ∈ 1:maxit+1
		Y = X
		for _ ∈    1:P
			Y = A*Y
		end
		(est,j) = findmax(sum(abs,Y;dims=1))
		if k ≥ 2 && est ≤ oldest
			est = oldest
			break
		end
		oldest = est
		if k > maxit
			break
		end
		S = map(_sgn,Y)
		# we forego parallel column detection for now
		Z = S
		for _ ∈ 1:P
			Z = A'*Z
		end
		h = maximum(abs,Z;dims=2)[:]
		if k ≥ 2 && maximum(h) == h[j]
			break
		end
		p = sortperm(h;order=Base.Order.Reverse)
		indhist[p[1]] && indhist[p[2]] && break
		count(!,indhist) < 2 && break
		i1 = findfirst(!,indhist)
		i2 = findnext(!,indhist,i1+1)
		indhist[i1] = true
		indhist[i2] = true
		X = SMatrix{N,2,T}(ifelse(i == b,one(T),zero(T)) for i ∈ 1:N,b ∈ (i1,i2))
	end
	return est
end