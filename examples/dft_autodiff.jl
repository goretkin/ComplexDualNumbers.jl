using DualNumbers
using ComplexDualNumbers

function dft_nocomplex{T<:Number}(x::Vector{T})
	#T should not have any imaginary part 

	N = length(x)
	
	#store DFT here
	C = zeros(T,N)
	S = zeros(T,N)

	for k=0:(N-1)
		for n=0:(N-1)
			cs = exp(-im*2pi* k*n/N)
			c,s = real(cs), imag(cs)
			@inbounds C[k+1] += c*x[n+1]
			@inbounds S[k+1] += s*x[n+1]
		end
	end
	return C,S 
end

function dft!(CS::Vector,x::Vector)
	N = length(x)
	@assert length(CS) == N
	for k=0:(N-1)
		for n=0:(N-1)
			arg = -im*2pi* k*n/N
			cs = exp(arg)
			@inbounds CS[k+1] += cs*x[n+1]
		end
	end
end

function dft{T<:Real}(x::Vector{T})
	N = length(x)
	CS = zeros(Complex{T},N)
	dft!(CS,x)
	return CS
end

function dft{T<:Complex}(x::Vector{T})
	N = length(x)
	CS = zeros(T,N)
	dft!(CS,x)
	return CS
end

function dft{S<:Real}(x::Vector{Dual{S}})
	N = length(x)
	CS = zeros(ComplexDual{S},N)
	dft!(CS,x)
	return CS
end

function dft{S<:Real}(x::Vector{ComplexDual{S}})
	N = length(x)
	CS = zeros(ComplexDual{S},N)
	dft!(CS,x)
	return CS
end

function idft(x::Vector)
	return (1/length(x)) * conj(dft(conj(x))) #return (1/length(x)) * dft(reverse(x)) should have worked too
end

#make sure our dft produces same result as built-in fft for real case.
r = rand(10,)
C,S = dft_nocomplex(r)
r_dft = C + im*S
r_fft = fft(r)#
@assert  maximum(abs(r_dft - r_fft)) <100eps()

#take gradient of a function R^n -> R at x
function grad{T<:Real}(fun::Function,x::Vector{T})
	N = length(x)
	G = zeros(T,N)
	for i = 1:N
		xd = [Dual(x[j],1.0*(j==i)) for j=1:N]
		G[i] = epsilon(fun(xd))
	end
	return G
end

function spec_nocomplex(x)
	C,S = dft_nocomplex(x)		#requires a non-standard form of the dft
	h = floor(length(x) / 2)
	w = cat(1, [1.0:h ], zeros(length(x)-int(h)) )
	sum( (C.^2 + S.^2) .* w ) #some weighting on different frequencies
end

function spec(x)
	CS = dft(x)					#same interface as built-in fft
	h = floor(length(x) / 2)
	w = cat(1, [1.0:h ], zeros(length(x)-int(h)) )
	sum( (abs2(CS)) .* w ) #some weighting on different frequencies
end

#spec and spec_nocomplex are mathematically identical
@assert maximum(abs(spec(r)-spec_nocomplex(r))) < 100eps()

g_nocomplex = grad(spec_nocomplex,r)	#can do this with Dual numbers only
g =  grad(spec,r)						#need ComplexDual

@assert maximum(abs(g_nocomplex-g)) < 100eps()