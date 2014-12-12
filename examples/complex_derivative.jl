using DualNumbers
using ComplexDualNumbers

function complex_derivative_ad{T<:Real}(f,z::Complex{T})
	h = rand(Complex{T},)	#choose a random direction. If derivative truly exists, this doesn't matter, aside from numerical error
	return complex_derivative_ad(f,z,h)
end

function complex_derivative_ad{V<:Real,S<:Real}(f,z::Complex{V},h::Complex{S})
	T = promote_type(V,S)
	zprobe = z + Dual(zero(T),one(T)) * h
	d = f(zprobe)
	notdualpart, dualpart = dual_over_complex(d)
	#in the real AD case, h is real and usually unity, so this division is unnecessary
	return dualpart / h
end

#take the derivative at a point in a bunch of directions, and measure how much it differs
function complex_derivative_soft_exists{T<:Real}(f,z::Complex{T})
	n = 1000
	a = zeros(Complex{T},n,)
	for i=1:n
		a[i] = complex_derivative_ad(f,z)
	end

	dispersion = (maximum(real(a)) - minimum(real(a)))^2  + (maximum(imag(a)) - minimum(imag(a)))^2 
	return dispersion
end


f(z) = z^3			#function from Complex |-> Complex
fp(z) = 3*z^2  		#its derivative
z = rand(Complex{Float64},)

abs(complex_derivative_ad(f,z) - fp(z)) #small

complex_derivative_soft_exists(z->z^3,Complex(1.,2.)) #small


#conj(z) is not differentiable
complex_derivative_soft_exists(z->conj(z),Complex(1.,2.)) #large. 
