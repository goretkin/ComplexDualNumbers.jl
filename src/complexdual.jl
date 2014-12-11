immutable ComplexDual{T<:Real} <: Number
    re::T
    du::T
    im::T
    imdu::T
end

ComplexDual(x::Real, y::Real, z::Real, w::Real) = ComplexDual(promote(x,y,z,w)...)
ComplexDual{T<:Real}(x::T) = ComplexDual(x, zero(T),zero(T),zero(T))
ComplexDual{T}(x::Complex{T}) = ComplexDual(real(x), zero(T),imag(x),zero(T))
ComplexDual{T}(x::Dual{T}) = ComplexDual(real(x), epsilon(x),zero(T),zero(T))

typealias ComplexDual256 ComplexDual{Float64}
typealias ComplexDual128 ComplexDual{Float32}

real(z::ComplexDual) = z.re
epsilon(z::ComplexDual) = z.du
imag(z::ComplexDual) = z.im
imagepsilon(z::ComplexDual)= z.imdu

imagfull(z::ComplexDual)= Dual(imag(z),imagepsilon(z)) #not sure about this
epsfull(z::ComplexDual)= Complex(epsilon(z),imagepsilon(z)) #not sure about this

#eps(z::ComplexDual) = eps(real(z))
#eps{T}(::Type{ComplexDual{T}}) = eps(T)
#one(z::Dual) = dual(one(real(z))) why is this present?
one{T}(::Type{ComplexDual{T}}) = complexdual(one(T))
inf{T}(::Type{ComplexDual{T}}) = complexdual(inf(T))
nan{T}(::Type{ComplexDual{T}}) = nan(T)
isnan(z::ComplexDual) = isnan(real(z))  || isnan(imag(z)) #match DualNumbers in treatment of NaNness of dual part.

convert{T<:Real}(::Type{ComplexDual{T}}, x::Real) =
  ComplexDual(convert(T, x))

convert{T<:Real}(::Type{ComplexDual{T}}, x::Complex) =
  ComplexDual(convert(T, real(x)), zero(T), convert(T, imag(x)), zero(T))

convert{T<:Real}(::Type{ComplexDual{T}}, x::Dual) =
  ComplexDual(convert(T, real(x)), convert(T, epsilon(x)),zero(T), zero(T))


convert{T<:Real}(::Type{ComplexDual{T}}, z::ComplexDual{T}) = z #why?
convert{T<:Real}(::Type{ComplexDual{T}}, z::ComplexDual) =
  ComplexDual(convert(T, real(z)), convert(T, epsilon(z)), convert(T,imag(z)), convert(T,imagepsilon(z)))

convert{T<:Real}(::Type{T}, z::ComplexDual) =
  ((epsilon(z)==0 && imag(z)==0 && imagepsilon==0) ? convert(T, real(z)) : throw(InexactError()))

promote_rule{T<:Real, S<:Real}(::Type{ComplexDual{T}}, ::Type{ComplexDual{S}}) =
    ComplexDual{promote_type(T, S)}

promote_rule{T<:Real, S<:Real}(::Type{Complex{T}}, ::Type{Dual{S}}) =
    ComplexDual{promote_type(T, S)}

promote_rule{T<:Real, S<:Real}(::Type{ComplexDual{T}}, ::Type{Dual{S}}) =
    ComplexDual{promote_type(T, S)}

promote_rule{T<:Real, S<:Real}(::Type{ComplexDual{T}}, ::Type{Complex{S}}) =
    ComplexDual{promote_type(T, S)}

promote_rule{T<:Real, S<:Real}(::Type{ComplexDual{T}}, ::Type{S}) =
    ComplexDual{promote_type(T, S)}


# these promotion rules shouldn't be used for scalar operations -- they're slow
#=
promote_rule{T<:Real}(::Type{Dual{T}}, ::Type{T}) = Dual{T}
promote_rule{T<:Real, S<:Real}(::Type{Dual{T}}, ::Type{S}) =
  Dual{promote_type(T, S)}
=#

complexdual(x, y, z, w) = ComplexDual(x, y, z, w)
complexdual(x) = ComplexDual(x)

@vectorize_1arg Real complexdual
#@vectorize_2arg Real complexdual
@vectorize_1arg ComplexDual epsilon
@vectorize_1arg ComplexDual imagepsilon

complexdual256(x::Float64, y::Float64, z::Float64, w::Float64) = ComplexDual{Float64}(x, y, z, w)
complexdual256(x::Real, y::Real) = complexdual128(float64(x), float64(y), float64(z), float64(w))
complexdual256(z) = complexdual256(real(z), epsilon(z), imag(z), imagepsilon(z))

iscomplexdual(x::ComplexDual) = true
iscomplexdual(x::Number) = false

real_valued{T<:Real}(z::ComplexDual{T}) = (epsilon(z) == 0 && imag(z) == 0 && imagepsilon(z)==0)
integer_valued(z::ComplexDual) = real_valued(z) && integer_valued(real(z))

isfinite(z::ComplexDual) = isfinite(real(z))
#reim(z::ComplexDual) = (real(z), epsilon(z)) 
flatten(z::ComplexDual) = (real(z), epsilon(z), imag(z), imagepsilon(z))

function complexdual_show(io::IO, z::ComplexDual, compact::Bool)
    x, y, z, w = flatten(z)
    print(io, "complexdual($x, $y, $z, $w)")
end
show(io::IO, z::ComplexDual) = complexdual_show(io, z, false)
showcompact(io::IO, z::ComplexDual) = complexdual_show(io, z, true)

function read{T<:Real}(s::IO, ::Type{ComplexDual{T}})
    x = read(s, T)
    y = read(s, T)
    z = read(s, T)
    w = read(s, T)

    ComplexDual{T}(x, y, z, w)
end

function write(s::IO, z::ComplexDual)
    write(s, real(z))
    write(s, epsilon(z))
    write(s, imag(z))
    write(s, imagepsilon(z))
end


## Generic functions of dual numbers ##

convert(::Type{ComplexDual}, z::ComplexDual) = z
convert(::Type{ComplexDual}, x::Real) = complexdual(x)

#=
==(z::Dual, w::Dual) = real(z) == real(w)
==(z::Dual, x::Real) = real(z) == x
==(x::Real, z::Dual) = real(z) == x
=#

#=
isequal(z::Dual, w::Dual) =
  isequal(real(z),real(w)) && isequal(epsilon(z), epsilon(w))
isequal(z::Dual, x::Real) = real_valued(z) && isequal(real(z), x)
isequal(x::Real, z::Dual) = real_valued(z) && isequal(real(z), x)

isless(z::Dual,w::Dual) = real(z) < real(w)
isless(z::Number,w::Dual) = z < real(w)
isless(z::Dual,w::Number) = real(z) < w


hash(z::Dual) =
  (x = hash(real(z)); real_valued(z) ? x : bitmix(x,hash(epsilon(z))))
=#

conj(z::ComplexDual) = ComplexDual(real(z), epsilon(z), -imag(z), -imagepsilon(z) )
abs(z::ComplexDual)  = abs(imagfull(z))
abs2(z::ComplexDual) = abs2(imagfull(z))

#= don't understand
# algebraic definitions
conjdual(z::ComplexDual) = ComplexDual(real(z),-epsilon(z))
absdual(z::ComplexDual) = abs(real(z))
abs2dual(z::ComplexDual) = abs2(real(z))
=#

+(a::ComplexDual, b::ComplexDual) = complexdual(real(a)+real(b), epsilon(a)+epsilon(b), imag(a)+imag(b), imagepsilon(a)+imagepsilon(b))
#+(a::Number, b::ComplexDual) = complexdual(a+real(b), epsilon(b), imag(b), imagepsilon(b))
#+(a::ComplexDual, b::Number) = complexdual(real(a)+b, epsilon(b), imag(b), imagepsilon(b))

-(z::ComplexDual) = complexdual(-real(z), -epsilon(z), -imag(z), -imagepsilon(z))
-(z::ComplexDual, w::ComplexDual) = z + (-w)
#-(z::Number, w::ComplexDual) = z + (-w)
#-(z::Dual, w::Number) = z + (-w)

# avoid ambiguous definition with Bool*Number
*(x::Bool, z::ComplexDual) = ifelse(x, z, ifelse(signbit(real(z))==0, zero(z), -zero(z)))
*(x::ComplexDual, z::Bool) = z*x

#16 terms total, but 4 of them are zero because eps^2 = 0
*(a::ComplexDual, b::ComplexDual) = complexdual(real(a)*real(b)     -imag(a)*imag(b), 
                                                epsilon(a)*real(b)  +real(a)*epsilon(b)  -imag(a)*imagepsilon(b)    -imagepsilon(a)*imag(b),
                                                real(a)*imag(b)     +imag(a)*real(b),
                                                real(a)*imagepsilon(b)  +imagepsilon(a)*real(b) +imag(a)*epsilon(b) +epsilon(a)*imag(b)
                                                )


#*(x::Real, z::ComplexDual) = complexdual(x*real(z), x*epsilon(z), x*imag(z), x*imagepsilon(z))
#*(z::ComplexDual, x::Real) = complexdual(real(z)*x, epsilon(z)*x, imag(z)*x, imagepsilon(z)*x)

#=

/(z::Real, w::Dual) = dual(z/real(w), -z*epsilon(w)/real(w)^2)
/(z::Dual, x::Real) = dual(real(z)/x, epsilon(z)/x)
=#


/(z::ComplexDual, w::ComplexDual) = z * conj(w) / abs2(z) 
#=
for f in [:^, :(NaNMath.pow)]
    @eval function ($f)(z::Dual, w::Dual)
        re = $f(real(z),real(w))

        du =
        epsilon(z)*real(w)*(($f)(real(z),real(w)-1))+epsilon(w)*($f)(real(z),real(w))*log(real(z))

        dual(re, du)
    end
end

# these two definitions are needed to fix ambiguity warnings
^(z::Dual, n::Integer) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))
^(z::Dual, n::Rational) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))

^(z::Dual, n::Real) = dual(real(z)^n, epsilon(z)*n*real(z)^(n-1))
NaNMath.pow(z::Dual, n::Real) = dual(NaNMath.pow(real(z),n), epsilon(z)*n*NaNMath.pow(real(z),n-1))

for (funsym, exp) in Calculus.derivative_rules
    @eval function $(funsym)(z::Dual)
        xp = epsilon(z)
        x = real(z)
        Dual($(funsym)(x),$exp)
    end
    # extend corresponding NaNMath methods
    if funsym in (:sin, :cos, :tan, :asin, :acos, :acosh, :atanh, :log, :log2, :log10,
          :lgamma, :log1p)
        funsym = Expr(:.,:NaNMath,Base.Meta.quot(funsym))
        @eval function $(funsym)(z::Dual)
            xp = epsilon(z)
            x = real(z)
            Dual($(funsym)(x),$exp)
        end
    end
end

=#
