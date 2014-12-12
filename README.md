### A ComplexDual Number

Is a four-dimensional algebra over Real numbers of the form
`a + b*ϵ + c*i + d*ϵi`

ComplexDual multiplication inherits from the Complex and Dual multiplication
|    | 1  | ϵ  | i  | ϵi |
|----|----|----|----|----|
| 1  | 1  | ϵ  | i  | ϵi |
| ϵ  | ϵ  | 0  | ϵi | 0  |
| i  | i  | ϵi | -1 | -ϵ |
| ϵi | ϵi | 0  | -ϵ | 0  |

A ComplexDual number can be interpreted as a Complex number over Dual numbers or as a Dual number over Complex numbers. These can be thought of as `Complex{Dual}`  (instead of something like `Complex{Float}`) or as a `Dual{Complex}`. These definitions are not possible with the implementation of `Complex` or `Dual`, because they are only defined over type `Real`.

There is some trouble in defining where in the type hierarchy `ComplexDual` goes. It should behave exactly like a `Complex` number, just as `Dual` behaves exactly like `Real`, but it is _not_ a `Complex` just like `Dual` is not a `Real`.

For example, a function of a `z::Complex{T<:Real}` might do something like `d = real(z)^2 + imag(z)^2`
