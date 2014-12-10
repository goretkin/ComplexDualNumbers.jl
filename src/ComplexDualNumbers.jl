module ComplexDualNumbers
  importall Base

  import NaNMath
  import Calculus

  include("complexdual.jl")

  export
    ComplexDual,
    ComplexDual128,
    ComplexDual64,
    ComplexDualPair,
    complexdual,
    complexdual128,
    complexdual64,
    iscomplexdual,
    complexdual_show,
    epsilon,
    conjcomplexdual,
    abscomplexdual,
    abs2complexdual
end
