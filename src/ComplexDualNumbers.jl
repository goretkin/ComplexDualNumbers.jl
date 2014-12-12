module ComplexDualNumbers
  importall Base

  import NaNMath
  import Calculus

  importall DualNumbers
  include("complexdual.jl")

  export
    ComplexDual,
    ComplexDual256,
    ComplexDual128,
    ComplexDualPair,
    complexdual,
    complexdual256,
    complexdual128,
    iscomplexdual,
    complexdual_show,
    epsilon,
    imagepsilon,
    imagfull,
    epsfull,
    conjcomplexdual,
    abscomplexdual,
    abs2complexdual,
    dual_over_complex,
    complex_over_dual
end
