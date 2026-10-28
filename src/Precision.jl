# Working precision. Every computation runs in the floating-point type T of its inputs, at least
# Float64; tolerances and truncation orders chosen for Float64 are carried to T by `_tol` and
# `_nterms`, which return the Float64 values unchanged. A harmonic built in BigFloat remembers the
# precision it was built with and evaluates under it (`_with_precision`).

# The working type for these inputs: their floating-point type, never narrower than Float64
# (Float32, integer and rational inputs compute in Float64)
_float_type(xs...) = float(promote_type(Float64, map(x -> real(typeof(x)), xs)...))
_float_type(::Type{T}) where {T} = float(promote_type(Float64, real(T)))

# An empirical Float64 tolerance sits between rounding (eps) and the physical scale (1): it
# keeps its place on the logarithmic scale between them, tol(T) = tol^(log eps(T)/log eps(Float64)).
_tol(::Type{T}, tol) where {T} = T(tol)^(log(eps(T)) / log(eps(Float64)))
_tol(::Type{Float64}, tol) = tol

# The number of terms of a series truncated after n64 terms in Float64, for the precision of T
# (a geometrically converging series needs a number of terms proportional to the number of digits)
_nterms(::Type{T}, n64) where {T} = cld(n64 * precision(T), 53)

# The precision a result in T is built and evaluated at (0 when T carries it in the type)
_precision_bits(::Type{BigFloat}) = precision(BigFloat)
_precision_bits(::Type) = 0

# Evaluate f() with BigFloat at precision p; other types carry their precision in the type
_with_precision(f, ::Type{BigFloat}, p) = setprecision(f, BigFloat, p)
_with_precision(f, ::Type, p) = f()

# The largest BigFloat precision among the inputs, under which a result is built
_input_precision(xs...) = maximum(_input_precision_bits, xs; init=0)
_input_precision_bits(x::BigFloat) = precision(x)
_input_precision_bits(x::Complex{BigFloat}) = max(precision(real(x)), precision(imag(x)))
_input_precision_bits(x::AbstractArray) = maximum(_input_precision_bits, x; init=0)
_input_precision_bits(x) = 0

#=
The public `precision` keyword: with an integer p, real and complex inputs are converted to
BigFloat at p bits (a Rational is rounded once); with `nothing` they are used as given.
=#
_with_precision_keyword(x, ::Nothing) = x
_with_precision_keyword(::Nothing, ::Integer) = nothing
_with_precision_keyword(x::Real, p::Integer) = setprecision(() -> BigFloat(x), BigFloat, p)
_with_precision_keyword(x::Complex, p::Integer) =
    setprecision(() -> Complex{BigFloat}(BigFloat(real(x)), BigFloat(imag(x))), BigFloat, p)
_with_precision_keyword(x::AbstractArray, p::Integer) = [_with_precision_keyword(v, p) for v in x]
# the precision keyword that rebuilds a part in the precision of the whole (0: the type carries it)
_precision_keyword(bits::Int) = bits == 0 ? nothing : bits
_check_precision_keyword(::Nothing) = nothing
_check_precision_keyword(p::Integer) =
    p >= 64 ? nothing : throw(ArgumentError("precision must be at least 64 bits, got $p."))
