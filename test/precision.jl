# Working precision: the precision helpers leave Float64 unchanged; BigFloat inputs (or the
# precision keyword) give BigFloat results whose accuracy follows the precision; harmonics remember
# the precision they were built with.
using SpinWeightedSpheroidalHarmonics, Test
const SWSH = SpinWeightedSpheroidalHarmonics

@testset "Float64 is unchanged by the precision helpers" begin
    for tol in (1e-3, 1e-6, 5e-12, 5e-13, 5e-14, 5e-15, 1e-14, 1e-16, 1e-30)
        @test SWSH._tol(Float64, tol) === tol
    end
    for n in (64, 4096, 6400)
        @test SWSH._nterms(Float64, n) === n
    end
    T = Float64
    @test T(π) / 2 === π / 2 && 2T(π) === 2π && 4 * T(π) === 4π
    @test SWSH._eigenvector_gap_tol(T) === eps(Float64) / sqrt(2 * 5e-13)
    @test SWSH._float_type(1) === Float64 && SWSH._float_type(1//2) === Float64
    @test SWSH._float_type(0.5f0) === Float64 && SWSH._float_type(0.5f0 + 1f0im) === Float64
    @test SWSH._float_type(big(1)) === BigFloat && SWSH._float_type(1.0 + big(2)im) === BigFloat
    @test SWSH._with_precision(() -> precision(BigFloat), Float64, 512) == precision(BigFloat)
    @test SWSH._with_precision(() -> precision(BigFloat), BigFloat, 512) == 512
end

@testset "Result types" begin
    # The eigenvalue follows the type of c, also at c = 0; the spherical eigenvalue is an exact integer
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.0) === 4.0
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0) === 4.0
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.0; method="leaver") === 4.0
    @test spin_weighted_spherical_eigenvalue(-2, 2) === 4
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.5f0) ===
        spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.5)
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 1//2) ===
        spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.5)
end

# Mathematica 14, BHPToolkit SpinWeightedSpheroidalHarmonics 1.0.1, 80-digit working precision:
# SpinWeightedSpheroidalEigenvalue[-2, 2, 2, 1/2] and SpinWeightedSpheroidalHarmonicS[-2, 2, 2, 1/2][11/10, 0]
const LAMBDA_REF = "0.725702761257765110195306308936146146138102156529773250230740"
const S_REF = "0.305813904079861803391934243129743924877990486622840116129065"

@testset "BigFloat results reach the precision" begin
    setprecision(BigFloat, 256) do
        λref, Sref = parse(BigFloat, LAMBDA_REF), parse(BigFloat, S_REF)
        θ = big(11) / 10
        for method in ("auto", "leaver")
            λ = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, big(1) / 2; method)
            @test λ isa BigFloat && abs(λ - λref) < 1e-58
        end
        for method in ("auto", "jacobi", "chebyshev", "leaver")
            h = spin_weighted_spheroidal_harmonic(-2, 2, 2, big(1) / 2; method)
            v = h(θ, big(0))
            @test v isa Complex{BigFloat} && abs(v - Sref) < 1e-58
        end
        # complex c: two precisions agree far beyond Float64
        c = big(1) / 2 - big(1)im / 2
        λ256 = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, c)
        λ384 = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, c; precision=384)
        @test λ256 isa Complex{BigFloat} && abs(λ256 - λ384) < 1e-70
    end
end

@testset "The precision keyword and remembered precision" begin
    h = spin_weighted_spheroidal_harmonic(-2, 2, 2, 1//2; precision=256)
    v = setprecision(() -> h(1.1, 0.0), BigFloat, 64)
    @test v isa Complex{BigFloat} && precision(real(v)) == 256
    @test abs(v - setprecision(() -> parse(BigFloat, S_REF), BigFloat, 256)) < 1e-15
    @test spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 1//2; precision=256) ==
        setprecision(() -> spin_weighted_spheroidal_eigenvalue(-2, 2, 2, big(1) / 2), BigFloat, 256)
    Y = spin_weighted_spherical_harmonic(-2, 3, 1; precision=200)
    y = setprecision(() -> Y(0.3, 0.2; theta_derivative=1), BigFloat, 64)
    @test y isa Complex{BigFloat} && precision(real(y)) == 200
    pair = continue_angular_mode(-2, 3, 2, 0.25 - 0.15im; precision=256, cache=nothing)
    @test pair isa AngularEigenpair{BigFloat} && pair.precision_bits == 256
    @test abs(pair.lambda - continue_angular_mode(-2, 3, 2, 0.25 - 0.15im; cache=nothing).lambda) < 1e-12
    @test_throws ArgumentError spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 0.5; precision=32)
end

@testset "Modes Float64 cannot resolve" begin
    # The l = 2 mode and its l = 3 neighbour are closer than Float64 resolves (relative gap ~1e-13 at
    # c = -20): Float64 throws, a wider precision selects the mode and resolves the harmonic
    @test_throws ErrorException spin_weighted_spheroidal_harmonic(-2, 2, 2, -20.0)
    for c in (-20, -20 - 1im//5)
        h = spin_weighted_spheroidal_harmonic(-2, 2, 2, c; precision=256)
        λ3 = spin_weighted_spheroidal_eigenvalue(-2, 3, 2, c; precision=256)
        @test abs(h.lambda - λ3) > 1e-14 * abs(λ3)
        @test abs(h(1.1, 0.3)) > 0
    end
    # the complex harmonic follows the precision (phase transported along the path in T)
    c = 1//2 - 1im//2
    v(bits) = (h = spin_weighted_spheroidal_harmonic(-2, 2, 2, c; precision=bits);
        setprecision(() -> h(big(11) / 10, big(3) / 10; theta_derivative=1), BigFloat, 400))
    a, b = v(256), v(384)
    @test setprecision(() -> abs(a - b) / abs(b), BigFloat, 400) < 1e-70
end
