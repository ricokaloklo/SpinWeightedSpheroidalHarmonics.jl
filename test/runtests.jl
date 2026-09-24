using SpinWeightedSpheroidalHarmonics
using LinearAlgebra
using Test

const SWSH = SpinWeightedSpheroidalHarmonics

function _theta_integral_abs2(harmonic; n=2001)
    theta = range(0.0, π; length=n)
    values = [abs2(harmonic(t, 0.0)) * sin(t) for t in theta]
    return step(theta) * (sum(values) - (first(values) + last(values)) / 2)
end

@testset "SpinWeightedSpheroidalHarmonics.jl" begin
    @testset "Spherical harmonics" begin
        for (s, l, m) in ((-2, 2, -2), (-1, 5, -2), (0, 0, 0), (2, 10, 7))
            direct = spin_weighted_spherical_harmonic(s, l, m; method="direct")
            jacobi = spin_weighted_spherical_harmonic(s, l, m; method="jacobi")
            for (theta, phi) in ((0.37, 0.91), (2.40, -0.20)), order in 0:2
                rtol = (2e-8, 1e-8, 1e-8)[order + 1]
                atol = (1e-10, 1e-10, 1e-9)[order + 1]
                @test jacobi(theta, phi; theta_derivative=order) ≈
                    direct(theta, phi; theta_derivative=order) rtol=rtol atol=atol
            end
        end
        @test spin_weighted_spherical_harmonic(-2, 10, 2).method == :direct
        @test spin_weighted_spherical_harmonic(-2, 40, 2).method == :jacobi
        @test spin_weighted_spherical_harmonic(-2, 40, 2; method="JaCoBi").method == :jacobi
        @test isfinite(spin_weighted_spherical_harmonic(-2, 200, 2)(1.1, 0.4))
        chebyshev = spin_weighted_spherical_harmonic(-2, 10, 2; method=:CHEBYSHEV)
        @test chebyshev(1.1, 0.4) ≈
            spin_weighted_spherical_harmonic(-2, 10, 2)(1.1, 0.4) atol=1e-10
    end

    @testset "Symmetry and normalization" begin
        negative = spin_weighted_spherical_harmonic(-2, 10, -3; method="jacobi")
        reflected = spin_weighted_spherical_harmonic(2, 10, 3; method="direct")
        for (theta, phi) in ((0.42, 0.30), (2.60, 2.20)), order in 0:1
            atol = (1e-12, 1e-10)[order + 1]
            @test negative(theta, phi; theta_derivative=order) ≈
                -conj(reflected(theta, phi; theta_derivative=order)) rtol=1e-9 atol=atol
        end
        direct = spin_weighted_spherical_harmonic(-2, 8, 2; method="direct")
        jacobi = spin_weighted_spherical_harmonic(-2, 8, 2; method="jacobi")
        @test jacobi(2π - 0.7, 0.4) ≈ direct(0.7, 0.4 + π) rtol=1e-9 atol=1e-12
        @test jacobi(-0.7, 0.4) ≈ direct(0.7, 0.4 + π) rtol=1e-9 atol=1e-12
        @test jacobi(2π - 0.7, 0.4; theta_derivative=1) ≈
            direct(0.7, 0.4 + π; theta_derivative=1) rtol=1e-9 atol=1e-10
        for (s, l, m) in ((-2, 6, 2), (0, 20, 0))
            direct = _theta_integral_abs2(spin_weighted_spherical_harmonic(s, l, m; method="direct"))
            jacobi = _theta_integral_abs2(spin_weighted_spherical_harmonic(s, l, m; method="jacobi"))
            @test jacobi ≈ direct rtol=1e-9 atol=1e-12
            @test direct ≈ 1 / (2π) rtol=5e-5 atol=5e-7
        end
    end

    @testset "Real spheroidicity" begin
        for (s, l, m, c) in ((-2, 2, 2, 0.68), (-2, 10, 2, 5.0), (-1, 8, -3, -4.0))
            fast = spin_weighted_spheroidal_harmonic(s, l, m, c)
            dense = spin_weighted_spheroidal_harmonic(
                s, l, m, c; N=fast.params.N, backend=:dense_reference)
            @test fast.lambda ≈ dense.lambda rtol=5e-14 atol=1e-12
            @test spin_weighted_spheroidal_eigenvalue(s, l, m, c) ≈
                dense.lambda rtol=5e-14 atol=1e-12
            @test fast.coeffs ≈ dense.coeffs rtol=2e-12 atol=2e-13
            for order in 0:1
                rtol = (2e-12, 3e-11)[order + 1]
                atol = (2e-13, 3e-12)[order + 1]
                @test fast(1.5, 0.8; theta_derivative=order) ≈
                    dense(1.5, 0.8; theta_derivative=order) rtol=rtol atol=atol
            end
        end
        for (s, l, m, c, reference) in (
                (1, 2, -2, -0.99999, -0.03510321969441066),
                (0, 2, 0, -9.9999, 54.50971251592067),
                (2, 2, 0, -45.0, 87.01704710570348))
            @test spin_weighted_spheroidal_eigenvalue(s, l, m, c) ≈
                reference rtol=3e-14 atol=5e-13
        end
        @test spin_weighted_spheroidal_eigenvalue(-2, 2, 0, -1e-7) ≈
            4.0000000000000047619047619047605 rtol=0 atol=1e-18
        @test spin_weighted_spheroidal_eigenvalue(1, 50, 0, -9.9999e-8) ≈
            2548.0 rtol=0 atol=1e-18
        @test spin_weighted_spheroidal_harmonic(-2, 6, 2, 0.0)(1.1, 0.7) ≈
            spin_weighted_spherical_harmonic(-2, 6, 2)(1.1, 0.7) rtol=1e-12 atol=1e-12
        @test_throws ArgumentError spin_weighted_spheroidal_harmonic(
            -2, 2, 2, 0.3; backend=:unknown)
    end

    @testset "Complex mode labels" begin
        c = 1.0 - 4.0im
        references = (-2.9585388471640925 + 24.185383504828405im,
                      -6.076746762274604 + 12.895628246290562im)
        for (l, reference) in zip(2:3, references)
            fast = spin_weighted_spheroidal_harmonic(-2, l, 2, c; N=33)
            dense = spin_weighted_spheroidal_harmonic(
                -2, l, 2, c; N=33, backend=:dense_reference)
            @test fast.lambda ≈ reference rtol=3e-13 atol=2e-12
            @test dense.lambda ≈ fast.lambda rtol=3e-13 atol=2e-12
            @test dense(1.1, 0.3) ≈ fast(1.1, 0.3) rtol=3e-12 atol=2e-12
            @test SWSH.spectral_coefficients(c, -2, l, 2, 33; backend=:dense_reference) ≈
                dense.coeffs rtol=3e-12 atol=2e-12
        end

        # Mathematica 14: SpheroidalEigenvalue[2, 2, I*c] - 4*c, at 40 digits.
        references = (37.768684073935812 - 40.343861362231600im,
                      50.371527920881364 - 27.546991700762978im,
                      53.065908639218583 - 12.568097753375649im)
        path = (-5.35 + 2.288im) .* exp.(im .* (-0.2, 0.0, 0.2))
        result = track_angular_mode(0, 2, 2, path)
        for (c, reference) in zip(path, references)
            state = result.states[argmin([abs(p.c - c) for p in result.states])]
            @test state.lambda ≈ reference rtol=5e-14 atol=2e-12
            @test spin_weighted_spheroidal_eigenvalue(0, 2, 2, c) ≈
                reference rtol=5e-14 atol=2e-12
        end
    end

    @testset "Angular equation" begin
        for (s, l, m, c) in ((-2, 2, 2, 1.0-4.0im), (-2, 3, 2, 1.0-4.0im),
                (0, 2, 2, -5.35+2.288im), (2, 3, -2, 0.25-0.15im))
            harmonic = spin_weighted_spheroidal_harmonic(s, l, m, c)
            A = harmonic.lambda - c^2 + 2m*c
            for theta in (0.3, 2.6)
                value = harmonic(theta, 0.0)
                first = harmonic(theta, 0.0; theta_derivative=1)
                second = harmonic(theta, 0.0; theta_derivative=2)
                potential = c^2*cos(theta)^2 - 2c*s*cos(theta) + s + A -
                    (m + s*cos(theta))^2 / sin(theta)^2
                terms = (second, cot(theta)*first, potential*value)
                @test abs(sum(terms)) / max(1, sum(abs, terms)) <= 1e-10
            end
        end
    end

    @testset "Continuation and cache" begin
        cache = AngularCache()
        positive = continue_angular_mode(-2, 3, 2, 0.4+0.15im; cache)
        negative = continue_angular_mode(-2, 3, 2, -0.4+0.15im; cache)
        dense = continue_angular_mode(-2, 3, 2, 0.4+0.15im; cache, backend=:dense_reference)
        @test positive.lambda != negative.lambda
        @test length(cache.values) == 3
        @test continue_angular_mode(-2, 3, 2, 0.4+0.15im; cache) === positive
        @test continue_angular_mode(
            -2, 3, 2, 0.4+0.15im; cache, backend=:dense_reference) === dense

        result = track_angular_mode(-2, 2, 2, [0.0, 0.15, 0.15+0.1im, 0.1im, 0.0])
        @test result.status == :closed
        @test result.lambda_closure_error <= 1e-12
        @test all(state.residual <= 5e-12 for state in result.states)
        @test all(abs(imag(dot(a.coefficients, b.coefficients))) <= 5e-13
            for (a, b) in zip(result.states, Iterators.drop(result.states, 1)))

        pair = continue_angular_mode(-2, 3, 2, 0.25-0.15im)
        harmonic = spin_weighted_spheroidal_harmonic(pair)
        observables = angular_observables(pair, 1.1, 0.3)
        @test harmonic(1.1, 0.3) ≈ observables.value
        @test harmonic(1.1, 0.3; theta_derivative=1) ≈ observables.first_derivative
        @test harmonic(1.1, 0.3; theta_derivative=2) ≈ observables.second_derivative
    end

    @testset "Lateral paths and mirror symmetry" begin
        lateral = continue_angular_lateral_pair(-2, 3, 1, 0.7, 0.4, 1e-4)
        right, left = last(lateral.right.states), last(lateral.left.states)
        @test (right.sheet_id, left.sheet_id) == (:right_lateral, :left_lateral)
        @test right.c == -conj(left.c)
        pair = continue_angular_mode(-2, 3, 1, 0.35-0.2im)
        parameters = mirror_angular_parameters(pair.s, pair.l, pair.m, pair.c)
        mirrored = continue_angular_mode(parameters.s, parameters.l, parameters.m, parameters.c)
        @test angular_mirror_residual(pair, mirrored) <= 2e-12
    end

    @testset "Precision and truncation" begin
        @test angular_precision_certificate(-2, 2, 2, 0.2-0.1im;
            truncation_orders=(24, 32), precision_bits=(128, 160)).accepted
    end
end
