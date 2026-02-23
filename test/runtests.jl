using SpinWeightedSpheroidalHarmonics
using Test

function _theta_integral_abs2(swsh; n::Int=2001, phi=0.0)
    thetas = range(0.0, π; length=n)
    h = step(thetas)
    vals = [abs2(swsh(θ, phi)) * sin(θ) for θ in thetas]
    # Trapezoidal rule on [0, π]
    return h * (sum(vals) - 0.5 * (vals[1] + vals[end]))
end

@testset "SpinWeightedSpheroidalHarmonics.jl" begin
    @testset "Jacobi spherical harmonic evaluation" begin
        test_spins = (-2, -1, 0, 1, 2)
        test_angles = [
            (0.37, 0.91),
            (1.10, 2.30),
            (2.40, -0.20),
        ]

        for s in test_spins
            for l in abs(s):20
                for m in -l:l
                    y_direct = spin_weighted_spherical_harmonic(s, l, m; method="direct")
                    y_jacobi = spin_weighted_spherical_harmonic(s, l, m; method="jacobi")

                    for (theta, phi) in test_angles
                        @test y_jacobi(theta, phi) ≈ y_direct(theta, phi) rtol=2e-8 atol=1e-10
                        @test y_jacobi(theta, phi; phi_derivative=2) ≈ y_direct(theta, phi; phi_derivative=2) rtol=2e-8 atol=1e-10
                        @test y_jacobi(theta, phi; theta_derivative=1) ≈ y_direct(theta, phi; theta_derivative=1) rtol=1e-8 atol=1e-10
                        @test y_jacobi(theta, phi; theta_derivative=2) ≈ y_direct(theta, phi; theta_derivative=2) rtol=1e-8 atol=1e-9
                        @test y_jacobi(theta, phi; theta_derivative=1, phi_derivative=1) ≈ y_direct(theta, phi; theta_derivative=1, phi_derivative=1) rtol=1e-8 atol=1e-10
                    end
                end
            end
        end
    end

    @testset "Auto method selection" begin
        @test spin_weighted_spherical_harmonic(-2, 10, 2; method="auto").method == :direct
        @test spin_weighted_spherical_harmonic(-2, 40, 2; method="auto").method == :jacobi
    end

    @testset "Jacobi symmetry tricks: m reflection" begin
        # Jacobi uses a symmetry relation for m < 0. Validate it against
        # direct evaluation, which is treated as the golden standard.
        symmetry_pairs = [
            (-2, 10, 3),
            (-1, 9, 2),
            (1, 8, 2),
            (2, 11, 1),
        ]
        symmetry_angles = [
            (0.42, 0.30),
            (1.70, -0.80),
            (2.60, 2.20),
        ]

        for (s, l, mpos) in symmetry_pairs
            mneg = -mpos
            y_jacobi_neg = spin_weighted_spherical_harmonic(s, l, mneg; method="jacobi")
            y_direct_neg = spin_weighted_spherical_harmonic(s, l, mneg; method="direct")
            y_direct_ref = spin_weighted_spherical_harmonic(-s, l, mpos; method="direct")
            phase = (-1)^(s - mneg)

            for (theta, phi) in symmetry_angles
                # Golden-standard agreement for the negative-m mode.
                @test y_jacobi_neg(theta, phi) ≈ y_direct_neg(theta, phi) rtol=1e-9 atol=1e-12
                @test y_jacobi_neg(theta, phi; theta_derivative=1) ≈ y_direct_neg(theta, phi; theta_derivative=1) rtol=1e-9 atol=1e-10
                @test y_jacobi_neg(theta, phi; theta_derivative=2) ≈ y_direct_neg(theta, phi; theta_derivative=2) rtol=1e-9 atol=1e-9

                # Explicitly check the symmetry trick against direct.
                @test y_jacobi_neg(theta, phi) ≈ phase * conj(y_direct_ref(theta, phi)) rtol=1e-9 atol=1e-12
                @test y_jacobi_neg(theta, phi; theta_derivative=1) ≈ phase * conj(y_direct_ref(theta, phi; theta_derivative=1)) rtol=1e-9 atol=1e-10
                @test y_jacobi_neg(theta, phi; theta_derivative=1, phi_derivative=1) ≈ phase * conj(y_direct_ref(theta, phi; theta_derivative=1, phi_derivative=1)) rtol=1e-9 atol=1e-10
            end
        end
    end

    @testset "Jacobi symmetry tricks: theta folding" begin
        # Jacobi maps theta to [0, pi] using (theta, phi) -> (2pi-theta, phi+pi)
        # when needed. Validate this mapping with direct evaluation.
        fold_modes = [
            (-2, 8, 2),
            (-2, 8, -2),
            (0, 9, 0),
            (1, 7, -1),
        ]
        principal_angles = [
            (0.70, 0.40),
            (1.40, -1.10),
            (2.20, 2.20),
        ]

        for (s, l, m) in fold_modes
            y_direct = spin_weighted_spherical_harmonic(s, l, m; method="direct")
            y_jacobi = spin_weighted_spherical_harmonic(s, l, m; method="jacobi")

            for (theta, phi) in principal_angles
                theta_folded = 2π - theta
                theta_negative = -theta

                @test y_jacobi(theta_folded, phi) ≈ y_direct(theta, phi + π) rtol=1e-9 atol=1e-12
                @test y_jacobi(theta_negative, phi) ≈ y_direct(theta, phi + π) rtol=1e-9 atol=1e-12
                @test y_jacobi(theta_folded, phi; phi_derivative=2) ≈ y_direct(theta, phi + π; phi_derivative=2) rtol=1e-9 atol=1e-10
                @test y_jacobi(theta_negative, phi; phi_derivative=2) ≈ y_direct(theta, phi + π; phi_derivative=2) rtol=1e-9 atol=1e-10
                @test y_jacobi(theta_folded, phi; theta_derivative=1) ≈ y_direct(theta, phi + π; theta_derivative=1) rtol=1e-9 atol=1e-10
                @test y_jacobi(theta_negative, phi; theta_derivative=1) ≈ y_direct(theta, phi + π; theta_derivative=1) rtol=1e-9 atol=1e-10
            end
        end
    end

    @testset "Case-insensitive method options" begin
        @test spin_weighted_spherical_harmonic(-2, 40, 2; method="JaCoBi").method == :jacobi
        @test spin_weighted_spherical_harmonic(-2, 40, 2; method=:CHEBYSHEV).method == :chebyshev
        @test spin_weighted_spherical_harmonic(-2, 10, 2; method="DiReCt").method == :direct
        @test spin_weighted_spherical_harmonic(-2, 10, 2; method=:AUTO).method == :direct
    end

    @testset "High-l Jacobi is finite" begin
        y = spin_weighted_spherical_harmonic(-2, 200, 2; method="jacobi")
        val = y(1.1, 0.4)
        @test isfinite(real(val))
        @test isfinite(imag(val))
    end

    @testset "Jacobi normalization convention" begin
        modes = [
            (-2, 6, 2),
            (-2, 12, 2),
            (-1, 15, 1),
            (0, 20, 0),
            (2, 18, -1),
        ]

        for (s, l, m) in modes
            y_direct = spin_weighted_spherical_harmonic(s, l, m; method="direct")
            y_jacobi = spin_weighted_spherical_harmonic(s, l, m; method="jacobi")

            norm_direct_phi0 = _theta_integral_abs2(y_direct; phi=0.0)
            norm_jacobi_phi0 = _theta_integral_abs2(y_jacobi; phi=0.0)
            norm_direct_phi1 = _theta_integral_abs2(y_direct; phi=1.3)
            norm_jacobi_phi1 = _theta_integral_abs2(y_jacobi; phi=1.3)

            # Jacobi should match direct normalization for the same mode.
            @test norm_jacobi_phi0 ≈ norm_direct_phi0 rtol=1e-9 atol=1e-12
            @test norm_jacobi_phi1 ≈ norm_direct_phi1 rtol=1e-9 atol=1e-12

            # Direct evaluation method follows the package convention:
            # ∫_0^π |sY_lm(θ, ϕ)|^2 sinθ dθ = 1/(2π), independent of ϕ.
            @test norm_direct_phi0 ≈ 1 / (2π) rtol=5e-5 atol=5e-7
            @test norm_direct_phi1 ≈ 1 / (2π) rtol=5e-5 atol=5e-7
        end
    end

    @testset "Spheroidal lambda consistency" begin
        s = -2
        l = 8
        m = 2
        c = 0.3
        swsh = spin_weighted_spheroidal_harmonic(s, l, m, c; method="JaCoBi")
        @test swsh.lambda ≈ spin_weighted_spheroidal_eigenvalue(s, l, m, c)
    end

    @testset "Spherical limit (c = 0)" begin
        s = -2
        l = 6
        m = 2
        spheroidal = spin_weighted_spheroidal_harmonic(s, l, m, 0.0; method="direct")
        spherical = spin_weighted_spherical_harmonic(s, l, m; method="direct")
        @test spheroidal(1.1, 0.7) ≈ spherical(1.1, 0.7) rtol=1e-12 atol=1e-12
    end
end
