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

    @testset "Leaver eigenvalue cross-check" begin
        # The spectral decomposition is treated as the golden standard here.
        cases = [
            (-2, 2, 2, 0.35),
            (-2, 4, 2, 0.28),
            (-1, 3, 1, 0.135),
            (0, 5, 2, 0.24),
            (2, 4, 1, 0.1375),
            (-2, 2, -2, -0.6),
            (-2, 8, 2, 1.4),
            (0, 2, 0, 2.5),
        ]

        for (s, l, m, c) in cases
            @test spin_weighted_spheroidal_eigenvalue(s, l, m, c; method="leaver") ≈
                  spin_weighted_spheroidal_eigenvalue(s, l, m, c) rtol=1e-9 atol=1e-10
        end

        # Complex spheroidicity, as needed for a complex frequency.
        c = 0.7 * (0.2 + 0.15im)
        @test spin_weighted_spheroidal_eigenvalue(-2, 3, 2, c; method="leaver") ≈
              spin_weighted_spheroidal_eigenvalue(-2, 3, 2, c) rtol=1e-8 atol=1e-9

        # A wider working precision survives a spheroidicity at which Float64 gives up,
        # with the tolerance following the precision of c on its own.
        setprecision(BigFloat, 128) do
            λ = spin_weighted_spheroidal_eigenvalue(-2, 2, 2, big(8.0); method="leaver")
            @test Float64(λ) ≈ spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 8.0; N=120) rtol=1e-12 atol=1e-12
            # ... whereas Float64 refuses rather than returning an eigenvalue it cannot resolve.
            @test_throws ErrorException spin_weighted_spheroidal_eigenvalue(-2, 2, 2, 8.0; method="leaver")
        end
    end

    @testset "Leaver harmonic cross-check" begin
        cases = [
            (-2, 2, 2, 0.35),
            (-2, 4, 2, 0.28),
            (-1, 3, 1, 0.135),
            (0, 5, 2, 0.24),
            (2, 4, 1, 0.1375),
            (-2, 2, -2, -0.6),
        ]
        sample_points = [
            (0.45, 0.30),
            (1.10, -0.70),
            (2.30, 1.80),
        ]

        for (s, l, m, c) in cases
            swsh_leaver = spin_weighted_spheroidal_harmonic(s, l, m, c; method="leaver")
            swsh_ref = spin_weighted_spheroidal_harmonic(s, l, m, c; method="jacobi")

            @test swsh_leaver.method == :leaver
            @test swsh_leaver.lambda ≈ swsh_ref.lambda rtol=1e-9 atol=1e-10

            for (theta, phi) in sample_points
                # Including the overall phase: leaver follows the same convention.
                @test swsh_leaver(theta, phi) ≈ swsh_ref(theta, phi) rtol=1e-7 atol=1e-9
                @test swsh_leaver(theta, phi; phi_derivative=2) ≈ swsh_ref(theta, phi; phi_derivative=2) rtol=1e-7 atol=1e-9
                @test swsh_leaver(theta, phi; theta_derivative=1) ≈ swsh_ref(theta, phi; theta_derivative=1) rtol=1e-7 atol=1e-9
                @test swsh_leaver(theta, phi; theta_derivative=2) ≈ swsh_ref(theta, phi; theta_derivative=2) rtol=1e-7 atol=1e-8
                @test swsh_leaver(theta, phi; theta_derivative=1, phi_derivative=1) ≈ swsh_ref(theta, phi; theta_derivative=1, phi_derivative=1) rtol=1e-7 atol=1e-9

                # theta outside [0, π] is folded back the same way.
                @test swsh_leaver(2π - theta, phi) ≈ swsh_ref(theta, phi + π) rtol=1e-7 atol=1e-9
                @test swsh_leaver(-theta, phi) ≈ swsh_ref(theta, phi + π) rtol=1e-7 atol=1e-9
            end

            # Package convention: ∫_0^π |S(θ, ϕ)|^2 sinθ dθ = 1/(2π).
            @test _theta_integral_abs2(swsh_leaver; n=2501) ≈ 1 / (2π) rtol=5e-5 atol=5e-7
        end

        # Spherical limit, where the branch is no longer visible in the recurrence.
        spheroidal = spin_weighted_spheroidal_harmonic(-2, 6, 2, 0.0; method="leaver")
        spherical = spin_weighted_spherical_harmonic(-2, 6, 2; method="direct")
        @test spheroidal(1.1, 0.7) ≈ spherical(1.1, 0.7) rtol=1e-10 atol=1e-12

        #=
        Powers of u are a poorly conditioned basis at large l, so the harmonic
        (unlike the eigenvalue, which stays accurate) runs out of precision.
        =#
        @test spin_weighted_spheroidal_eigenvalue(-2, 20, 2, 0.6; method="leaver") ≈
              spin_weighted_spheroidal_eigenvalue(-2, 20, 2, 0.6) rtol=1e-8 atol=1e-8
        @test_throws ErrorException spin_weighted_spheroidal_harmonic(-2, 20, 2, 0.6; method="leaver")

        # A wider working precision buys back the harmonic that Float64 refuses.
        setprecision(BigFloat, 192) do
            wide = spin_weighted_spheroidal_harmonic(-2, 20, 2, big(0.6); method="leaver")
            reference = spin_weighted_spheroidal_harmonic(-2, 20, 2, 0.6; method="jacobi")
            @test ComplexF64(wide(big(1.1), big(0.4))) ≈ reference(1.1, 0.4) rtol=1e-10 atol=1e-12
        end

        # Supplying the eigenvalue skips the continued-fraction solve.
        s, l, m, c = -2, 3, 2, 0.4
        λ = spin_weighted_spheroidal_eigenvalue(s, l, m, c)
        given = spin_weighted_spheroidal_harmonic(s, l, m, c; method="leaver", lambda=λ)
        @test given.lambda == λ
        @test given(1.1, 0.4) ≈ spin_weighted_spheroidal_harmonic(s, l, m, c; method="leaver")(1.1, 0.4) rtol=1e-8 atol=1e-10
    end

    @testset "Leaver accepts a non-integer harmonic index" begin
        # l is only restricted to an integer by the spectral decomposition.
        @test spin_weighted_spheroidal_harmonic(-2, 3.0 + 0.1im, 2, 0.28; branch_n=1).method == :leaver
        @test spin_weighted_spheroidal_harmonic(-2, 3.5, 2, 0.28; branch_n=1).method == :leaver
        @test spin_weighted_spheroidal_harmonic(-2, 3, 2, 0.28; method="LeAvEr").method == :leaver

        value = spin_weighted_spheroidal_harmonic(-2, 3.0 + 0.1im, 2, 0.28; branch_n=1)(1.1, 0.4)
        @test isfinite(real(value))
        @test isfinite(imag(value))

        λ = spin_weighted_spheroidal_eigenvalue(-2, 3.0 + 0.1im, 2, 0.28; branch_n=1)
        @test isfinite(real(λ))
        @test isfinite(imag(λ))

        #=
        The continued fraction has a discrete spectrum, so a non-integer l does
        not give a new eigenvalue: it only moves the initial guess, which here
        lands on the neighboring l = 4 mode. Every inversion of the continued
        fraction shares the same roots, so branch_n does not override that.
        =#
        @test spin_weighted_spheroidal_eigenvalue(-2, 3.7, 2, 0.28; branch_n=1) ≈
              spin_weighted_spheroidal_eigenvalue(-2, 4, 2, 0.28) rtol=1e-9 atol=1e-10
        # Starting the solver elsewhere is what picks out another branch.
        @test spin_weighted_spheroidal_eigenvalue(-2, 3.7, 2, 0.28; branch_n=1, lambda0=8.5) ≈
              spin_weighted_spheroidal_eigenvalue(-2, 3, 2, 0.28) rtol=1e-9 atol=1e-10

        @test_throws ErrorException spin_weighted_spheroidal_harmonic(-2, 3.5, 2, 0.28; method="jacobi")
        @test_throws ErrorException spin_weighted_spheroidal_eigenvalue(-2, 3.5, 2, 0.28; method="chebyshev")
        @test_throws ErrorException spin_weighted_spherical_harmonic(-2, 3, 2; method="leaver")
    end
end
