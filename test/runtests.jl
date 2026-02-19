using SpinWeightedSpheroidalHarmonics
using Test

@testset "SpinWeightedSpheroidalHarmonics.jl" begin
    @testset "Jacobi spherical harmonic evaluation" begin
        test_modes = [
            (-2, 6, 2),
            (-2, 6, -2),
            (-2, 8, 0),
            (0, 8, -3),
            (1, 6, -1),
            (2, 7, 0),
        ]
        test_angles = [
            (0.37, 0.91),
            (1.10, 2.30),
            (2.40, -0.20),
        ]

        for (s, l, m) in test_modes
            y_direct = spin_weighted_spherical_harmonic(s, l, m; method="direct")
            y_jacobi = spin_weighted_spherical_harmonic(s, l, m; method="jacobi")

            for (theta, phi) in test_angles
                @test y_jacobi(theta, phi) ≈ y_direct(theta, phi) rtol=1e-10 atol=1e-12
                @test y_jacobi(theta, phi; phi_derivative=2) ≈ y_direct(theta, phi; phi_derivative=2) rtol=1e-10 atol=1e-12
                @test y_jacobi(theta, phi; theta_derivative=1) ≈ y_direct(theta, phi; theta_derivative=1) rtol=1e-8 atol=1e-10
                @test y_jacobi(theta, phi; theta_derivative=2) ≈ y_direct(theta, phi; theta_derivative=2) rtol=1e-8 atol=1e-9
                @test y_jacobi(theta, phi; theta_derivative=1, phi_derivative=1) ≈ y_direct(theta, phi; theta_derivative=1, phi_derivative=1) rtol=1e-8 atol=1e-10
            end
        end
    end

    @testset "Auto method selection" begin
        @test spin_weighted_spherical_harmonic(-2, 10, 2; method="auto").method == :direct
        @test spin_weighted_spherical_harmonic(-2, 40, 2; method="auto").method == :jacobi
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
