var documenterSearchIndex = {"docs":
[{"location":"APIs/#APIs","page":"APIs","title":"APIs","text":"","category":"section"},{"location":"APIs/","page":"APIs","title":"APIs","text":"spin_weighted_spheroidal_harmonic","category":"page"},{"location":"APIs/#SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_harmonic","page":"APIs","title":"SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_harmonic","text":"spin_weighted_spheroidal_harmonic(s::Int, l::Int, m::Int, c, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0, N::Int=10)\n\nCompute the spin-weighted spheroidal harmonic with spin weight s, harmonic index l, azimuthal index m, spheroidicity c (c = aomega), and coordinates theta and phi.\n\nThe optional arguments theta_derivative and phi_derivative specify the order of partial derivatives to take with respect to theta and phi, respectively.\n\nThe optional argument N specifies the number of terms to use in the spectral decomposition. The default value is N=10.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"spin_weighted_spheroidal_eigenvalue","category":"page"},{"location":"APIs/#SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_eigenvalue","page":"APIs","title":"SpinWeightedSpheroidalHarmonics.spin_weighted_spheroidal_eigenvalue","text":"spin_weighted_spheroidal_eigenvalue(s::Int, l::Int, m::Int, c; N::Int=10)\n\nCompute the eigenvalue of the spin-weighted spheroidal harmonic with spin weight s, harmonic index l, azimuthal index m, and spheroidicity c (c = aomega).\n\nThe optional argument N specifies the number of terms to use in the spectral decomposition. The default value is N=10.\n\nThis function is simply a wrapper to Teukolsky_lambda_const for backward compatbility.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"spin_weighted_spherical_harmonic","category":"page"},{"location":"APIs/#SpinWeightedSpheroidalHarmonics.spin_weighted_spherical_harmonic","page":"APIs","title":"SpinWeightedSpheroidalHarmonics.spin_weighted_spherical_harmonic","text":"spin_weighted_spherical_harmonic(s::Int, l::Int, m::Int, theta, phi; theta_derivative::Int=0, phi_derivative::Int=0)\n\nCompute the spin-weighted spherical harmonic with spin weight s, harmonic index l, azimuthal index m, and coordinates theta and phi.\n\nThe optional arguments theta_derivative and phi_derivative specify the order of partial derivatives to take with respect to theta and phi, respectively.\n\n\n\n\n\n","category":"function"},{"location":"APIs/","page":"APIs","title":"APIs","text":"spin_weighted_spherical_eigenvalue","category":"page"},{"location":"APIs/#SpinWeightedSpheroidalHarmonics.spin_weighted_spherical_eigenvalue","page":"APIs","title":"SpinWeightedSpheroidalHarmonics.spin_weighted_spherical_eigenvalue","text":"spin_weighted_spherical_eigenvalue(s::Int, l::Int, m::Int=0)\n\nCompute the eigenvalue of the spin-weighted spherical harmonic with spin weight s, harmonic index l, and azimuthal index m (but the eigenvalue is independent of m).\n\n\n\n\n\n","category":"function"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"#Home","page":"Home","title":"Home","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SpinWeightedSpheroidalHarmonics.jl computes spin-weighted spheroidal harmonics and eigenvalues using a spectral decomposition method.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The two main features are implemented as","category":"page"},{"location":"","page":"Home","title":"Home","text":"spin_weighted_spheroidal_harmonic for computing the harmonic, and\nspin_weighted_spheroidal_eigenvalue for computing the eigenvalue","category":"page"},{"location":"","page":"Home","title":"Home","text":"and both supporting complex spheroidicity (and hence frequency omega). See Quick-start below for some simple examples.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Additionally, we provide two similar functions","category":"page"},{"location":"","page":"Home","title":"Home","text":"spin_weighted_spherical_harmonic, and\nspin_weighted_spherical_eigenvalue","category":"page"},{"location":"","page":"Home","title":"Home","text":"that return the exact harmonic and eigenvalue respectively.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Exact partial derivatives (with respect to either theta and/or phi) can be evaluated by specifying the derivative order with theta_derivative and phi_derivative respectively when calling the functions for a harmoic.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package using the Julia package manager, simply type the following in the Julia REPL:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"SpinWeightedSpheroidalHarmonics\")","category":"page"},{"location":"#Quick-start","page":"Home","title":"Quick-start","text":"","category":"section"},{"location":"#Computing-the-spin-weighted-spheroidal-eigenvalue","page":"Home","title":"Computing the spin-weighted spheroidal eigenvalue","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For example, to compute the spin-weighted spheroidal eigenvalue lambda for the mode s = -2 ell = 2 m = 2 a = 07 omega = 05, simply do","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SpinWeightedSpheroidalHarmonics\ns=-2; l=2; m=2; a=0.7; omega=0.5;\nspin_weighted_spheroidal_eigenvalue(s, l, m, a*omega)","category":"page"},{"location":"#Computing-the-spin-weighted-spheroidal-harmonic","page":"Home","title":"Computing the spin-weighted spheroidal harmonic","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For example, to compute the spin-weighted spheroidal harmonic for the mode s = -2 ell = 2 m = 2 a = 07 omega = 05 at theta = pi6 phi = pi3, simply do","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SpinWeightedSpheroidalHarmonics\ns=-2; l=2; m=2; a=0.7; omega=0.5;\ntheta=π/6; phi=π/3;\nspin_weighted_spheroidal_harmonic(s, l, m, a*omega, theta, phi)","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The package is licensed under the MIT License.","category":"page"}]
}
