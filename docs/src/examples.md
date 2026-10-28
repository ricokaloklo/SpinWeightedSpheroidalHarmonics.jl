# Examples

Two Jupyter notebooks in the [`examples`](https://github.com/ricokaloklo/SpinWeightedSpheroidalHarmonics.jl/tree/main/examples)
folder of the repository go with this page: `plot-SWSH.ipynb` makes the plots below, and
`arbitrary-precision.ipynb` walks through computing in `BigFloat`, complex frequencies, the accuracy of
each method and the cost.

# Example 1: Solving and visualizing a spin-weighted spheroidal harmonic
In this example, we solve for the spin-weighted spheroidal harmonic with $s=-2, \ell =2, m = 2, c \equiv a\omega = 0.35$.
```julia
using SpinWeightedSpheroidalHarmonics
using Plots, LaTeXStrings

s=-2; l=2; m=2; a=0.7; omega=0.5;
swsh = spin_weighted_spheroidal_harmonic(s, l, m, a*omega)
```
That's it! Now let us visualize the harmonic itself and the first and second partial derivative with respect to $\theta$:
```julia
thetas = collect(0:0.01:1)
plot(thetas, [real(swsh(theta*pi, 0)) for theta in thetas], linewidth=2, label=L"{}_{-2}S_{2,2}(\theta, 0; c = 0.35)", dpi=300)
plot!(thetas, [real(swsh(theta*pi, 0; theta_derivative=1)) for theta in thetas], linewidth=2, label=L"\partial_{\theta} {}_{-2}S_{2,2}(\theta, 0; c = 0.35)")
plot!(thetas, [real(swsh(theta*pi, 0; theta_derivative=2)) for theta in thetas], linewidth=2, label=L"\partial_{\theta}^2 {}_{-2}S_{2,2}(\theta, 0; c = 0.35)")
plot!(
    legendfontsize=14,
    xguidefontsize=14,
    yguidefontsize=14,
    xtickfontsize=14,
    ytickfontsize=14,
    foreground_color_legend=nothing,
    background_color_legend=nothing,
    legend=:bottomright,
    formatter=:latex,
    xlabel=L"\theta/\pi",
    left_margin = 2Plots.mm,
    right_margin = 3Plots.mm,
)
```
![SWSH.png](SWSH.png)

# Example 2: Solving and visualizing a spin-weighted spherical harmonic
In this example, we solve for the spin-weighted _spherical_ harmonic with $s=-2, \ell =200, m = 10$. The highly oscillatory nature could be challenging for some codes. Let us see if this is the case here. 
```julia
using SpinWeightedSpheroidalHarmonics
using Plots, LaTeXStrings

s=-2; l=200; m=10;
Y = spin_weighted_spherical_harmonic(s, l, m)
```
Now we visualize the harmonic:
```julia
thetas = collect(0:0.001:1)
plot(thetas, [real(Y(float(theta*π), 0)) for theta in thetas], linewidth=2, label=nothing, dpi=300)
ylabel!(L"{}_{-2}Y_{200,10}(\theta, 0)")
plot!(
    legendfontsize=14,
    xguidefontsize=14,
    yguidefontsize=14,
    xtickfontsize=14,
    ytickfontsize=14,
    foreground_color_legend=nothing,
    background_color_legend=nothing,
    legend=:bottomright,
    formatter=:latex,
    xlabel=L"\theta/\pi",
    left_margin = 2Plots.mm,
    right_margin = 3Plots.mm,
)
```
![Y.png](Y.png)
Looking good!

# Example 3: A complex frequency
For a complex frequency nothing changes in the call. Here ``a = 0.7`` and ``\omega`` is close to the
fundamental quasinormal-mode frequency of the ``\ell = m = 2`` mode:
```julia
using SpinWeightedSpheroidalHarmonics

a = 0.7; omega = 0.5326 - 0.0808im
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, a * omega)
S.lambda
S(1.1, 0.3), S(1.1, 0.3; theta_derivative=2)
```
The harmonic is normalized with its modulus, and its phase is fixed by making its overlap with
``{}_{-2}Y_{22}`` real and positive. For a quasinormal-mode residue, which needs the analytic
continuation of the harmonic in `c`, normalize ``\int S^2 \, d\Omega = 1`` instead (see [APIs](@ref)):
```julia
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, a * omega; normalization=:analytic)
```

# Example 4: More digits
The same harmonic to 256 bits. A rational input is rounded once, in the target precision, and the
harmonic evaluates at 256 bits wherever it is called:
```julia
S = spin_weighted_spheroidal_harmonic(-2, 2, 2, 7//10 * (5326//10000 - 808im//10000); precision=256)
S(big(11) / 10, 0.3)
```
See [Arbitrary precision](@ref) for the details.
