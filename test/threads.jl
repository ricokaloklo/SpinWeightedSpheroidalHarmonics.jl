# Concurrent evaluation: building and evaluating harmonics, eigenvalues and continued modes from
# several threads at once gives the same bits as doing it serially, in Float64 and at 256 bits.
using SpinWeightedSpheroidalHarmonics, Test

const THETAS = (0.0, 0.3, 1.1, π / 2, 2.9, π, 4.0)

# The outputs of one job, each built afresh so that lazily built parts start empty; an exception is
# recorded by its type
function _outputs(job)
    try
        return job()
    catch err
        return Any[typeof(err)]
    end
end

function harmonic_job(s, l, m, c, method; precision=nothing)
    return () -> begin
        h = spin_weighted_spheroidal_harmonic(s, l, m, c; method, precision)
        Any[h.lambda; [h(θ, 0.7; theta_derivative=d) for θ in THETAS for d in 0:2]]
    end
end
spherical_job(s, l, m, method; precision=nothing) = () -> begin
    Y = spin_weighted_spherical_harmonic(s, l, m; method, precision)
    Any[Y(θ, 0.4; theta_derivative=d, phi_derivative=p) for θ in THETAS for d in 0:2 for p in 0:1]
end
eigenvalue_job(s, l, m, c; kwargs...) = () -> Any[spin_weighted_spheroidal_eigenvalue(s, l, m, c; kwargs...)]
continuation_job(s, l, m, c; precision=nothing) = () -> begin
    pair = continue_angular_mode(s, l, m, c; cache=AngularCache(), precision)
    path = track_angular_mode(s, l, m, [0, c / 2, c]; precision)
    Any[pair.lambda; pair.coefficients; [st.lambda for st in path.states]]
end

# bit for bit: === for Float64; for BigFloat the same value and precision
_same(p, q) = p === q || (p isa Number && q isa Number && typeof(p) == typeof(q) &&
    isequal(p, q) && (!(real(p) isa BigFloat) || precision(real(p)) == precision(real(q))))
same(x, y) = length(x) == length(y) && all(((p, q),) -> _same(p, q), zip(x, y))

function concurrent_matches_serial(jobs)
    serial = [_outputs(job) for job in jobs]
    for _ in 1:2
        tasks = [Threads.@spawn _outputs(job) for job in jobs]
        all(same.(fetch.(tasks), serial)) || return false
    end
    return true
end

@testset "Concurrent evaluation ($(Threads.nthreads()) threads)" begin
    float64_jobs = [
        harmonic_job(-2, 2, 2, 0.7, "auto"), harmonic_job(-2, 3, 1, 1.5, "chebyshev"),
        harmonic_job(-2, 2, 2, 0.7, "leaver"), harmonic_job(0, 4, -1, 0.3 - 0.4im, "auto"),
        harmonic_job(1, 31, 2, 2.0, "jacobi"), harmonic_job(-2, 2, 2, -20.0, "auto"),
        spherical_job(-2, 3, 1, "direct"), spherical_job(2, 33, -1, "jacobi"),
        spherical_job(-1, 2, 0, "chebyshev"),
        eigenvalue_job(-2, 2, 2, 3.0), eigenvalue_job(-2, 2, 2, 1.0 - 4.0im),
        eigenvalue_job(-2, 3, 2, 0.5; method="leaver"),
        continuation_job(-2, 3, 2, 0.4 + 0.15im), continuation_job(0, 2, 2, -5.35 + 2.288im),
    ]
    @test concurrent_matches_serial(float64_jobs)
    big_jobs = [
        harmonic_job(-2, 2, 2, 7//10, "auto"; precision=256), harmonic_job(-2, 2, 2, 7//10, "leaver"; precision=256),
        harmonic_job(0, 3, -1, 3//10 - 2im//5, "auto"; precision=256),
        spherical_job(-2, 3, 1, "direct"; precision=256), spherical_job(2, 31, -1, "jacobi"; precision=256),
        eigenvalue_job(-2, 2, 2, 3; precision=256), eigenvalue_job(-2, 2, 2, 1 - 4im; precision=256),
        continuation_job(-2, 3, 2, 2//5 + 3im//20; precision=256),
    ]
    @test concurrent_matches_serial(big_jobs)
    # the shared default cache, filled from several threads at once
    clear_angular_cache!()
    cs = [0.1im * k + 0.05 for k in 1:8]
    serial = [continue_angular_mode(-2, 2, 2, c; cache=nothing).lambda for c in cs]
    @test fetch.([Threads.@spawn continue_angular_mode(-2, 2, 2, c).lambda for c in cs]) == serial
end
