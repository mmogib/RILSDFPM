include("dependencies.jl")
include("included_files.jl")

settings = Dict(
    :dim => 100_000,
    :itrs => 2000,
    :param_size => 100,
    :tol => 1e-7,
)
ps = getAllTestProblem(settings[:dim])
# ; exclude=["NonsmoothSine","Logarithmic"]
dfs = pickRandomParams(RILSDFPMParameters, ps, settings[:itrs], settings[:param_size], settings[:tol])
vscodedisplay(dfs)

# (α, β, ξ, μ, ζ, γ, σ, ρ)
# params = RILSDFPMParameters(0.5, 0.8, 1.0, 0.5, 0.1, 1.3, 0.001, 0.46)
params = RILSDFPMParameters(0.24561134935215856, 0.8564877572100413, 0.20166231652151723, 0.3270752453777446, 0.8111659397644992, 0.6153130406009055, 0.7110051332386418, 1.7965)
# α = 0.24561134935215856, β = 0.8564877572100413, ξ = 0.20166231652151723, μ = 0.3270752453777446, ζ = 0.8111659397644992, γ = 0.6153130406009055, σ = 0.7110051332386418, ρ = 1.7965

options = AlgorithmOptions(3000, 1000, 1e-7)
dfs = runExperiment(RILSDFPMAlgorithm; dim=100_000, options=options, params=params, save=true)
vscodedisplay(dfs)