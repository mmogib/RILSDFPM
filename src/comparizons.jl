include("dependencies.jl")
include("included_files.jl")

dims = [1000, 10_000, 100_000]
options = AlgorithmOptions(2000, 1000, 1e-11, Inf, (fx0, fx1) -> norm(fx0 - fx1) / norm(fx0) <= 1e-7)
saveit = true
algorithms = [
    RILSDFPMAlgorithm,
    BBSTTDFPAlgorithm,
    BBSTTDFPAlgorithmNoIN,
    MOPCGAlgorithm,
    CGDFPAlgorithm
]



function runComparizon(dims::Vector{Int64}, options::AlgorithmOptions, algorithms::Vector{DataType}, saveit::Bool)
    dfs = runExperiment(algorithms, dims, options; saveit=saveit)

    pyplot()
    printstyled("Ploting ...\n", color=:blue)
    plts = plotResultsProfiles(algorithms, dfs)
    p1 = plot(plts..., layout=(3, 1), size=(900, 800))
    output_folder = outputfolder("./results/profiles")
    map(["png", "svg", "pdf", "eps"]) do ext
        savefig(p1, "$output_folder/profile.$ext")
    end
    printstyled("Saving plots in $output_folder \n", color=:blue)

    numericalResultDF2CSV(vcat(dfs...), "$output_folder/profile.csv")
    printstyled("Saving data in $output_folder/profile.csv \n", color=:blue)
end
# pgfplotsx()
# p2 = plot(plts..., layout=(3, 1), size=(800, 900), tex_output_standalone=true)
# savefig(p, "./results/latex/profiles.tex")
runComparizon(dims, options, algorithms, saveit)