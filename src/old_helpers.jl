
# val = Float64(maximum(values[findall(!isempty, values)]))
# values |> x -> cleanVectorFromEmpty(x, 2 * val) .|> Float64



# function runNumericalExperiments(::Type{BBSTTDFPProblem},
#     ps::Vector{TestProblem},
#     params::BBSTTDFPParams,
#     options::AlgorithmOptions;
#     save::Bool=false,
#     experiment_name::Union{Nothing,String}=nothing
# )

#     println("Running Experiments $(isnothing(experiment_name) ? "" : "($experiment_name)") testing new algorithm.")
#     dfs = ps .|> problem -> begin
#         println("Solving ... $(problem.name)")
#         results = BBSTTDFP(problem; params, options=options)
#         numericalResult2DataFrame(results, algorithm="BBSTTDFPAlgorithm")
#     end
#     df = vcat(dfs...)
#     if (save)
#         println("Saving to file ... ")
#         filename = isnothing(experiment_name) ? "" : replace(experiment_name, " " => "_")
#         file_name = outputfilename(filename, "csv")
#         numericalResultDF2CSV(df, file_name)
#         println("DONE!")
#     end
#     df
# end
# function runNumericalExperiments(
#     ::Type{BBSTTDFPParams},
#     problems::Vector{TestProblem};
#     save::Bool=false,
#     experiments_count::Int64=5,
#     paramstr::Vector{String}=["t=0.13824781669963615, β=0.5909641129216927,  σ=0.25514946641464453, γ=1.4117390335497881, αmin=0.912323896610999, αmax=1.693242503065831, r=0.8025450374268133, ψ=0.6642550140797701"],
#     randomized::Union{Nothing,Vector{Symbol}}=nothing,
#     options::AlgorithmOptions=AlgorithmOptions(100_000, 1e-5),
#     rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
#     newparams::Bool=false
# )
#     params = Iterators.flatten(str2params.(paramstr, rng; randomized=randomized, newparams=newparams) for _ in 1:(experiments_count)) |> collect


#     dfs = enumerate(params) .|> p -> begin
#         i, param = p
#         runNumericalExperiments(BBSTTDFPProblem, problems, param, options; experiment_name="experiment $i")
#     end
#     df = vcat(dfs...)
#     if save
#         file_name = outputfilename("testing_bbsttdp", "csv")
#         numericalResultDF2CSV(df, file_name)
#     end
#     dfs
# end

# function runNumericalExperiments(
#     ::Type{BBSTTDFPParams},
#     ps::Vector{TestProblem},
#     itrs::Int64,
#     szs::Int64;
#     rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
#     newparams::Bool=false
# )
#     println("Running Experiments to pick best parameters.")

#     params = createParams(BBSTTDFPProblem, szs; rng=rng, newparams=newparams)
#     options = AlgorithmOptions(itrs, 1e-5)
#     dfs = enumerate(ps) .|> p -> begin
#         i, problem = p
#         println("Solving ... $(problem.name)")
#         results = BBSTTDFP(problem, params; options=options)

#         numericalResult2DataFrame(results)
#     end
#     df = vcat(dfs...)

#     println("Saving to file ... ")
#     file_name = outputfilename("paramters", "csv")
#     numericalResultDF2CSV(df, file_name)
#     println("DONE!")
#     return df
# end
