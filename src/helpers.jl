"""
New functions for paper 2
"""
function getParametrizer()
    zdt = now(tz"Asia/Riyadh")
    t = parse(Int, Dates.format(zdt, "YYYYMMddHH"))
    rng = MersenneTwister(t)
    return rng
end
function getNewParams(::Type{RILSDFPMParameters}, sz::Int)
    rng = getParametrizer()
    params = map(1:sz) do _
        getNewParams(RILSDFPMParameters; rng=rng)
    end
    params
end
function getNewParams(::Type{RILSDFPMParameters}; rng::Union{Nothing,MersenneTwister}=nothing)
    rng = isnothing(rng) ? getParametrizer() : rng
    (α, β, ξ, μ, ζ, γ, σ) = rand(rng, 7)
    γ = rand(rng, [0, 1], 1)[1] + γ
    ρ_numerator = (16 / γ) * (α - 1)^2
    ρ_denominator = 16 * (α - 0.25)^2 + 7
    upper_bound_of_ρ = ρ_numerator / ρ_denominator
    ρ = rand(rng, 0:0.0001:upper_bound_of_ρ)
    RILSDFPMParameters(α, β, ξ, μ, ζ, γ, σ, ρ)
end

function pickRandomParams(
    ::Type{RILSDFPMParameters},
    ps::Vector{TestProblem},
    itrs::Int64,
    szs::Int64,
    tol::Float64;
    rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM")))
)
    println("Running Experiments to pick best parameters for RILSDFPMParameters.")

    params = getNewParams(RILSDFPMParameters, szs)
    options = AlgorithmOptions(itrs, tol)
    dfs = map(enumerate(params)) do (pindx, param)
        pdfs = map(ps) do problem
            println("Solving ... $(problem.name)  parameter set: $pindx out of $szs")
            nr = try

                result, elapsed_time = @timed RILSDFPM(problem.f, problem.x0, param, options)
                NumericalResult(problem, length(problem.x0), options, param, SuccessfullResult(result, elapsed_time))
            catch err
                NumericalResult(problem, length(problem.x0), options, param, FailedResult("Opps $(problem.name) -- $err"))
            end
            numericalResult2DataFrame(nr, algorithm="RILSDFPM")
        end
        vcat(pdfs...)
    end
    # |> sols -> filter(sol -> sol.result.flag != :maxiter_reached, sols)

    df = vcat(dfs...)

    println("Saving to file ... ")
    file_name = outputfilename("paramters", "csv")
    numericalResultDF2CSV(df, file_name)
    println("DONE!")
    return df
end


function runExperiment(::Type{RILSDFPMAlgorithm}; dim::Int64, options::AlgorithmOptions, params::RILSDFPMParameters, save::Bool=false)
    ps = getAllTestProblem(dim)
    println("Running our RILSDFPM Algorithm")
    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("RILSDFPMAl ($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        RILSDFPM(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="RILSDFPM")


    df = vcat(dfs...)
    if save
        filename = outputfilename("RILSDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end



function runExperiment(::Type{RILSDFPMAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)
    ps = createTestProlems(dim)
    params = getAlgorithmParams(RILSDFPMAlgorithm, experiment_name="numerical")
    println("Running our new RILSDFPM Algorithm")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("RILSDFPM($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        RILSDFPM(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="RILSDFPM")


    df = vcat(dfs...)
    if save
        filename = outputfilename("RILSDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end


"""
Old function from paper 1

"""


function cleanVectorFromEmpty(x, v::Float64)
    x .|> d -> isempty(d) ? v : d
end
function createTestProlem(name::String, F::Function, P::Function, Pcheck::Function, x0::Vector{<:Real}, x0label::String)
    TestProblem(
        name,
        F,
        P,
        Pcheck,
        x0,
        x0label
    )
end

function createTestProlem(fns::Vector{Tuple{String,Function,Function,Function}}, x::Tuple{String,Vector{Float64}})
    map(fns) do fn
        name, F, P, Pcheck = fn
        x0label, x0 = x
        createTestProlem(name, F, P, Pcheck, x0, x0label)
    end
end

createParams(::Type{T}) where {T<:GenericProblem} = createParams(T, 1000)
function createParams(::Type{BBSTTDFPProblem}, sz::Int64, rng::MersenneTwister)
    ϵ = eps(1.0)
    ηs = rand(rng, ϵ:ϵ:0.5, sz, 1)
    ξs = ηs .+ rand(rng, ϵ:ϵ:(0.5-ϵ), sz, 1)
    hcat(ηs, ξs)
end
function createParams(::Type{BBSTTDFPProblem},
    sz::Int64;
    rng::MersenneTwister=MersenneTwister(parse(Int, Dates.format(now(), "ddHHMM"))),
    newparams::Bool=false
)
    rand_init = rand(rng, sz, 8)
    rand_init[:, 4] = rand_init[:, 4] + rand(rng, sz)
    rand_init[:, 6] = rand_init[:, 5] + rand(rng, sz)
    rand_init = newparams ? hcat(rand_init, createParams(BBSTTDFPProblem, sz, rng)) : rand_init
    params = map(p -> BBSTTDFPParams(p...), eachrow(rand_init))
    params
end

function numericalResult2DataFrame(result::NumericalResult; algorithm::String="")
    algorithm = algorithm
    problem_name = result.problem.name
    x0_label = result.problem.x0label
    dim = result.dim
    if (isa(result.result, SuccessfullResult))
        iters = result.result.iterations
        evals = result.result.functionEvalauations
        fNorm = @sprintf("%.5e", result.result.Fnorm)
        exe_time = result.result.exec_time
        message = isnothing(result.result.flag) ? "problem solved" : result.result.flag
        return DataFrame(
            algorithm=[algorithm],
            problem_name=[problem_name],
            x0=[x0_label],
            dim=[dim],
            iters=[iters],
            evals=[evals],
            time=isnothing(exe_time) ? ["--"] : [exe_time],
            fnorm=[result.result.Fnorm],
            fNorm=[fNorm],
            message=[message],
            options=[result.options],
            params=[result.params]
        )

    else
        return DataFrame(
            algorithm=[algorithm],
            problem_name=[problem_name],
            x0=[x0_label],
            dim=[dim],
            iters=[""],
            evals=[""],
            time=[""],
            fnorm=[""],
            fNorm=[""],
            message=["❌😒" * result.result.message],
            options=[result.options],
            params=[result.params]
        )
    end
end


function numericalResult2DataFrame(results::Vector{NumericalResult}, algorithm::String="")
    df = DataFrame()
    dfs = results .|> result -> numericalResult2DataFrame(result, algorithm=algorithm)
    vcat(df, dfs...)
end

function numericalResultDF2CSV(df::DataFrame, filename::String)
    df |> CSV.write(filename)
end



function optimizeParameters(p::TestProblem)
    optimizeParameters([p])
end
function optimizeParameters(ps::Vector{TestProblem})
    function goodParams(vparams::Vector{Float64}, ps::Vector{TestProblem})
        params = BBSTTDFPParams(vparams...)
        options = AlgorithmOptions(50_000, 1e-5)

        results = ps .|> p -> BBSTTDFP(p; params=params, options=options)
        if all(map(result -> isa(result.result, SuccessfullResult), results))
            return max(map(result -> result.result.iterations, results)...)
        else
            return options.maxiters
        end
    end
    names_of_all = ps .|> (p -> p.name) |> names -> join(names, ",")
    println("Working on problems $(names_of_all), please wait..")
    ga = GA(populationSize=100,
        selection=susinv,
        crossover=DC,
        mutation=PLM()
    )
    lower = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0002, 0.0001, 0.0001]
    upper = [0.9999, 0.9999, 0.9999, 1.9999, 0.9999, 0.9999, 0.9999, 0.9999]
    result = Evolutionary.optimize(
        x -> goodParams(x, ps),
        BoxConstraints(lower, upper),
        ones(8),
        ga
    )
    root_dir = outputfolder()
    params = BBSTTDFPParams(Evolutionary.minimizer(result)...)
    options = AlgorithmOptions(50000, 1e-5)
    for p in ps
        println("Solving problem $p.name with optimized parameters...")
        result_p = BBSTTDFP(p; params=params, options=options)
        dfs = numericalResult2DataFrame(result_p)
        println("Saving to file ... ")
        numericalResultDF2CSV(dfs, "$root_dir/optimized_params_$(p.name).csv")
    end
end


function getNewParams(oldparams::BBSTTDFPParams, paramsToRenew::Vector{Symbol})
    t, β, σ, γ, αmin = oldparams.t, oldparams.β, oldparams.σ, oldparams.γ, oldparams.αmin
    αmax, r, ψ, η, ξ = oldparams.αmax, oldparams.r, oldparams.ψ, oldparams.η, oldparams.ξ
    # t::Float64 ∈ (0,1) ~ 0.11
    t = (:t in paramsToRenew) ? rand() : t
    # β::Float64 ∈ (0,1) ~ 0.5
    β = (:β in paramsToRenew) ? rand() : β
    # σ::Float64 ∈ (0,1) ~ 0.01
    σ = (:σ in paramsToRenew) ? rand() : σ
    # γ::Float64 ∈ (0,2) ~ 1.8 or 1.72
    γ = (:γ in paramsToRenew) ? rand([0, 1]) + rand() : γ
    # αmin::Float64 0 < αmin < αmax < 1 ~ 0
    αmin = (:αmin in paramsToRenew) ? rand() : αmin
    # αmax::Float64 0 < αmin < αmax < 1 ~ Inf64
    αmax = (:αmax in paramsToRenew) ? αmin + rand() : αmax
    # r::Float64  ∈ ( 0,1) ~ 0.1
    r = (:r in paramsToRenew) ? rand() : r
    # ψ::Float64 ∈ (0 ,1) ~ 0.5 
    ψ = (:ψ in paramsToRenew) ? rand() : ψ
    # η::Union{Nothing,Float64} optional 0 < η <= ξ  ~ 0.001 == μ in a2023LiuWuShaoZhangCao
    η = (:η in paramsToRenew) ? rand() : η
    # ξ::Union{Nothing,Float64} optional 0 < η <= ξ  ~ 0.6 == ν in a2023LiuWuShaoZhangCao
    ξ = (:ξ in paramsToRenew) ? rand() : ξ


    BBSTTDFPParams(t, β, σ, γ, αmin, αmax, r, ψ, η, ξ)
end

function str2params(strparams::String, rng::MersenneTwister;
    randomized::Union{Nothing,Vector{Symbol}}=nothing,
    newparams::Bool=false
)

    parts = Dict(split(strparams, ",") .|> x -> split(x, "=") .|> x -> strip.(x))

    namedtuples = (; zip(Symbol.(keys(parts)), parse.(Float64, (values(parts))))...)


    getnewparam(s::Symbol, passedparams::NamedTuple, randomized::Union{Nothing,Vector{Symbol}}) = begin
        if isnothing(randomized)
            get(passedparams, s, rand(rng))
        else
            s in randomized ? rand(rng) : get(passedparams, s, rand(rng))
        end
    end
    pms = BBSTTDFPParams(
        getnewparam(:t, namedtuples, randomized),
        getnewparam(:β, namedtuples, randomized),
        getnewparam(:σ, namedtuples, randomized),
        getnewparam(:γ, namedtuples, randomized),
        getnewparam(:αmin, namedtuples, randomized),
        getnewparam(:αmax, namedtuples, randomized),
        getnewparam(:r, namedtuples, randomized),
        getnewparam(:ψ, namedtuples, randomized),
    )
    npms = if newparams
        ϵ = eps(1.0)
        η = rand(rng, ϵ:ϵ:0.2)
        ξ = rand(rng, (η+0.3):ϵ:(1-ϵ))
        BBSTTDFPParams(pms.t,
            pms.β,
            pms.σ,
            pms.γ,
            pms.αmin,
            pms.αmax,
            pms.r,
            pms.ψ,
            η,
            ξ
        )
    else
        pms
    end
    npms
end

function outputfolder(root::String="./results")
    zdt = now(tz"Asia/Riyadh")
    dayfolder = Dates.format(zdt, "yyyy_mm_dd")
    hourfolder = Dates.format(zdt, "HH_MM")
    root_dir = mkpath("$root/$dayfolder/$hourfolder")
    return root_dir
end

function outputfilename(name::String, extension::String; root::String="./results", suffix::Union{Nothing,String}=nothing)
    root_dir = outputfolder(root)
    filename = if isnothing(suffix)
        "$root_dir/$name.$extension"
    else
        "$root_dir/$(name)_$suffix.$extension"
    end
    filename
end



function runExperiment(::Type{MOPCGAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    ps = createTestProlems(dim)
    params = getAlgorithmParams(MOPCGAlgorithm, experiment_name="numerical")
    println("Running MOPCG - SabiuShahStanimirovicIvanovWaziri (2023) Experiments for comparing.")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("MOPCG($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        MOPCG(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="MOPCG")


    df = vcat(dfs...)
    if save
        filename = outputfilename("MOPCG", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end



function runExperiment(::Type{CGDFPAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)

    ps = createTestProlems(dim)
    params = getAlgorithmParams(CGDFPAlgorithm, experiment_name="numerical")
    println("Running CGDFP Algorithm - 2020ZhengYangLiang (2020) Experiments for comparing.")

    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("CGDFP($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        CGDFP(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="CGDFPM")


    df = vcat(dfs...)
    if save
        filename = outputfilename("CGDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df

end


function runExperiment(::Type{AHDFPMAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)
    ps = createTestProlems(dim)
    params = getAlgorithmParams(AHDFPMAlgorithm, experiment_name="numerical")

    println("Running AHDFPM Algorithm Algorithm - 2023LiuWuShaoZhangCao (2023) Experiments for comparing.")
    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("AHDFPM($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        AHDFPM(problem; options=options, params=params)
    end |> res -> numericalResult2DataFrame(res, algorithm="AHDFPM")


    df = vcat(dfs...)
    if save
        filename = outputfilename("AHDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end
    return df
end

# 

function runExperiment(::Type{BBSTTDFPAlgorithmNoIN}, dim::Int64, options::AlgorithmOptions; save::Bool=false)
    ps = createTestProlems(dim)
    params = getAlgorithmParams(BBSTTDFPAlgorithmNoIN, experiment_name="numerical", ϵ=options.tol)
    println("Running our STTDFPM  Algorithm")
    counter = 1
    total = length(ps)

    dfs = ps .|> problem -> begin
        println("STTDFPM ($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        BBSTTDFP(problem, params=params, options=options; inertial=false)

    end |> res -> numericalResult2DataFrame(res, algorithm="STTDFPM")
    df = vcat(dfs...)
    if save
        filename = outputfilename("STTDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df
end

function runExperiment(::Type{BBSTTDFPAlgorithm}, dim::Int64, options::AlgorithmOptions; save::Bool=false)
    ps = createTestProlems(dim)
    params = getAlgorithmParams(BBSTTDFPAlgorithm, experiment_name="numerical", ϵ=options.tol)
    println("Running our ISTTDFPM Algorithm")
    counter = 1
    total = length(ps)
    dfs = ps .|> problem -> begin
        println("ISTTDFPM ($(length(problem.x0))) - ($counter/$total) Solving $(problem.name)...")
        counter += 1
        BBSTTDFP(problem, params=params, options=options)

    end |> res -> numericalResult2DataFrame(res, algorithm="ISTTDFPM")

    df = vcat(dfs...)
    if save
        filename = outputfilename("ISTTDFPM", "csv")
        numericalResultDF2CSV(df, filename)
    end

    return df
end


function runExperiment(algorithms::Vector{T} where {T<:DataType}, dims::Vector{Int64}, options::AlgorithmOptions;
    saveit::Bool=false)
    ## [*] Run our experiments  
    alldfs = map(dims) do dim
        df = map(algorithms) do algo
            runExperiment(algo, dim, options, save=saveit)
        end
        df
    end
    # dfs = 

    vcat(alldfs...)

end
function plotResultsProfiles(algorithms::Vector{T} where {T<:DataType}, dfs::Vector{DataFrame}; k=11)

    labels = algorithms .|> t -> replace("$t", "Algorithm" => "")
    # labels = collect(1:length(algorithms)) .|> i -> algorithms[i][2] ? "$(labels[i])_NEW_PARAMS" : labels[i]
    # xlabel = join(labels, " vs ")
    xlabel = L"\tau"
    legendfontsize, leg = 8, :bottomright
    plts = map([
        (:iters, "Iterations"),
        (:evals, "Evalauations"),
        (:time, "CPU Time")]) do (item, title)
        T = map(dfs) do df
            values = df[!, item]
            val = Float64(maximum(values[findall(!isempty, values)]))
            values |> x -> cleanVectorFromEmpty(x, 2 * val) .|> Float64
        end
        theme(:default)
        performance_profile(
            PlotsBackend(),
            hcat(T...),
            labels,
            title=title,
            xlabel=xlabel,
            ylabel=L"\rho(\tau)",
            legendfontsize=legendfontsize,
            leg=leg,
            palette=:Dark2_5
        )
    end
    plts

end
function getAllTestProblem(dim::Int64; exclude::Vector{String})
    getAllTestProblem(dim) |> ps -> filter(p -> !(p.name ∈ exclude), ps)
end
function getAllTestProblem(dim::Int64)
    # (name::String, F::Function, P::Function, Pcheck::Function, x0::Vector{<:Real}, x0label::String)
    problem2 = createTestProlem(
        "PolynomialI",
        x -> PolynomialI(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        ones(dim),
        "ones"
    )
    problem3 = createTestProlem(
        "SmoothSine",
        x -> SmoothSine.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        ones(dim),
        "ones"
    )

    # problem4 = createTestProlem(
    #     "PolynomialSineCosine",
    #     x -> PolynomialSineCosine.(x),
    #     x -> projectOnBox(x; bounds=(0.0, Inf64)),
    #     x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
    #     ones(dim),
    "ones"
    # )
    problem5 = createTestProlem(
        "ExponetialI",
        x -> ExponetialI.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        ones(dim),
        "ones"
    )
    problem6 = createTestProlem(
        "NonsmoothSine",
        x -> NonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=0),
        x -> projectOnTriangleCheck(x; lower=0),
        ones(dim),
        "ones"
    )
    problem7 = createTestProlem(
        "ModifiedNonsmoothSine",
        x -> ModifiedNonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        ones(dim),
        "ones"
    )
    problem8 = createTestProlem(
        "ModifiedNonsmoothSine2",
        x -> ModifiedNonsmoothSine2.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        ones(dim),
        "ones"
    )
    problem9 = createTestProlem(
        "ExponetialSineCosine",
        x -> ExponetialSineCosine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        ones(dim),
        "ones"
    )
    problem10 = createTestProlem(
        "ModifiedTrigI",
        x -> ModifiedTrigI(x),
        x -> projectOnBox(x; bounds=(-3.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-3.0, Inf64)),
        ones(dim),
        "ones"
    )
    # problem11 = createTestProlem(
    #     "ModifiedTrigII",
    #     x -> ModifiedTrigII(x),
    #     x -> projectOnBox(x; bounds=(-2.0, Inf64)),
    #     dim
    # ones(# )),
    # "ones"
    problem12 = createTestProlem(
        "ModifiedTridiagonal",
        x -> ModifiedTridiagonal(x),
        x -> projectOnTriangle(x; lower=0, β=1),
        x -> projectOnTriangleCheck(x; lower=0, β=1),
        ones(dim),
        "ones"
    )
    problem13 = createTestProlem(
        "Logarithmic",
        x -> Logarithmic(x),
        x -> projectOnBox(x; bounds=(-1.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-1.0, Inf64)),
        ones(dim),
        "ones"
    )
    problem14 = createTestProlem(
        "NonmoothLogarithmic",
        x -> NonmoothLogarithmic(x),
        x -> projectOnBox(x; bounds=(0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0, Inf64)),
        ones(dim),
        "ones"
    )

    # problem15 = createTestProlem(
    #     "ARWHEADGrad",
    #     x -> ARWHEADGrad(x),
    #     x -> projectOnTriangle(x; lower=0),
    #     x -> projectOnTriangleCheck(x; lower=0),
    #     ones(dim),
    # "ones"
    # )

    problem16 = createTestProlem(
        "ENGVAL1Grad",
        x -> ENGVAL1Grad(x),
        x -> x,
        x -> true,
        ones(dim),
        "ones"
    )

    pseudo1 = createTestProlem(
        "Pseudo1",
        x -> PseudoMonotone1(x),
        x -> projectOnBox(x; bounds=(-10, 10)),
        x -> true,
        ones(2),
        "ones"
    )
    pseudo2 = createTestProlem(
        "Pseudo2",
        x -> PseudoMonotone2(x),
        x -> projectOnBox(x; bounds=(-5, 5)),
        x -> true,
        ones(3),
        "ones"
    )

    return vcat(
        pseudo1,
        pseudo2,
        problem2,
        problem3,
        # problem4,
        problem5,
        problem6,
        problem7,
        problem8,
        problem9,
        problem10,
        problem12,
        problem13,
        problem14,
        problem16
    )

end
function createTestProlems(dim::Int64)
    # PExponetialIII(x0=x0),
    # problem1 = createTestProlems(
    #     "ExponetialIII",
    #     x -> ExponetialIII(x),
    #     x -> projectOnBox(x; bounds=(0.0, Inf64)),
    #     x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
    #     dim
    # )
    problem2 = createTestProlems(
        "PolynomialI",
        x -> PolynomialI(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )
    problem3 = createTestProlems(
        "SmoothSine",
        x -> SmoothSine.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )

    # problem4 = createTestProlems(
    #     "PolynomialSineCosine",
    #     x -> PolynomialSineCosine.(x),
    #     x -> projectOnBox(x; bounds=(0.0, Inf64)),
    #     x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
    #     dim
    # )
    problem5 = createTestProlems(
        "ExponetialI",
        x -> ExponetialI.(x),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        dim
    )
    problem6 = createTestProlems(
        "NonsmoothSine",
        x -> NonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=0),
        x -> projectOnTriangleCheck(x; lower=0),
        dim
    )
    problem7 = createTestProlems(
        "ModifiedNonsmoothSine",
        x -> ModifiedNonsmoothSine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem8 = createTestProlems(
        "ModifiedNonsmoothSine2",
        x -> ModifiedNonsmoothSine2.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem9 = createTestProlems(
        "ExponetialSineCosine",
        x -> ExponetialSineCosine.(x),
        x -> projectOnTriangle(x; lower=-1),
        x -> projectOnTriangleCheck(x; lower=-1),
        dim
    )
    problem10 = createTestProlems(
        "ModifiedTrigI",
        x -> ModifiedTrigI(x),
        x -> projectOnBox(x; bounds=(-3.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-3.0, Inf64)),
        dim
    )
    # problem11 = createTestProlems(
    #     "ModifiedTrigII",
    #     x -> ModifiedTrigII(x),
    #     x -> projectOnBox(x; bounds=(-2.0, Inf64)),
    #     dim
    # )
    problem12 = createTestProlems(
        "ModifiedTridiagonal",
        x -> ModifiedTridiagonal(x),
        x -> projectOnTriangle(x; lower=0, β=1),
        x -> projectOnTriangleCheck(x; lower=0, β=1),
        dim
    )
    problem13 = createTestProlems(
        "Logarithmic",
        x -> Logarithmic(x),
        x -> projectOnBox(x; bounds=(-1.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(-1.0, Inf64)),
        dim
    )
    problem14 = createTestProlems(
        "NonmoothLogarithmic",
        x -> NonmoothLogarithmic(x),
        x -> projectOnBox(x; bounds=(0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0, Inf64)),
        dim
    )

    # problem15 = createTestProlems(
    #     "ARWHEADGrad",
    #     x -> ARWHEADGrad(x),
    #     x -> projectOnTriangle(x; lower=0),
    #     x -> projectOnTriangleCheck(x; lower=0),
    #     dim
    # )

    problem16 = createTestProlems(
        "ENGVAL1Grad",
        x -> ENGVAL1Grad(x),
        x -> x,
        x -> true,
        dim
    )

    return vcat(
        problem2,
        problem3,
        # problem4,
        problem5,
        problem6,
        problem7,
        problem8,
        problem9,
        problem10,
        problem12,
        problem13,
        problem14,
        problem16
    )

end
function createTestProlems(name::String, F::Function, P::Function, Pcheck::Function, dim::Int64)
    onez = ones(dim)
    indxs = 1:dim
    startin_points_1 = [
        ("x1: zeros", zeros(dim)),
        ("x2: 0.2", 0.2 * onez),
        ("x3: 0.4", 0.4 * onez),
        ("x4: 0.5", 0.5 * onez),
        ("x5: 0.6", 0.6 * onez),
        ("x6: 0.8", 0.8 * onez),
        # ("x6: -1", -1 * onez),
        ("x7: 1", onez),
        ("x8: 1.1", 1.1 * onez),
    ]
    startin_points_2 = [
        # ("x9: 1/2ᵏ", indxs .|> t -> 1 / 2^t),
        ("x10: 1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
        ("x11: 1/k", indxs .|> t -> 1 / t),
        ("x12: (k-1)/n", indxs .|> t -> (t - 1) / dim),
        ("x13: 1/n", indxs .|> t -> 1 / dim),
        ("x14: 1/3ᵏ", indxs .|> t -> 1 / BigFloat(3^t)),
        ("x15: k/n", indxs .|> t -> t / dim)
    ]
    startin_points = [
        startin_points_1..., startin_points_2...
    ]
    # startin_points = [
    #     ("x1: 1", onez),
    #     ("x2: 0.2", 0.2 * onez),
    #     ("x3: 1/2ᵏ", indxs .|> t -> 1 / 2^t),
    #     ("x4: (k-1)/n", indxs .|> t -> (t - 1) / dim),
    #     ("x5: 1/k", indxs .|> t -> 1 / t),
    #     ("x6: 1/n", indxs .|> t -> 1 / dim),
    #     ("x7: 1 - 1/n", indxs .|> t -> 1 - (1 / dim)),
    #     ("x8: 1.1", 1.1 * onez),
    # ]
    map(startin_points) do (label, x0)
        createTestProlem(name, F, P, Pcheck, x0, label)
    end
end


"""
CS:Compressed Sensing helpers
"""

function createCSData(m, n, r, σ; d::Distribution=Normal())
    x = zeros(n)
    q = randperm(n)
    x[q[1:r]] = rand(d, r)
    A = rand(d, m, n)
    A = Matrix(qr(A').Q)'
    noise = σ * rand(d, m)
    b = A * x + noise
    x0 = A' * b
    τ = 0.01 * norm(x0, Inf)
    c = τ * ones(2n) + vcat(-x0, x0)
    z0 = vcat(max.(x0, 0), max.(-x0, 0))
    A' * A, z0, c, x, b
end

function createCSPlots(original, observed, reconstructed)
    p1 = bar(original, ylims=(-1.1, 1.1), bar_width=0.001, title="Original Signal", label=nothing)
    p2 = plot(observed, title="Observed Signal", label=nothing)
    p3 = bar(reconstructed, ylims=(-1.1, 1.1), bar_width=0.001, title="Reconstrcucted Signal", label=nothing)
    p4 = plot(p1, p2, p3, layout=(3, 1), size=(900, 300))
    p1, p2, p3, p4
end

# l₁ Regularized Least Squares 
function L1LS(z, ATA, c, n)
    u = z[1:n]
    v = z[n+1:2n]
    # z = vcat(u, v)
    Bu = ATA * (u - v)
    Hz = vcat(Bu, -Bu)
    min.(z, Hz + c)
end
function createCSProblems(ATA, z0, c, n)

    problem1 = TestProblem(
        "SignalProcessing - no projection",
        x -> L1LS(x, ATA, c, n),
        x -> x,
        x -> true,
        z0,
        "A'b"
    )
    problem2 = TestProblem(
        "SignalProcessing - with projection",
        x -> L1LS(x, ATA, c, n),
        x -> projectOnBox(x; bounds=(0.0, Inf64)),
        x -> projectOnBoxCheck(x; bounds=(0.0, Inf64)),
        z0,
        "A'b"
    )
    [
        problem1,
        problem2
    ]
end
function saveDataToFile(flder_to_save::String,
    ATA, initial_point, c, original_signal, observed_signal)

    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_ATA.txt", "w") do io
        writedlm(io, ATA)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_initial_point.txt", "w") do io
        writedlm(io, initial_point)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_c.txt", "w") do io
        writedlm(io, c)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_original_signal.txt", "w") do io
        writedlm(io, original_signal)
    end
    open("$flder_to_save/matrix_$(m)_$(n)_$(r)_observed_signal.txt", "w") do io
        writedlm(io, observed_signal)
    end

end

function getSavedData(flder_where_saved::String, m, n, r)
    ATA = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_ATA.txt")
    initial_point = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_initial_point.txt")
    c = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_c.txt")
    original_signal = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_original_signal.txt")
    observed_signal = readdlm("$flder_where_saved/matrix_$(m)_$(n)_$(r)_observed_signal.txt")
    ATA, initial_point, c, original_signal, observed_signal
end

function MSI(original, reconstructed)
    (1 / length(original)) * norm(reconstructed - original)^2
end


function PSNR(original, reconstructed)
    10 * log10(norm(reconstructed, Inf)^2 / MSI(original, reconstructed))
end
function readPaperData(dfs::DataFrame, names::Vector{String}; newnames::Union{Nothing,Vector{String}}=nothing)
    new_names = isnothing(newnames) ? names : newnames
    dfs = transform(dfs, :algorithm => ByRow(strip) => :algorithm)
    data = map(enumerate(names)) do (i, name)
        filter(row -> row[:algorithm] == name, dfs) |> df -> transform(df, :algorithm => ByRow(_ -> new_names[i]) => :algorithm)
    end
    data

end
function readPaperData(input_file::String, names::Vector{String}; newnames::Union{Nothing,Vector{String}}=nothing)
    dfs = CSV.File(input_file) |> DataFrame
    readPaperData(dfs, names; newnames=newnames)
end
function plotPerformanceProfile(T::Matrix{<:Number}, title::String, labels::Vector{String}; colors::Union{ColorPalette,Vector{Symbol}}=palette(:tab10), logscale::Bool=true)
    (w, h) = Plots._plot_defaults[:size]
    (_, _, max_ratio) = performance_profile_data(T)
    p = performance_profile(PlotsBackend(),
        T, labels;
        size=(1.2w, h),
        logscale,
        title=title,
        xlabel=L"\tau",
        ylabel=L"\rho(\tau)",
        legendfontsize=8,
        linestyle=:dash,
        palette=colors,
        linewidth=2.5,
        minorgrid=true, leg=:bottomright
    )
    p = if (logscale)
        xs = 1:ceil(max_ratio + 0.5)
        ys = 0:0.1:1.01
        plot(p,
            xticks=(xs, map(x -> "$(Int(x))", xs)),
            yticks=(ys, map(x -> "$x", ys)),
            framestyle=:origin
        )
    else
        p
    end
    p

end

function producePaperProfileImages(input_file::String)
    colors = [:red, :black, :green, :blue, :purple]

    output_folder = "./results/paper/imgs"

    newnames = [
        "RILSDFPM",
        "ISTTDFPM",
        "STTDFPM",
        "MOPCG",
        "CGDFPM"
    ]
    printstyled("reading stored data in "; color=:blue)
    printstyled("$input_file\n"; color=:blue, bold=true)
    data = readPaperData(input_file, newnames)
    printstyled("creating plots ...\n"; color=:green)


    ns = length(data)
    np, = size(data[1])
    # 
    plts = map([
        (:iters, "Iterations"),
        (:evals, "Function Evaluations"),
        (:time, "CPU Time")]) do (item, title)
        T = dataFrames2Matrix(data, item, np, ns)
        printstyled("creating plots for $title \n"; color=:green)
        p = plotPerformanceProfile(T, title, newnames; colors)
        p, title
    end
    map(plts) do (p, title)
        file_name = replace(title, " " => "_") |> lowercase
        printstyled("saving plots for $title to file \n"; color=:reverse)
        png_file = outputfilename(file_name, "png"; root=output_folder, suffix="performance_profiles")
        savefig(p, png_file)
        svg_file = outputfilename(file_name, "svg"; root=output_folder, suffix="performance_profiles")
        savefig(p, svg_file)
        pdf_file = outputfilename(file_name, "pdf"; root=output_folder, suffix="performance_profiles")
        savefig(p, pdf_file)
        eps_file = outputfilename(file_name, "eps"; root=output_folder, suffix="performance_profiles")
        savefig(p, eps_file)
    end
    printstyled("Finished saving in $output_folder\n"; color=:green)
    plts[1]
end
function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol, np::Int64, ns::Int64)

    T = map(data) do df
        values = df[!, field]
        values .|> x -> isempty(x) ? NaN : Float64(x)
    end |> d -> vcat(d...) |> d -> reshape(d, np, ns)
    T
end

function dataFrames2Matrix(data::Vector{DataFrame}, field::Symbol)
    ns = length(data)
    np, = size(data[1])
    dataFrames2Matrix(data, field, np, ns)
end

