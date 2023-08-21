PExponetialIII(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "PExponetialIII",
    x -> ExponetialIII(x),
    x -> projectOnBox(x; bounds=(0.0, Inf64)),
    x0
)

PPolynomialI(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "PPolynomialI",
    x -> PolynomialI(x),
    x -> projectOnBox(x; bounds=(0.0, Inf64)),
    x0
)


PSmoothSine(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "SmoothSine",
    x -> SmoothSine.(x),
    x -> projectOnBox(x; bounds=(0.0, Inf64)),
    x0
)





PPolynomialSineCosine(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "PolynomialSineCosine",
    x -> PolynomialSineCosine.(x),
    x -> projectOnBox(x; bounds=(0.0, Inf64)),
    x0
)

PExponetialI(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "PExponetialI",
    x -> ExponetialI.(x),
    x -> projectOnBox(x; bounds=(0.0, Inf64)),
    x0
)

PNonsmoothSine(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "NonsmoothSine",
    x -> NonsmoothSine.(x),
    x -> projectOnTriangle(x; lower=0),
    x0
)


PModifiedNonsmoothSine(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ModifiedNonsmoothSine",
    x -> ModifiedNonsmoothSine.(x),
    x -> projectOnTriangle(x; lower=-1),
    x0
)

PModifiedNonsmoothSine2(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ModifiedNonsmoothSine2",
    x -> ModifiedNonsmoothSine2.(x),
    x -> projectOnTriangle(x; lower=-1),
    x0
)

PExponetialSineCosine(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ExponetialSineCosine",
    x -> ExponetialSineCosine.(x),
    x -> projectOnTriangle(x; lower=-1),
    x0
)


PModifiedTrigI(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ModifiedTrigI",
    x -> ModifiedTrigI(x),
    x -> projectOnBox(x; bounds=(-3.0, Inf64)),
    x0
)


PModifiedTrigII(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ModifiedTrigII",
    x -> ModifiedTrigII(x),
    x -> projectOnBox(x; bounds=(-2.0, Inf64)),
    x0
)


# PTridiagonal(; x0::Vector{<:Real}=ones(100)) = TestProblem(
#     "Tridiagonal",
#     x -> Tridiagonal(x),
#     x -> x,
#     x0
# )


PModifiedTridiagonal(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ModifiedTridiagonal",
    x -> ModifiedTridiagonal(x),
    x -> projectOnTriangle(x; lower=0, β=1),
    x0
)



PLogarithmic(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "Logarithmic",
    x -> Logarithmic(x),
    x -> projectOnBox(x; bounds=(-1.0, Inf64)),
    x0
)

PNonmoothLogarithmic(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "NonmoothLogarithmic",
    x -> NonmoothLogarithmic(x),
    x -> projectOnBox(x; bounds=(0, Inf64)),
    x0
)


ARWHEADProblem(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ARWHEAD",
    x -> ARWHEADGrad(x),
    x -> x,
    x0
)

ARWHEADProblemProjectedOnBox(; x0::Vector{<:Real}=ones(100)) = TestProblem(
    "ARWHEADWithPrjection",
    x -> ARWHEADGrad(x),
    x -> projectOnTriangle(x; lower=0),
    x0
)
ENGVAL1Problem(; x0::Vector{<:Real}=2.0 * ones(100)) = TestProblem(
    "ENGVAL1",
    x -> ENGVAL1Grad(x),
    x -> x,
    x0
)

# t=0.026290993205301927, 
# β=0.6330801269981523,  
# σ=0.13316373378861357, 
# γ=1.449630720751143, 
# αmin=0.7680819467829025, 
# αmax=0.879578054884707, 
# r=0.38701943850338805, 
# ψ=0.1729093962921029

Ding2017PIIProblem(; dim=100) = TestProblem(
    "ding2017PII",
    x -> ding2017FunII.(x),
    x -> Pj3(x),
    ones(dim)
)

Ding2017DIAGONAL6Problem(; dim=100) = TestProblem(
    "DIAGONAL6",
    x -> ding2017Diagonal6Fun.(x),
    x -> x,
    ones(dim)
)
# t=0.026290993205301927, 
# β=0.6330801269981523,  
# σ=0.13316373378861357, 
# γ=1.449630720751143, 
# αmin=0.7680819467829025, 
# αmax=0.879578054884707, 
# r=0.38701943850338805, 
# ψ=0.1729093962921029

PENALTY1Problem(; dim=100) = TestProblem(
    "PENALTY1",
    x -> PENALTY1Grad(x),
    x -> x,
    1.0 * [i for i in 1:dim] # -ones(100)
)
# using Flux
# dixon3dq_flux(x) = gradient(DIXON3DQFun, x)[1]
DIXON3DQProblem(; dim=100) = TestProblem(
    "DIXON3DQ",
    x -> DIXON3DQGrad(x),
    x -> x,
    -ones(dim)
)

GENHUMPSProblem(; dim=100) = TestProblem(
    "GENHUMPS",
    x -> GENHUMPSGrad(x; ζ=2),
    x -> x,
    -[506.0, 506.2 * ones(dim - 1)...]
)



DIXMAANHProblem(; dim=300) = TestProblem(
    "DIXMAANH",
    x -> DIXMAANHGrad(x),
    x -> x,
    2.0 * ones(dim)
)


DIXMAANIProblem(; dim=300) = TestProblem(
    "DIXMAANI",
    x -> DIXMAANIGrad(x),
    x -> x,
    2.0 * ones(dim)
)






