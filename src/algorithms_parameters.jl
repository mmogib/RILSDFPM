function getAlgorithmParams(::Type{BBSTTDFPAlgorithm}; experiment_name::String="numerical", ϵ::Float64=1e-30)
    getAlgorithmParams(BBSTTDFPAlgorithmNoIN; experiment_name=experiment_name, ϵ=ϵ)
end

function getAlgorithmParams(::Type{BBSTTDFPAlgorithmNoIN}; experiment_name::String="numerical", ϵ::Float64=1e-30)
    if experiment_name == "numerical"

        params_numerical_exprs = BBSTTDFPParams(
            0.11, # 0.98
            0.5,
            0.01,
            1.8,
            10 * ϵ,# 1e-30,
            Inf64, #1e+30,
            0.1,
            0.2,#0.5
            0.001,
            0.6
        )
        return params_numerical_exprs

    elseif experiment_name == "cs"
        #used for cs_proccessing
        params2 = BBSTTDFPParams(
            0.9,
            0.85,
            0.9,
            1.75,
            10 * ϵ,
            Inf,
            0.94,
            0.1,
            0.001,
            0.6
        )

        return params2
    end
end


function getAlgorithmParams(::Type{AHDFPMAlgorithm}; experiment_name::String="numerical")
    experiment_name == "numerical" ? AHDFPMParameters(0.6, 1.72, 0.001, 0.001, 0.6, 2.3, 0.7, 0.01, 1) : AHDFPMParameters(0.6, 1.75, 0.0001, 0.001, 0.4, 2.4, 0.8, 0.1, 1)
end


#used for cs_proccessing
function getAlgorithmParams(::Type{CGDFPAlgorithm}; experiment_name::String="numerical")
    experiment_name == "numerical" ? CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3) : CGDFPParameters(0.6, 1, 0.001, 0.7, 0.3)
end

#used for cs_proccessing
function getAlgorithmParams(::Type{MOPCGAlgorithm}; experiment_name::String="numerical")
    experiment_name == "numerical" ? MOPCGParameters(0.1, 0.9, 0.0001, 0.1) : MOPCGParameters(1, 0.6, 0.1, 0.2)
end

function getAlgorithmParams(::Type{RILSDFPMAlgorithm}; experiment_name::String="numerical")
    experiment_name == "numerical" ? RILSDFPMParameters(0.24561134935215856, 0.8564877572100413, 0.20166231652151723, 0.3270752453777446, 0.8111659397644992, 0.6153130406009055, 0.7110051332386418, 1.7965) : RILSDFPMParameters()
end
