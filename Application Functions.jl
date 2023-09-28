include("BlackScholes.jl")


function apply_BS_model()
    
    if model == "1976"
        //

        return BlackScholes(param)
    elseif model == "Asay"
        return BlackScholes(param)
    elseif model == "Garman-Kohlhagen"
        return BlackScholes(param)
    else
        return "Error: Model not recognized"
    end
end
    

