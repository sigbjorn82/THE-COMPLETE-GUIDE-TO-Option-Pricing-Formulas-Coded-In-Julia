using Distributions, BenchmarkTools

abstract type BlackScholes end


mutable struct BS1976 <: BlackScholes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
end


mutable struct BSMertonExt <: BlackScholes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    q::Float64
end


mutable struct BSFutures <: BlackScholes
    option::String
    F::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
end



mutable struct BSAsay <: BlackScholes
    option::String
    F::Float64
    X::Float64
    T::Float64
    σ::Float64
end


mutable struct BSCurrency <: BlackScholes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    r_f::Float64
end

mutable struct BSGeneralized <: BlackScholes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    q::Float64
end








function BlackScholes(param::BS1976)
    
    d1 = (log(param.S/param.X) + (param.r + 0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*exp(-param.r*param.T)*cdf(Normal(0,1),-d2) - param.S*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
end

function BlackScholes(param::BSMertonExt)
    
    d1 = (log(param.S/param.X) + (param.r - param.q + 0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S * exp(-param.q * param.T)*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X * exp(-param.r*param.T)*cdf(Normal(0,1),-d2) - param.S * exp(-param.q * param.T)*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
end



function BlackScholes(param::BSFutures)
    
    d1 = (log(param.F/param.X) + 0.5*param.σ^2*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return exp(-param.r*param.T)*(param.F*cdf(Normal(0,1),d1) - param.X*cdf(Normal(0,1),d2))
    elseif param.option == "Put"
        return exp(-param.r*param.T)*(param.X*cdf(Normal(0,1),-d2) - param.F*cdf(Normal(0,1),-d1))
    else
        return "Error: Option type not recognized"
    end
end 

function BlackScholes(param::BSAsay)
    
    d1 = (log(param.F/param.X) + (0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.F*cdf(Normal(0,1),d1) - param.X*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*cdf(Normal(0,1),-d2) - param.F*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
end

function BlackScholes(param::BSCurrency)
    
    d1 = (log(param.S/param.X) + (param.r - param.r_f + 0.5* param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S*exp(-param.r_f * param.T)*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*exp(-param.r * param.T)*cdf(Normal(0,1),-d2) - param.S*exp(-param.r_f*param.T)*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
end



function BlackScholes(param::BSGeneralized)
    b = param.r - param.q

    d1 = (log(param.S/param.X) + (b + 0.5*param.σ^2) * param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S*exp((b - param.r)*param.T)*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*exp(-param.r * param.T)*cdf(Normal(0,1),-d2) - param.S*exp((b - param.r) * param.T)*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
end
    
#option_BS = OptionParams("Call", S=60.00, X=65.00, T=0.25, r=0.08, σ=0.3)
#option_BSM = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)
#option_BSF = OptionParams("Call", F=19.00, X=19.00, T=0.75, r=0.1, σ=0.28)
#option_BSA = OptionParams("Put", F=4200.00, X=3800.00, T=0.75, σ=0.15)
#option_BSC = OptionParams("Call", S=1.56, X=1.60, T=0.5, r=0.06, r_f=0.08, σ=0.12)
#option_GBS = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)

#Currency (BlackScholes, option::String, S::Float64, X::Float64, T::Float64, r::Float64, σ::Float64, r_f::Float64)

#generalized (option::String, S::Float64, X::Float64, T::Float64, r::Float64, σ::Float64, q::Float64)


    test1  = BS1976("Call", 60.0, 65.0, 0.25, 0.08, 0.3)
    test2 = BSMertonExt("Put", 100.0, 95.0, 0.5, 0.1, 0.2, 0.05)
    test3 = BSFutures("Call", 19.0, 19.0, 0.75, 0.1, 0.28)
    test4 = BSAsay("Put", 4200.00, 3800.00, 0.75, 0.15)
    test5 = BSCurrency("Call", 1.56, 1.6, 0.5, 0.06, 0.08, 0.12)
    test6 = BSGeneralized("Put", 100.00, 95.00, 0.5, 0.1, 0.2, 0.05)

    BS = BlackScholes(test1)
    BSM = BlackScholes(test2)
    BSF = BlackScholes(test3)
    BSA = BlackScholes(test4)
    BSC = BlackScholes(test5)
    BSG = BlackScholes(test6)

   println(BS, "\n", BSM, "\n", BSF, "\n", BSA, "\n", BSC, "\n", BSG)

# Path: THE COMPLETE GUIDE TO Option Pricing Formulas Coded/BlackScholes.jl