using Distributions

export BlackScholes, BlackScholesTypes, BS1976, BSMertonExt, BSFutures, BSAsay, BSCurrency, BSGeneralized

"""
`BlackScholesTypes`

Abstract type representing the generic concept of a Black-Scholes model for option pricing.
"""
abstract type BlackScholesTypes end

"""
`BS1976`

Mutable struct representing the Black-Scholes model as proposed in 1976. This model is used for the pricing of options.
The struct holds the following fields:

* `option`: the type of the option, can be "Call" or "Put"
* `S`: the current price of the underlying asset
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `r`: the risk-free interest rate
* `σ`: the volatility of the underlying asset
"""
mutable struct BS1976 <: BlackScholesTypes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
end

"""
`BSMertonExt`

Mutable struct representing the Merton extension to the Black-Scholes model, which includes a dividend yield.
The struct holds the same fields as `BS1976`, and an additional `q` field for the dividend yield.

* `option`: the type of the option, can be "Call" or "Put"
* `S`: the current price of the underlying asset
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `r`: the risk-free interest rate
* `σ`: the volatility of the underlying asset
* `q`: the dividend yield
"""
mutable struct BSMertonExt <: BlackScholesTypes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    q::Float64
end


"""
`BSFutures`

Mutable struct representing the Black-Scholes model for futures options.
The struct holds the following fields:

* `option`: the type of the option, can be "Call" or "Put"
* `F`: the current price of the futures contract
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `r`: the risk-free interest rate
* `σ`: the volatility of the underlying futures contract
"""
mutable struct BSFutures <: BlackScholesTypes
    option::String
    F::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
end


"""
`BSAsay`

Mutable struct representing the Asay model, a variant of the Black-Scholes model used for futures options.
The struct holds the following fields:

* `option`: the type of the option, can be "Call" or "Put"
* `F`: the current price of the futures contract
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `σ`: the volatility of the underlying futures contract
"""
mutable struct BSAsay <: BlackScholesTypes
    option::String
    F::Float64
    X::Float64
    T::Float64
    σ::Float64
end

"""
`BSCurrency`

Mutable struct representing the Garman-Kohlhagen model, a variant of the Black-Scholes model used for currency options.
The struct holds the following fields:

* `option`: the type of the option, can be "Call" or "Put"
* `S`: the current exchange rate of the foreign currency
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `r`: the domestic risk-free interest rate
* `σ`: the volatility of the exchange rate
* `r_f`: the foreign risk-free interest rate
"""
mutable struct BSCurrency <: BlackScholesTypes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    r_f::Float64
end

"""
`BSGeneralized`

Mutable struct representing a generalized Black-Scholes model that allows for a cost-of-carry rate `b`.
The struct holds the following fields:

* `option`: the type of the option, can be "Call" or "Put"
* `S`: the current price of the underlying asset
* `X`: the strike price of the option
* `T`: the time to maturity of the option
* `r`: the risk-free interest rate
* `σ`: the volatility of the underlying asset
* `q`: the dividend yield
"""
mutable struct BSGeneralized <: BlackScholesTypes
    option::String
    S::Float64
    X::Float64
    T::Float64
    r::Float64
    σ::Float64
    q::Float64
end







"""
    BlackScholes(param::BS1976)

Computes the Black-Scholes price of a call or put option with parameters provided by the `param` struct.
The model used is the original Black-Scholes model from 1976. 
Throws an error if the option type is not recognized.
"""
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

"""
    BlackScholes(param::BSMertonExt)

Computes the Black-Scholes price of a call or put option with parameters provided by the `param` struct.
The model used is the Merton extension to the Black-Scholes model, which includes a dividend yield.
Throws an error if the option type is not recognized.
"""
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


"""
    BlackScholes(param::BSFutures)

Computes the Black-Scholes price of a call or put futures option with parameters provided by the `param` struct.
The model used is the Black-Scholes model for futures options. 
Throws an error if the option type is not recognized.
"""
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

"""
    BlackScholes(param::BSAsay)

Computes the price of a call or put futures option with parameters provided by the `param` struct using the Asay model.
Throws an error if the option type is not recognized.
"""
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

"""
    BlackScholes(param::BSCurrency)

Computes the price of a call or put currency option 
with parameters provided by the `param` struct using 
the Garman-Kohlhagen model. Throws an error if the option 
type is not recognized.
"""
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


"""
    BlackScholes(param::BSGeneralized)

Computes the price of a call or put option with parameters provided 
by the `param` struct using a generalized Black-Scholes model that 
allows for a cost-of-carry rate `b`. Throws an error if the option 
type is not recognized.
"""
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

#Testing functions:

    #test1  = BS1976("Call", 60.0, 65.0, 0.25, 0.08, 0.3)
    #test2 = BSMertonExt("Put", 100.0, 95.0, 0.5, 0.1, 0.2, 0.05)
    #test3 = BSFutures("Call", 19.0, 19.0, 0.75, 0.1, 0.28)
    #test4 = BSAsay("Put", 4200.00, 3800.00, 0.75, 0.15)
    #test5 = BSCurrency("Call", 1.56, 1.6, 0.5, 0.06, 0.08, 0.12)
    #test6 = BSGeneralized("Put", 100.00, 95.00, 0.5, 0.1, 0.2, 0.05)

    #BS = BlackScholes(test1)
    #BSM = BlackScholes(test2)
    #BSF = BlackScholes(test3)
    #BSA = BlackScholes(test4)
    #BSC = BlackScholes(test5)
    #BSG = BlackScholes(test6)

   #println(BS, "\n", BSM, "\n", BSF, "\n", BSA, "\n", BSC, "\n", BSG)

# Path: THE COMPLETE GUIDE TO Option Pricing Formulas Coded/BlackScholes.jl