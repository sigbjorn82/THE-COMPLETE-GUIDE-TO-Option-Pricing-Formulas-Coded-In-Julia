using Distributions

# Option Parameters Struct
mutable struct OptionParams
    option::String  # Type of option (Call or Put)
    S::Float64  # Stock price or currency pair spot price
    F::Float64  # Futures price
    X::Float64  # Strike price of an option
    T::Float64  # Time to maturity
    r::Float64  # Risk-free interest rate
    r_f::Float64  # Foreign risk-free interest rate
    σ::Float64  # Volatility of the underlying stock price
    q::Float64  # Dividend yield

    OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
        new(option, S, F, X, T, r, r_f, σ, q)
end



function BlackScholes(param::OptionParams)
    
    d1 = (log(param.S/param.X) + (param.r + 0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*exp(-param.r*param.T)*cdf(Normal(0,1),-d2) - param.S*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
    """
    This is the standard Black-Schoules formula for pricing 
    European options without Merton extention.

    input:
    param:: OptionParams [mutable struct with parameters of option applicable to function]

    OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
        new(option, S, F, X, T, r, r_f, σ, q)
        
    example: 
            option_BS = OptionParams("Call", S=60.00, X=65.00, T=0.25, r=0.08, σ=0.3)
            BlackScholes(option_BS) = 2.13337...
    """
end

function BlackScholesMerton(param::OptionParams)
    
    d1 = (log(param.S/param.X) + (param.r - param.q + 0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S * exp(-param.q * param.T)*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X * exp(-param.r*param.T)*cdf(Normal(0,1),-d2) - param.S * exp(-param.q * param.T)*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
    """
    This is Black-Schoules-Merton formula for pricing 
    European options with Merton extention for dividende yield.

        input:
        param:: OptionParams [mutable struct with parameters of option applicable to function]
    
        OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
            new(option, S, F, X, T, r, r_f, σ, q)
            
    example: 
            option_BSM = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)
            BlackScholesMerton(option_BSM) = 2.46478...
    """
end



function BlackScholesFutures(param::OptionParams)
    
    d1 = (log(param.F/param.X) + 0.5*param.σ^2*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return exp(-param.r*param.T)*(param.F*cdf(Normal(0,1),d1) - param.X*cdf(Normal(0,1),d2))
    elseif param.option == "Put"
        return exp(-param.r*param.T)*(param.X*cdf(Normal(0,1),-d2) - param.F*cdf(Normal(0,1),-d1))
    else
        return "Error: Option type not recognized"
    end
    """
    This is the Black-Schoules (1976) formula for pricing options, 
    when underlying asset is Futures or Fwd contract.

    input:

    input:
    param:: OptionParams [mutable struct with parameters of option applicable to function]

    OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
        new(option, S, F, X, T, r, r_f, σ, q)
        
    example:
            option_BSF = OptionParams("Call", F=19.00, X=19.00, T=0.75, r=0.1, σ=0.28)
            BlackScholesFutures(option_BSF) = 1.70105...
    """
end 


function BlackScholesAsay(param::OptionParams)
    
    d1 = (log(param.F/param.X) + (0.5*param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.F*cdf(Normal(0,1),d1) - param.X*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*cdf(Normal(0,1),-d2) - param.F*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
    """
    This is the standard Black-Schoules Asay modificated formula for pricing 
    Futures options. Extention is same as Black-76 formula without intrest 
    rate term.

    input:
    param:: OptionParams [mutable struct with parameters of option applicable to function]

    OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
        new(option, S, F, X, T, r, r_f, σ, q)
        
    example: 
            option_BSA = OptionParams("Put", F=4200.00, X=3800.00, T=0.75, σ=0.15)
            BlackScholesAsay(option_BSA) = 65.61854...
    """
end

function BlackScholesCurrency(param::OptionParams)
    
    d1 = (log(param.S/param.X) + (param.r - param.r_f + 0.5* param.σ^2)*param.T)/(param.σ*sqrt(param.T))
    d2 = d1 - param.σ*sqrt(param.T)

    if param.option == "Call"
        return param.S*exp(-param.r_f * param.T)*cdf(Normal(0,1),d1) - param.X*exp(-param.r * param.T)*cdf(Normal(0,1),d2)
    elseif param.option == "Put"
        return param.X*exp(-param.r * param.T)*cdf(Normal(0,1),-d2) - param.S*exp(-param.r_f*param.T)*cdf(Normal(0,1),-d1)
    else
        return "Error: Option type not recognized"
    end
    """
    This is Garman and Kohlhagen (1983) modified Black-Scholes model can
    be used to price European currency options. equivalent to Black-Scholes-Merton
    with foreign risk-free interest rate replacing the dividend yield.

    NB! since for currency pair (i.e USD/EUR) options a call on the quote is equivalent 
    to a put on the base currency (USD-call/EUR-put) and vice versa, the formula input 
    refers to the base currency (from example above "call" gives USD-Call)


    input:
    param:: OptionParams [mutable struct with parameters of option applicable to function]

    OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
        new(option, S, F, X, T, r, r_f, σ, q)
        
    example: 
                option_BSC = OptionParams("Call", S=1.56, X=1.60, T=0.5, r=0.06, r_f=0.08, σ=0.12)
                BlackScholesCurrency(option_BSC)= 0.02909...
    """
end



function GBlackScholes(param::OptionParams)
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


    """
    This is Black-Schoules-Merton formula for pricing 
    European options with Merton extention for dividende yield.

        input:
        param:: OptionParams [mutable struct with parameters of option applicable to function]
    
        OptionParams(option::String; S::Float64=0.0, F::Float64=0.0, X::Float64=0.0, T::Float64=0.0, r::Float64=0.0, r_f::Float64=0.0, σ::Float64=0.0, q::Float64=0.0) = 
            new(option, S, F, X, T, r, r_f, σ, q)
           

    b = r:               gives the Black and Scholes (1973) stock 
                         option model.

    b = r — q:           gives the Merton (1973) stock option model 
                         with continuous dividend yield q.

    b = 0:               gives the Black (1976) futures option model.

    b = 0 and r = 0:     gives the Asay (1982) margined futures option
                         model.

    b= r — r_f:          gives the Garman and Kohlhagen (1983) currency
                         option model.

    example: 
            option_GBS = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)
            GBlackScholes(option_GBS) = GBlackScholes(option_GBS)= 2.46478...
    """
end

# Rekkefølge (option, S, F, X, T, r, r_f, σ, q)

option_BS = OptionParams("Call", S=60.00, X=65.00, T=0.25, r=0.08, σ=0.3)
option_BSM = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)
option_BSF = OptionParams("Call", F=19.00, X=19.00, T=0.75, r=0.1, σ=0.28)
option_BSA = OptionParams("Put", F=4200.00, X=3800.00, T=0.75, σ=0.15)
option_BSC = OptionParams("Call", S=1.56, X=1.60, T=0.5, r=0.06, r_f=0.08, σ=0.12)
option_GBS = OptionParams("Put", S=100.00, X=95.00, T=0.5, r=0.1, σ=0.2, q=0.05)


println("BS returns ", BlackScholes(option_BS))
println("BSM returns ", BlackScholesMerton(option_BSM))
println("BSF returns ", BlackScholesFutures(option_BSF))
println("BSA returns ", BlackScholesAsay(option_BSA))
println("BSC returns ", BlackScholesCurrency(option_BSC))
println("GBS returns ", GBlackScholes(option_GBS))



#1.2 PARITIES AND SYMMETRIES

