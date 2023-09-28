using Plots, Statistics,Random, Distributions, MarketData,DataFrames, TimeSeries, SpecialFunctions

#Data
SPX = yahoo("^GSPC");
SPX[end]["Close"];
Ma = moving(mean, SPX["Close"], 50);

plot(from(SPX["Close"],Date(2022, 01, 01)));

df = DataFrame(from(yahoo("^GSPC")["Close"], Date(2022, 01, 01)));


#parameters
begin
	μ = mean(df[end-2:end,2]);
	σ = std(df[!,2]);
	S₀ = df[end,2];
end


#Montecarlo Simulation Algorithm 1

function GBM_sim(S0::Float64, r::Float64, mu::Float64, sigma::Float64, T::Float64, n::Int64)
    dt = T/n
    W = cumsum([sqrt(dt)*randn() for i in 1:n])
    St = S0* exp(-r * T) * exp.(-r)*exp.((mu-0.5*sigma^2)*dt .+ sigma*W)
    return St
end


GBM = GBM_sim(S₀, 0.02, 1.08, 0.3, 1.0, 100);
plot(GBM);


#Montecarlo Simulation Algorithm 1

function brownian_motion_geometric(S, r, sigma, T, N)
    dt = T / N
    sdt = sqrt(dt)
    log_S = log(S)
    W = randn(N) * sdt
    for i = 2:N
        log_S += (r - 0.5 * sigma^2) * dt + sigma * W[i]
        W[i] += W[i - 1]
    end
    S_out = exp(log_S)
    return S_out
end

X = zeros(100000)
	
for i in 1:100000

X[i] = brownian_motion_geometric(S₀, 0.02, 0.3, 1, 100)
end

Φₓ = zeros(1000)
			
for x in 1:1000
    Φₓ[x] = Φ(X[x], 4000)
end

Φₘ = mean(Φₓ)



Y = zeros(100000, 100)
        
Y[1,:] = GBM_sim(S₀, 0.02, 1.08, 0.3, 1.0, 100)
    
for i in 2:100000
    Y[i, :] = GBM_sim(S₀, 0.02, 1.08, 0.3, 1.0, 100)
end
    Φₘ₂ = mean([max(i, 0) for i in Y[:, end].-4000]);
end


mean([max(i, 0) for i in Y[:, end].-4000])



#Black Scholes Method

function black_scholes_call(S::Float64, K::Float64, r::Float64, sigma::Float64, T::Float64)
    d1 = (log(S/K) + (r + 0.5*sigma^2)*T)/(sigma*sqrt(T))
    d2 = d1 - sigma*sqrt(T)
    N(d1)*S - N(d2)*K*exp(-r*T)
end

function N(x::Float64)
    (1.0 + erf(x/sqrt(2.0))) / 2.0
end

Xᵦ = black_scholes_call(S₀, 4000.00, 0.02, 0.3, 1.0)

#Vector of all the tree estimates
Φᵥ = [Xᵦ, Φₘ, Φₘ₂]

#Standard deviation of the result between the three estimates
std(Φᵥ[1:2])


mean(Φᵥ)