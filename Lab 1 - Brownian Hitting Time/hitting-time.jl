using SpecialFunctions: gamma
using Random: randn, AbstractRNG, MersenneTwister
using GLM, StatsBase, DataFrames, Plots
plotlyjs()   # set the backend to Plotly

"""
Compute the true value of E[τ ∧ T]
"""
function trueE(T, a)
    x = a * a / (2. * T)
    return T * (1 - 1 / sqrt(pi) * gamma(0.5, x)) + a * a / (2 * sqrt(pi)) * gamma(-0.5, x)
end

"""
Sample the Brownian motion on the kT/N for k=1,...,N. The value at time 0 is not included
The Brownian path is stored in `W` on exit
"""
function sampleW!(W::Vector{<:Real}, T::Real, N::Integer, rng::AbstractRNG)
    @assert length(W) == N
    sqrt_dt = sqrt(T/N)
    prev = 0.
    for i in 1:length(W)
        W[i] = prev + sqrt_dt * randn(rng)
        prev = W[i]
    end
end

"""
Naive approximation of τ
"""
function computeTau(W::Vector{<:Real}, a, nDates, T)
    tau = findfirst(w -> w > a, W)
    if isnothing(tau)
        return T
    else
        return tau * T / nDates
    end
end


"""
Refine the computation of τ when the Brownian motion exceeds `a` - `threshold` using `nSubTicks` more intermediate dates.
"""
function computeRefinedTau(W::Vector{<:Real}, a::Real, nDates::Integer, T::Real, nSubTicks::Integer, threshold::Real, rng::AbstractRNG)
    dt = T / nDates
    sub_dt = T / (nSubTicks * nDates)
    t = dt
    for i in 1:length(W)
        if abs(W[i] - a) < threshold
            Wtj = W[i]
            for j in 1:nSubTicks
                Wtj += randn(rng) * sqrt(sub_dt)
                if Wtj >= a
                    return t + j * sub_dt
                end
            end
        end
        t += dt
    end
    return computeTau(W, a, nDates, T)
end


function mcTau(a, T, nDates, mSamples, rng::AbstractRNG, tauFunction::Function)
    W = Vector{Float64}(undef, nDates)
    average = 0.
    squares = 0.
    for m in 1:mSamples
            sampleW!(W, T, nDates, rng)
            tau = tauFunction(W, a, nDates, T)
            average += tau
            squares += tau * tau
    end
    average /= mSamples
    stdDev = sqrt((squares / mSamples - average * average) / mSamples)
    return average
end

function main()
    T = 3.
    a = 1.
    M = 1E7
    rng = MersenneTwister()
    true_value = trueE(T, a)
    println("True value: ", true_value)
    N = [10, 20, 50, 100, 125, 250, 500]
    logN = log.(N)

    # Regression for the naive approach
    println("-- Méthode naive")
    taus = Vector{Float64}(undef, length(N))
    taus .= mcTau.(a, T, N, M, rng, computeTau)
    logError = log.(taus .- true_value)
    data = DataFrame(X=logN, Y=logError) # .- true_value)
    tauRegression = lm(@formula(Y ~ X), data)
    # println(tauRegression)
    plot(logN, GLM.coef(tauRegression)[1] .+ logN .* GLM.coef(tauRegression)[2], label="naive regression", legend = :topright)
    display(scatter!(log.(N), logError, label="naive"))
    println("R2 pour la méthode naive ", r2(tauRegression))
    println("Ordre de la méthode naive = $(GLM.coef(tauRegression)[2])")


    # Regression for the refined approach
    println("-- Méthode par ponts browniens")
    nSubTicks = 4
    threshold = 0.1
    refinedTauFunction(W, a, nDates, T) = computeRefinedTau(W, a, nDates, T, nSubTicks, threshold, rng)
    taus2 = Vector{Float64}(undef, length(N))
        taus2 .= mcTau.(a, T, N, M, rng, refinedTauFunction)
    logError = log.(abs.(taus2 .- true_value))
    data = DataFrame(X=logN, Y=logError) # .- true_value)
    tauRegression = lm(@formula(Y ~ X), data)
    # println(tauRegression)
    plot!(logN, GLM.coef(tauRegression)[1] .+ logN .* GLM.coef(tauRegression)[2], label="bridge regression", legend = :topright)
    display(scatter!(log.(N), logError, label="bridge error"))
    println("R2 pour la méthode par ponts Browniens ", r2(tauRegression))
    println("Ordre de la méthode par ponts Browniens = $(GLM.coef(tauRegression)[2])")
end

main()
