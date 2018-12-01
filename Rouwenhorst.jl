"""
Rouwenhorst's method to approximate an AR(1) process that follows
    y_t = μ + ρ y_{t-1} + ε_t,
where ε_t ~ N (0, σ^2)
##### Arguments
- `N::Integer` : Number of points in Markov process
- `ρ::Real` : Persistence parameter in AR(1) process
- `σ::Real` : Standard deviation of random component of AR(1) process
- `μ::Real(0.0)` : Constant of AR(1) process
##### Returns
- `Λ::LinRange{Float64}` : nx1 vector of equally spaced state states centered around y_bar = μ/(1-ρ)
- `P::Array{Float64,2}` : Markov transition matrix
"""
function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1-ρ^2)
    p  = (1+ρ)/2
    Θ = [p 1-p; 1-p p]
    ψ = sqrt(N-1) * σ_y
    m = μ / (1 - ρ) # Not sure if I should center the space around m or μ
    Λ, P = _rouwenhorst(p, p, m, ψ, N)
    return Λ, P
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m-Δ, m+Δ],  [p 1-p; 1-q q]
    else
        _, θ_nm1 =     _rouwenhorst(p, q, m, Δ, n-1)
              θN =     p*[θ_nm1 zeros(n-1, 1); zeros(1, n)] +
                   (1-p)*[zeros(n-1, 1) θ_nm1; zeros(1, n)] +
                       q*[zeros(1, n); zeros(n-1, 1) θ_nm1] +
                   (1-q)*[zeros(1, n); θ_nm1 zeros(n-1, 1)]
        θN[2:end-1, :] ./= 2
        return LinRange(m-Δ, m+Δ, n), θN
    end
end
