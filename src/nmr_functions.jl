export J, calculateG, delta_omega

"""
fucntion to calculate spectral density
"""
function J(G::Vector, t::Vector, ω::Real)::Real
    return 2 * trapz(t, (G .* cos.(ω .* t)))
end

function calculateG(F::Matrix)
    G_ens_av = mean(ACF.(eachrow(F))) # Ensemble average (no prefactors)
    return prefactor * G_ens_av * 1e60 
end

function delta_omega(F::Matrix)
    return 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
end

