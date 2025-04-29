export J, calculateG, delta_omega

"""
fucntion to calculate spectral density
"""
function J(G::Vector, t::Vector, ω::Real)::Real
    return 2 * trapz(t, (G .* cos.(ω .* t)))
end

function calculateG(rij::Matrix)

    zvecs = [
        [0.0, 0.0, 1.0],
        #=[√(8/9), 0, -1/3],=#
        #=[-√(2/9), √(2/3), -1/3],=#
        #=[-√(2/9), -√(2/3), -1/3]=#
    ]

    nframes = size(rij,1)
    phi_ij = zeros(Float64, nframes)
    tmp = zeros(nframes)
    rijn = 0.0

    for p in axes(rij, 2) # for each pair

        for z in zvecs # for each z vector

            for f in 1:nframes 

                rijn = norm(rij[f,p])
                phi_ij[f] = (3 * (dot(rij[f,p],z)/rijn)^2 - 1) / rijn^3

            end

            tmp += ACF(phi_ij)
        end

    end

    G_ens_av = tmp / (size(rij,2) > 2000 ? 1024 : 512)

    return prefactor * G_ens_av * 1e60 
end

function delta_omega(F::Matrix)
    return 3 * prefactor * mean(mean(F.^2, dims=2)) * 1e60
end

