"""
Radial distribution function
(to see if the model makes sense)


(NEEDS TO BE UPDATED FOR SVECTORS)


"""
function calculate_rdf(dumpfilepath)

    n_lines = countlines(dumpfilepath)

    n_atoms, r_max = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        n_atoms = parse(Int, readline(io))
        readline(io)
        r_max = (parse.(Float32, split(readline(io), ' ')))[2] * 0.7 # nk = less than half box length
        return n_atoms, r_max 
    end

    n_oxygen = n_atoms ÷ 3 
    totalsteps = n_lines ÷ (n_atoms + 9)
    n_pairs = n_oxygen * (n_oxygen - 1) ÷ 2

    # Initialise arrays
    boxlengths::Vector{Float32} = zeros(3)
    positions::Vector{Vector{Float32}} = [zeros(3) for _ in 1:n_oxygen]
    o_xyz::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:n_pairs]
    o_rtp::Vector{Vector{Float32}} = [[1.0, 1.0, 1.0] for _ in 1:n_pairs]
    h::Vector{Int} = zeros(Int, 1000)

    nk = length(h)
    dr = r_max / nk
    r = [dr * i for i in 0:nk-1]

    # Open the dump file
    open(dumpfilepath) do io
        
        # Loop over time steps
        for s in 1:totalsteps

            # Skip headers
            for _ in 1:5
                readline(io)
            end

            # Read box bounds
            for i in 1:3
                boxlengths[i] .= sum(abs.(parse.(Float32, split(readline(io), ' '))))
            end

            # Skip header
            readline(io)

            # Read positions
            for i in 1:n_oxygen
                positions[i] = parse.(Float32, split(readline(io), ' '))[3:5]
                readline(io) # skip the hydrogen
                readline(io) # skip the hydrogen
            end

            # Do calculations
            getpairs!(o_xyz, positions, "all")
            periodicboundary!(o_xyz, boxlengths)
            cart2sph!(o_rtp, o_xyz)

            for rij in getindex.(o_rtp, 1)

                k = floor(Int, rij/dr) + 1
                if k <= nk
                    h[k] += 2
                end

            end

            # Go to next timestep
        end

        # Close the IO file
    end
    
    ri = r # lower shells
    rj = push!(r[2:end], r[end]+dr) # upper shells
    g = h ./(901 * 501 * (4π /3) * (6.022/180) .* (rj.^3 .- ri.^3)) 

    return r, g  
    # Exit the function
end

