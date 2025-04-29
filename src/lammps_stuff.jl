export time_array, collect_rij

"""
Collect all possible distances between all pairs (rows), for all timesteps (columns)
"""
function collect_rij(dumpfilepath, contributions::String)

    num_lines::Int = countlines(dumpfilepath)

    natoms::Int = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    totalsteps::Int = num_lines ÷ (natoms + 9) 

    #=exclude_atoms = natoms ÷ 2=#
    exclude_atoms = 0
    natoms -= exclude_atoms
    nhydrogens::Int = 2 * natoms ÷ 3

    npairs::Int64 = 0 

    if contributions == "all"
        npairs += nhydrogens * (nhydrogens - 1) ÷ 2
    elseif contributions == "inter"
        npairs += nhydrogens * (nhydrogens - 1) ÷ 2 - nhydrogens ÷ 2 
    elseif contributions == "intra"
        npairs += nhydrogens ÷ 2
    end

    @show npairs
    @show nhydrogens

    # Initialise arrays
    boxlengths::MVector{3, Float32} = @SVector zeros(Float32, 3)
    positions::Vector{SVector{3, Float32}} = [@SVector zeros(Float32, 3) for _ in 1:nhydrogens]
    Hpairs::Vector{SVector{3, Float32}} = [@SVector zeros(Float32, 3) for _ in 1:npairs]
    P::Matrix{SVector{3, Float32}} = fill(SA_F32[0,0,0], totalsteps, npairs)

    # Open the dump file
    open(dumpfilepath) do io

        # Loop over time steps
        for s in 1:totalsteps

            # Print progress (optional)
            if isinteractive()
                if s in floor.(Int, collect((totalsteps/10):(totalsteps/10):totalsteps))
                    progresspercent = ceil(s * 100 / totalsteps)
                    display("Calculation progress: $progresspercent %")
                end
            end

            # Skip headers
            for _ in 1:5
                readline(io)
            end

            # Read box bounds
            for i in 1:3
                boxlengths[i] = sum(abs.(parse.(Float32, split(readline(io), ' '))))
            end

            # Skip header
            readline(io)

            # Read positions
            for i in 1:2:nhydrogens
                readline(io) # skip the oxygen
                positions[i] = SVector( parse.(Float32, split(readline(io), ' '))[3:5]...)
                positions[i+1] = SVector( parse.(Float32, split(readline(io), ' '))[3:5]...)
            end

            readuntil(io,"TIME")

            # Do calculations
            getpairs!(Hpairs, positions, contributions)
            periodicboundary!(Hpairs, boxlengths)

            P[s, :] = Hpairs

            # Go to next timestep
        end

        # Close the IO file
    end

    return P

    # Exit the function
end


"""
extract a time array from the dump file
"""
function time_array(dumpfilepath, timestep)

    num_lines::Int32 = countlines(dumpfilepath)

    natoms = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    totalsteps = num_lines ÷ (natoms + 9) 

    dump_every_how_many_steps = open(dumpfilepath) do io
        readuntil(io, "TIMESTEP")
        readuntil(io, "TIMESTEP")
        readline(io)
        return  parse(Int,readline(io))
    end
    
    t = collect(Float32, 0:totalsteps-1) .* (dump_every_how_many_steps * timestep) .* 1e-15;

    return t
end

