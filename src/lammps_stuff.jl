export time_array, calculateF
"""
F in this case is the quantity (3cos(θ)^2 - 1) / r^3
The output of this function is a matrix which contains this quantity 
for all pairs of hydrogens (columns)
and each time step (rows).
"""
function calculateF(dumpfilepath, contributions::String)

    num_lines::Int = countlines(dumpfilepath)

    natoms::Int = open(dumpfilepath) do io
        for _ in 1:3
            readline(io)
        end
        return parse(Int, readline(io))
    end

    totalsteps::Int = num_lines ÷ (natoms + 9) 

    exclude_atoms = natoms ÷ 2
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
    F::Matrix{Float32} = zeros(Float32, npairs, totalsteps)
    zvec::SVector{3, Float32} = SA_F32[0.0, 0.0, 1.0]
    vecnorm::Float32 = 1.0

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

            for (i, p) in enumerate(Hpairs)
                vecnorm = norm(p)
                F[i, s] = ( 3 * ( dot(zvec, p)/vecnorm)^2 -1 ) / vecnorm ^3
            end


            # Go to next timestep
        end

        # Close the IO file
    end

    return F

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

