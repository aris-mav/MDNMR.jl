export getpairs!, periodicboundary!
"""
function which takes an array of positions and modifies the Hpairs array to
store all the combinations of vectors connecting the hydrogens
(Hpairs is allocated within calculateF function)
"""
function getpairs!(Hpairs::Vector{SVector{3, Float32}},
                   positions::Vector{SVector{3, Float32}},
                   combinations::String)

    counter::Int64 = 1

    if combinations == "all"
        for k in axes(positions, 1)
            for j in axes(positions, 1)
                j > k || continue
                Hpairs[counter] = positions[k] - positions[j]
                counter += 1
            end
        end

    elseif combinations == "intra"
        for i in 1:2:length(positions)
            Hpairs[counter] = positions[i] - positions[i+1]
            counter += 1
        end

    elseif combinations == "inter"
        for i in axes(positions, 1)
            for j in axes(positions, 1)
                (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue
                Hpairs[counter] = positions[i] - positions[j]
                counter += 1
            end
        end

        ## Just a sanity check of how the loop above works, not actual part of the code
        #=for i in 1:10=#
        #=     for j in 1:10=#
        #=        (isodd(i) && j > (i + 1)) || (iseven(i) && j > i) || continue=#
        #=        display("$i with $j")=#
        #=     end=#
        #=end=#
    end

end

"""
applies periodic boundary to the coordinates of the atoms
according to the box size
"""
function periodicboundary!(xyz::Vector{SVector{3, Float32}}, box::MVector{3, Float32})
    displacement_vector::MVector{3, Float32} = @MVector zeros(Float32, 3)
    for i in 1:length(xyz)
        for l in 1:3
            if xyz[i][l] > box[l] / 2
                displacement_vector[l] -= box[l]
            elseif xyz[i][l] < -box[l] / 2
                displacement_vector[l] += box[l]
            end
        end
        xyz[i] += displacement_vector
        displacement_vector .= 0.0
    end
end


"""
converts [x, y, z] to [r, θ, φ]
"""
function cart2sph!(rtp::Vector{SVector{3, Float32}}, xyz::Vector{SVector{3, Float32}})
    for i in eachindex(rtp)
        rtp[i] = @SVector [ norm(xyz[i]), acos(xyz[i][3] / norm(xyz[i])), atan(xyz[i][2], xyz[i][1]) ]
    end
end

