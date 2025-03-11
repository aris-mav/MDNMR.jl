module MDNMR

using DelimitedFiles
using LinearAlgebra
using FFTW
using Statistics
using Trapz
using StaticArrays

# Define physical constants
const γ = 267.52218744e6; # (rad s^-1 T^-1)
const ħ = 1.054571817e-34;  # (J s)
const μ₀ = 1.25663706212e-6; #N A^-2
const prefactor = (3 / 16) * (μ₀/(4π))^2 * ħ^2 * γ^4
export γ, ħ, μ₀, prefactor

include("./rdf.jl")
include("./acf.jl")
include("./nmr_functions.jl")
include("./coordinates.jl")
include("./lammps_stuff.jl")

end # module MDNMR
