using Downloads: download
using JLD2: load_object, save_object
using LinearAlgebra
using LinearSolve: LinearProblem, KrylovJL_GMRES, RFLUFactorization, solve
using OrdinaryDiffEq
using RecursiveArrayTools: VectorOfArray
using SparseArrays: SparseMatrixCSC, spzeros
using SPICE
using StaticArrays

# Plotting package
using GLMakie

# Include utilities
include(joinpath(@__DIR__, "utils", "spice.jl"))
include(joinpath(@__DIR__, "utils", "math.jl"))
include(joinpath(@__DIR__, "utils", "spline.jl"))
include(joinpath(@__DIR__, "utils", "gravity.jl"))

# Include astrodynamics functions
include(joinpath(@__DIR__, "state_representations.jl"))
include(joinpath(@__DIR__, "ephemerides.jl"))
include(joinpath(@__DIR__, "force_models.jl"))
include(joinpath(@__DIR__, "equations_of_motion.jl"))

# Include Q-Law specifics
#include(joinpath(@__DIR__, "QLawParameters.jl"))