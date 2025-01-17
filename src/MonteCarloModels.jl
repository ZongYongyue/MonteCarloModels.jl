module MonteCarloModels
using LinearAlgebra: diagm
using MonteCarlo: Model, UnitCell, Bond, Lattice, bonds, bonds_open
using MonteCarlo: StructArray, CMat64, CVec64, BlockDiagonal, matrix_type
using MonteCarlo: Model, DQMC, AbstractField, MagneticHirschField, DensityHirschField
using MonteCarlo: @with_kw_noshow
import MonteCarlo: choose_lattice, lattice, hopping_eltype, hopping_matrix_type, hopping_matrix, greens_eltype, greens_matrix_type, choose_field, unique_flavors, total_flavors
import Base.summary, Base.show


include("lattices.jl")
export BTLattice

include("models/BTHubbardModel.jl")
export BTHubbardModel

include("correlator.jl")
export current_current_correlator


end#module