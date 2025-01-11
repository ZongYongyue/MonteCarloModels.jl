@with_kw_noshow struct BTHubbardModel <: Model
    mu::Float64 = 0.0
    U::Float64 = 1.0
    @assert U >= 0. "U must be positive."
    tre::Float64 = 1.0
    tim::ComplexF64 = 0.0
    l::Lattice{2}
    boundary::Symbol = :periodic
    @assert l.unitcell.name == "Bond-depended Triangular"
end

# Constructors
BTHubbardModel(params::Dict{Symbol}) = BTHubbardModel(; params...)
BTHubbardModel(params::NamedTuple) = BTHubbardModel(; params...)
function BTHubbardModel(lattice::Lattice{2}; kwargs...)
    BTHubbardModel(l = lattice; kwargs...)
end
function BTHubbardModel(L, dims; kwargs...)
    l = choose_lattice(BTHubbardModel, dims, L)
    BTHubbardModel(l = l; kwargs...)
end
choose_lattice(::Type{<:BTHubbardModel}, dims, L) = BTLattice(L)
unique_flavors(::BTHubbardModel) = 2
total_flavors(::BTHubbardModel) = 2
lattice(m::BTHubbardModel) = m.l

choose_field(m::BTHubbardModel) = m.U < 0.0 ? MagneticHirschField : DensityHirschField

hopping_eltype(::BTHubbardModel) = ComplexF64
hopping_matrix_type(::AbstractField, ::BTHubbardModel) = BlockDiagonal{ComplexF64, 2, CMat64}
greens_eltype(::AbstractField, ::BTHubbardModel) = ComplexF64
greens_matrix_type( ::AbstractField, ::BTHubbardModel) = BlockDiagonal{ComplexF64, 2, CMat64}

# cosmetics
import Base.summary
import Base.show
Base.summary(model::BTHubbardModel) = "BTHubbard model"
function Base.show(io::IO, model::BTHubbardModel)
    print(io, "BTHubbard model, $(length(model.l)) sites")
end
Base.show(io::IO, m::MIME"text/plain", model::BTHubbardModel) = print(io, model)

@inline parameters(m::BTHubbardModel) = (N = length(m.l), t = tre + tim*im, U = -m.U, mu = m.mu, boundary=m.boundary)

function hopping_matrix(m::BTHubbardModel)
    N = length(m.l)
    tup = diagm(0 => fill(-ComplexF64(m.mu), N))
    tdown = diagm(0 => fill(-ComplexF64(m.mu), N))
    tp = m.tre + m.tim
    tm = m.tre - m.tim 
    if m.boundary == :periodic
        bds = bonds(m.l, Val(true))
    elseif m.boundary == :open
        bds = bonds_open(m.l, Val(true))
    else
        throw(ArgumentError("boundary condition should be :periodic or :open"))
    end
    for b in bds
        # NN paper direction
        if b.label == 1 
            tup[b.from, b.to] = tp
            tdown[b.from, b.to] = tm
        # NN reverse direction
        else b.label == 2
            tup[b.from, b.to] = tm
            tdown[b.from, b.to] = tp
        end
    end

    return BlockDiagonal(StructArray(tup), StructArray(tdown))
end
