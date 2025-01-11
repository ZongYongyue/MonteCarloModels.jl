function BTLattice(Lx, Ly = Lx)
    uc = UnitCell(
        "Bond-depended Triangular",
        (Float64[sqrt(3.0)/2, -0.5], Float64[sqrt(3.0)/2, +0.5]),
        [Float64[0.0, 0.0]],
        [
            Bond(1, 1, ( 0,  1), 1),
            Bond(1, 1, (-1,  0), 1),
            Bond(1, 1, (+1, -1), 1),

            Bond(1, 1, (+1,  0), 2),
            Bond(1, 1, (-1, +1), 2),
            Bond(1, 1, ( 0, -1), 2),
        ]
    )
    return Lattice(uc, (Lx, Ly))
end
