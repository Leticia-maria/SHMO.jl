__precompile__()

module SHMO

import Pkg
Pkg.add("LinearAlgebra")
Pkg.add("MolecularGraph")

using MolecularGraph, LinearAlgebra

#--- GUI → file 1

#--- TODO: build AdjacencyMatrix and AdjacencyList → file 2

Molecule = smilestomol(add_widget(smiles))

"""The function '''adjmatrix''' gives an adjacency matrix of a molecule (considered an undirected graph)"""
function adjmatrix(smiles_code::String)
    Molecule = smilestomol(smiles_code)
    bonds = Molecule.edges
    append!(bonds, reverse.(bonds))
    Natoms = atomcount(Molecule)
    AdjMatrix = zeros(Int64, Natoms, Natoms)

    for i in 1:length(bonds)
        AdjMatrix[bonds[i][1], bonds[i][2]] = 1
    end
    return AdjMatrix
end

"""The function '''get_eigvals''' returns the eigenvalues of an adjacency matrix"""
function get_eigvals(adj_matrix::Matrix{Int64})
    vals = eigvals(adj_matrix)
    return vals
end


end # module
