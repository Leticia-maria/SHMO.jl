using MolecularGraph
using LinearAlgebra
using Statistics
using DataStructures

function compound(SMILES::String) 
    mol = smilestomol(SMILES)
    mol_ed = [] # molecule_edges extracted from MolecularGraph
    i_vert = []
    j_vert = []
    for edge in mol.edges
        push!(mol_ed, edge, reverse(edge))
    end
    for ed in mol_ed
        push!(i_vert, ed[1])
        push!(j_vert, ed[2])
    end

    neighbour_list = tuple(i_vert, j_vert)
    return neighbour_list
end

function edges2neighbors(xs, ys)
    max_vertex = max(maximum(xs), maximum(ys))
    neighbours = Vector{Int}[Int[] for i in 1:max_vertex]
    for (i_v, j_v) in zip(xs, ys)
        push!(neighbours[i_v], j_v)
    end
    return neighbours
end

a = edges2neighbors(compound("C1=CC=CC=CC=C1")[1], compound("C1=CC=CC=CC=C1")[2])

print(a)

function inner_dfs_cycle_detection!(cycles, candidate, neighbours, cycle_length)
    @assert length(candidate) <= cycle_length
    for next_vertex in neighbours[last(candidate)]
        if next_vertex == first(candidate)
            if length(candidate) < cycle_length
                continue
            else
                push!(cycles, copy(candidate))
            end
        elseif next_vertex in candidate
            # Do not cycle back to the middle of candidate.
            continue
        elseif length(candidate) < 1
            # If the candidate does not have a subcycle and it is smaller
            # than the desired cycle length, then recurse.
            push!(candidate, next_vertex)
            inner_dfs_cycle_detection!(cycles, candidate, neighbours, cycle_length)
        end
    end
    return
end

function dfs_cycle_detection(neighbours, cycle_length)
    cycles = Vector{eltype(neighbours)}()
    for initial_vertex in eachindex(neighbours)
        inner_dfs_cycle_detection!(cycles, [initial_vertex], neighbours, cycle_length)
    end
    return cycles
end

function unique_cycles(cycles)
    (x -> first.(x)).(unique!(map(c -> Set(c .=> circshift(c, 1)), cycles)))
end

cycles_base_1 = dfs_cycle_detection(a, 8)
unique_cycles_base_1 = unique_cycles(cycles_base_1)
println(unique_cycles_base_1)


