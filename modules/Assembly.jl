# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

using Memoize
using SparseArrays

"""
    initialize_assembly!(mesh::Mesh)

Initialize the assembly process for the given mesh by computing the necessary geometric quantities.

# Arguments
- `mesh::Mesh`: The mesh object for which the assembly is initialized.
"""
function initialize_assembly!(mesh::Mesh)
    get_Bk!(mesh)
    get_detBk!(mesh)
    get_invBk!(mesh)
end

########################### GLOBAL ASSEMBLER ########################### 
"""
    assemble_global(mesh::Mesh, local_assembler!)

Assemble the global stiffness matrix and force vector for the given mesh using the provided local assembler function.

# Arguments
- `mesh::Mesh`: The mesh object.
- `local_assembler!`: A function that assembles the local stiffness matrix and force vector.

# Returns
- `K::SparseMatrixCSC`: The global stiffness matrix.
- `f::Vector`: The global force vector.
"""
function assemble_global(mesh::Mesh, local_assembler!)
    T = mesh.T
    p = mesh.p
    n_points = size(p,2)
    n_triangles = size(T, 2)
    rows = []
    cols = []
    data = Float64[]
    
    rows_f = []
    data_f = Float64[]

    A_loc = zeros(3,3)
    f_loc = zeros(3)
    for k in 1:n_triangles
        local_assembler!(A_loc, f_loc, mesh, k) 
        indices = T[:, k]
    
        for i in 1:3
            i_glob = indices[i]
            for j in 1:3
                j_glob = indices[j]
                append!(rows, i_glob)
                append!(cols, j_glob)
                append!(data, A_loc[i,j])
            end
            append!(rows_f, i_glob)
            append!(data_f, f_loc[i])
        end    
    end
    A_glob = sparse(rows, cols, data, n_points, n_points)
    F_glob = sparse(rows_f, ones(size(rows_f)), data_f)
    return A_glob, F_glob
end

"""
    impose_dirichlet(A, b, g, mesh)

Impose Dirichlet boundary conditions on the system.

# Arguments
- `A`: The global stiffness matrix.
- `b`: The global force vector.
- `g`: The Dirichlet boundary condition function.
- `mesh::Mesh`: The mesh object.

# Returns
- `A_cond`: The modified stiffness matrix with Dirichlet conditions imposed.
- `b_cond`: The modified force vector with Dirichlet conditions imposed.
- `uh`: The solution vector with Dirichlet conditions applied.
"""
function impose_dirichlet(A, b, g, mesh)
    ######################
    ### COMPLETARE QUI ###
    ######################
end

########################################################################
########################### LOCAL ASSEMBLERS ###########################
########################################################################

########################### POISSON PROBLEM ###########################
"""
    shapef_2DLFE(quadrule::TriQuad)

Compute the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `shapef`: The shape functions evaluated at the quadrature points.
"""
@memoize function shapef_2DLFE(quadrule::TriQuad)
    p = quadrule.points
    phi_1 = (x) -> - x[1] -x[2] +1
    phi_2 = (x) -> x[1]
    phi_3 = (x) -> x[2] 
    if size(p,1) == 2
        # println(size(p))
        values = zeros(3, size(p,2))
        for i in 1:size(p,2)
            # println(i)
            values[:, i] = [phi_1(p[:, i]), phi_2(p[:, i]), phi_3(p[: , i])] 
        end
        return values
    end
end

"""
    ∇shapef_2DLFE(quadrule::TriQuad)

Compute the gradients of the shape functions for the Poisson problem.

# Arguments
- `quadrule::TriQuad`: The quadrature rule.

# Returns
- `∇shapef`: The gradients of the shape functions evaluated at the quadrature points.
"""
@memoize function ∇shapef_2DLFE(quadrule::TriQuad)
    ∇_phi_1 = (x) ->  [-1; -1]
    ∇_phi_2 = (x) ->  [1; 0]
    ∇_phi_3 = (x) ->  [0; 1]
    p = quadrule.points
    if size(p,1) == 2
        gradients = zeros(2,3,size(p,2))
        # println([∇_phi_1(p[:, 1])  ∇_phi_2(p[:, 1])  ∇_phi_3(p[:, 1])])
        # println(gradients)
        for i in 1:size(p, 2)
            gradients[:, :, i] =  [∇_phi_1(p[:, i])  ∇_phi_2(p[:, i]) ∇_phi_3(p[:, i])]
        end
        return gradients
    end
end

"""
    poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)

Assemble the local stiffness matrix and force vector for the Poisson problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
function poisson_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f)
    B, a = get_Bk!(mesh)
    detB = get_detBk!(mesh)
    invB= get_invBk!(mesh)
    Bk = B[:, :, cell_index]
    ak = a[:, cell_index]
    detBk = detB[cell_index]
    invBk = invB[:, :, cell_index]


    quadr_matrix = Q2_ref
    phi_grad = ∇shapef_2DLFE(quadr_matrix)
    
    quadr_vect = Q2_ref
    phi_val_vector = shapef_2DLFE(quadr_vect)
    points_vector = quadr_vect.points

    fill!(Ke, 0)
    fill!(fe, 0)

    weights_vector = quadr_vect.weights
    weights_matrix = quadr_matrix.weights

    for i = 1:3
        for j = 1:3
            # trasposta di invBk per gradiente della iesima phi in TUTTI i punti
            # quindi contiene ... nel primo, ... nel secondo e poi nel terzo
            bktm1_∇phi_i = transpose(invBk)*phi_grad[:, i, :]
            bktm1_∇phi_j = transpose(invBk)*phi_grad[:, j, :]
            for s in 1:size(quadr_matrix.points, 2) 
                Ke[i, j] += bktm1_∇phi_i[:, s] ⋅ bktm1_∇phi_j[:, s]*detBk*weights_matrix[s]    
            end
        end 
        f_cap = (x) -> f(Bk*x+ak)
        # quadratura Q2
        for l in 1:size(points_vector, 2)
            fe[i] += weights_vector[l]*f_cap(points_vector[:, l])*phi_val_vector[i,  l] * detBk
        end
    end 
    return Ke, fe
end