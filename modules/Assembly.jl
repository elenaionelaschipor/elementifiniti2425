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
            println(i)
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
        for i in 1:size(p, 2)
            gradients[:, :, i] =  [∇_phi_1(p[:, i]); ∇_phi_2(p[:, i]); ∇_phi_3(p[:, i])]
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
    Bk = B[cell_index]
    ak = a[cell_index]
    detBk = detB[cell_index]
    invBk = invB[cell_index]

    phi_grad = ∇shapef_2DLFE(Q0_ref)
    phi_val = shapef_2DLFE(Q2_ref)
    points_Q2 = Q2_ref.points
    Ke = zeros(3,3)
    fe = zeros(3)
    for i = 1:3
        for j = 1:3
            K_ij = 0.5 *(Transpose(invBk)*phi_grad[:, :, j]) * (Transpose(invBk)*phi_grad[:, :, i]) * detBk
            Ke[i, j] = Ke[i, j] + K_ij
        end 
        f_cap = (x) -> f(Bk*x+ak)
        int_part = 0.5/3*(Transpose(invBk)*f_cap(points_Q2[:, i])) * (Transpose(invBk)*phi_val[:,  i]) * detBk
            
        fe[i] = fe[i] +  int_part
         
    end 
    return Ke, fe
end