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
    F_glob = Matrix(sparse(rows_f, ones(size(rows_f)), data_f))
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
    # Get tags of dirichlet dofs and free dofs
    ndofs = get_ndofs(mesh)
    freedofs, dirichletdofs = get_freedofs(mesh), get_dirichletdofs(mesh)
    # Impose Dirichlet BCs by lifting
    uh = zeros(ndofs)
    uh[dirichletdofs] = dropdims(mapslices(g, mesh.p[:, dirichletdofs]; dims=1); dims=1)
                                    # applico g sui punti
                        # è una matrice 2x1x1... diventa 2x1 (Esempio)
    
    # Modify the system to include Dirichlet BCs
    A_cond = A[freedofs, freedofs]
    b_cond = b[freedofs] - A[freedofs, dirichletdofs] * uh[dirichletdofs]

    return A_cond, b_cond, uh
end

"""
homogeneous neumann, imposed dirichlet
"""
function impose_dirichlet_neumann(A,b,mesh, g, D)
    
    F = setdiff(1:get_ndofs(msh), D)

    A_cond = A[F,F]
    # assembly of vector G
    p = mesh.p
    G = zeros(size(p,2))
    for i in 1:size(p,2)
        G[i] = g(p[:, i])
    end    


    b_cond = Vector(b[F]) - A[F, D]*G[D] 

    uF = A_cond\b_cond
    uD = G[D]

    u_h = zeros(size(p, 2))
    u_h[F] = uF
    u_h[D] = uD

    return A_cond, b_cond, u_h
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
- beta: the transport term

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
    phi_val_matrix= shapef_2DLFE(quadr_matrix)
    
    quadr_vect = Q0_ref
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
            phi_j = phi_val_matrix[j, :]
            for s in 1:size(quadr_matrix.points, 2) # sommo sui punti di quadratura
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

########################### TRANSPORT PROBLEM ###########################
"""
    transport_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k, β; stab = nothing, δ = 0.5)

Assemble the local stiffness matrix and force vector for the transport problem.

# Arguments
- `Ke::Matrix`: The local stiffness matrix to be assembled.
- `fe::Vector`: The local force vector to be assembled.
- `mesh::Mesh`: The mesh object.
- `cell_index::Integer`: The index of the current cell.
- `f`: The source term function.
- `k`: The diffusion coefficient function.
- `β`: The advection velocity function.
- `stab`: The stabilization method (optional).
- `δ`: The stabilization parameter (optional).

# Returns
- `Ke`: The assembled local stiffness matrix.
- `fe`: The assembled local force vector.
"""
# function transport_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k, β; stab = nothing, δ = 0.5)
#     B, a = get_Bk!(mesh)
#     detB = get_detBk!(mesh)
#     invB= get_invBk!(mesh)
#     Bk = B[:, :, cell_index]
#     ak = a[:, cell_index]
#     detBk = detB[cell_index]
#     invBk = invB[:, :, cell_index]


#     quadr_matrix = Q2_ref
#     phi_grad = ∇shapef_2DLFE(quadr_matrix)
#     phi_val_matrix= shapef_2DLFE(quadr_matrix)
    
#     quadr_vect = Q2_ref
#     phi_val_vector = shapef_2DLFE(quadr_vect)
#     points_vector = quadr_vect.points


#     fill!(Ke, 0)
#     fill!(fe, 0)

#     weights_vector = quadr_vect.weights
#     weights_matrix = quadr_matrix.weights

#     for i = 1:3
#         for j = 1:3
#             K_cap = (x) -> K(Bk*x+ak) 
#             # trasposta di invBk per gradiente della iesima phi in TUTTI i punti
#             # quindi contiene ... nel primo, ... nel secondo e poi nel terzo
#             bktm1_∇phi_i = transpose(invBk)*phi_grad[:, i, :]
#             bktm1_∇phi_j = transpose(invBk)*phi_grad[:, j, :]


#             phi_j = phi_val_matrix[j, :]
#             phi_i = phi_val_matrix[i, :]
#             for s in 1:size(quadr_matrix.points, 2) # sommo sui punti di quadratura
#                 beta_p = beta(quadr_matrix.points[:, s])
#                 if isnothing(stab)
#                     Ke[i, j] += (K_cap(quadr_matrix.points[:, s])*bktm1_∇phi_i[:, s]) ⋅ bktm1_∇phi_j[:, s]*detBk*weights_matrix[s] +  beta_p ⋅ bktm1_∇phi_j[:, s] ⋅ phi_i[s]*detBk*weights_matrix[s] 
#                 end
                
#                 h_T = maximum([norm(Bk[:, 1]), norm(Bk[:, 2]), norm(Bk[:, 1] - Bk[:, 2])])
#                 βnormInf_T = maximum(norm.(beta_p, Inf))
#                 eps_h = 0.5*βnormInf_T*h_T
                    
#                 if stab == "NCAD"
#                     Ke[i, j] += (eps_h*bktm1_∇phi_i[:, s]) ⋅ bktm1_∇phi_j[:, s]*detBk*weights_matrix[s] +  beta_p ⋅ bktm1_∇phi_j[:, s] ⋅ phi_i[s]*detBk*weights_matrix[s] 
#                 end
                
#                 if stab == "NCSD"
#                     Ke[i, j] += (K_cap(quadr_matrix.points[:, s])*bktm1_∇phi_i[:, s]) ⋅ bktm1_∇phi_j[:, s]*detBk*weights_matrix[s] +  beta_p ⋅ bktm1_∇phi_j[:, s] ⋅ phi_i[s]*detBk*weights_matrix[s] 
#                     n_beta_p = beta_p ./ norm(beta_p)
#                     Ke[i,j] += eps_h * (n_beta_p ⋅ bktm1_∇phi_i[:, s])  * (n_beta_p ⋅ bktm1_∇phi_j[:, s]) *detBk *weights_matrix[s]
#                 end

#                 if stab == "SUPG"
#                     Ke[i, j] += (K_cap(quadr_matrix.points[:, s])*bktm1_∇phi_i[:, s]) ⋅ bktm1_∇phi_j[:, s]*detBk*weights_matrix[s] +  beta_p ⋅ bktm1_∇phi_j[:, s] ⋅ phi_i[s]*detBk*weights_matrix[s] 
#                     tau_h = δ * h_T/βnormInf_T  # assuming beta_p constant over all omega
#                     Ke[i,j] += tau_h * (beta_p ⋅ bktm1_∇phi_j[:, s])  * (beta_p ⋅ bktm1_∇phi_i[:, s])*detBk*weights_matrix[s] 


#                 end
                
#             end
#         end 
#         f_cap = (x) -> f(Bk*x+ak)
#         # quadratura Q2
#         for l in 1:size(points_vector, 2)
#             fe[i] += weights_vector[l]*f_cap(points_vector[:, l])*phi_val_vector[i,  l] * detBk
#             if stab=="SUPG"
#                 h_T = maximum([norm(Bk[:, 1]), norm(Bk[:, 2]), norm(Bk[:, 1] - Bk[:, 2])])
#                 eps_h = 0.5*norm(beta_p)*h_T
#                 bktm1_∇phi_i = transpose(invBk)*phi_grad[:, i, :]
#                 tau_h = δ * h_T/norm(beta_p)  # assuming beta_p constant over all omega    
#                 fe[i] += tau_h * f_cap(points_vector[:, l]) * beta_p ⋅ bktm1_∇phi_i[:, l]* detBk
#             end
#         end
#     end 
#     return Ke, fe






# end

function transport_assemble_local!(Ke::Matrix, fe::Vector, mesh::Mesh, cell_index::Integer, f, k, β; stab = nothing, δ = 0.5)

    n_basefuncs = 3
    # Reset to 0
    fill!(Ke, 0)
    fill!(fe, 0)
    # We use Q0 quadrule to assemble the stiffness matrix exactly
    quadrule = Q2_ref
    points_e = mesh.Bk[:, :, cell_index] * quadrule.points .+ mesh.ak[:, cell_index]
    # Evaluate basis functions and their gradient
    shapef = shapef_2DLFE(quadrule)
    invBk = mesh.invBk[:, :, cell_index]
    ∇shapef = mapslices(x -> invBk' * x, ∇shapef_2DLFE(quadrule), dims=(1, 2))
    if ~isnothing(stab)
        hₜ = maximum(norm.(eachcol(hcat(mesh.Bk[:, :, cell_index], mesh.Bk[:, 2, cell_index]-mesh.Bk[:, 1, cell_index]))))
    end
    if stab == "SUPG"
        βnormInf_T = maximum(norm.(β.(eachcol(points_e)), Inf))
        τₕ = δ * hₜ / βnormInf_T
    end
    # Loop over quadrature points
    for (q_index, q_point) in enumerate(eachcol(points_e))
        # Get the quadrature weight
        dΩ = quadrule.weights[q_index] * mesh.detBk[cell_index]
        # Get the value of k and β at q_point
        k_eval = k(q_point)
        β_eval = β(q_point) 
        f_eval = f(q_point)
        # Loop over test shape functions
        for i in 1:n_basefuncs
            v = shapef[i, q_index]
            ∇v = ∇shapef[:, i, q_index]
            # Add contribution to fe
            fe[i] += f_eval * v * dΩ
            if stab == "SUPG"
                fe[i] += τₕ * (∇v ⋅ β_eval) * f_eval * dΩ
            end
            # Loop over trial shape functions
            for j in 1:n_basefuncs
                ∇u = ∇shapef[:, j, q_index]
                # Add contribution to Ke
                if stab != "NCAD"
                    Ke[i, j] += (∇v ⋅ (k_eval * ∇u) + (β_eval ⋅ ∇u) * v) * dΩ
                    if ~isnothing(stab)
                        if stab == "NCSD"
                            nβ = β_eval / norm(β_eval)
                            Ke[i, j] += 0.5 * norm(β_eval) * hₜ * (∇v ⋅ nβ) * (∇u ⋅ nβ) * dΩ
                        elseif stab == "SUPG"
                            Ke[i, j] += τₕ * (∇v ⋅ β_eval) * (∇u ⋅ β_eval) * dΩ
                        else
                            error("Unknown stabilization")
                        end
                    end
                else # stab == "NCAD"
                    Ke[i, j] += (0.5 * norm(β_eval) * hₜ * (∇v ⋅ ∇u) + (β_eval ⋅ ∇u) * v) * dΩ
                end
            end
        end
    end
    return Ke, fe
end