# Author: Ivan Bioli (https://github.com/IvanBioli)
# Inspired by code written by Jochen Hinz (https://github.com/JochenHinz) for MATH-451 @ EPFL

"""
    struct TriQuad

A structure representing a triangular quadrature rule.

# Fields
- `name::String`: The name of the quadrature rule.
- `order::Integer`: The order of the quadrature rule.
- `points::Matrix`: The quadrature points.
- `weights::Array`: The quadrature weights.
"""
struct TriQuad
    name::String
    order::Integer
    points::Matrix
    weights::Array
end

Q0_ref = TriQuad("Q0", 1, reshape([1/3; 1/3], 2,1), [0.5] )

Q1_ref = TriQuad("Q1", 1, reshape([0,1,0, 0,0,1], 2, 3), [1/6, 1/6, 1/6])

Q2_ref = TriQuad("Q2", 2, reshape([0, 0.5, 0.5, 0.5, 0, 0.5], 2, 3), [1/6, 1/6, 1/6])


"""
    Quadrature(u, mesh::Mesh, ref_quad::TriQuad)

Perform numerical integration of a function over a mesh using a given quadrature rule.

# Arguments
- `u`: The function to be integrated.
- `mesh::Mesh`: The mesh object.
- `ref_quad::TriQuad`: The reference quadrature rule.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""
function Quadrature(u, mesh::Mesh, ref_quad::TriQuad)
    w_cap = ref_quad.weights
    p_cap = ref_quad.points
    # println(p_cap)
    q = size(p_cap, 2)
    Bk, ak = get_Bk!(mesh)
    detbk = get_detBk!(mesh)
    
    Q_k = zeros(size(mesh.T, 2))

    for k in 1:size(mesh.T, 2)
        
        A = [Bk[:, :, k]*p_cap[:,i] + ak[:, k] for i in 1:size(p_cap,2)]
        p =  hcat(first.(A), last.(A))'
        println(p)
        # fix here
        # p = reshape(F_k(p_cap), 2, size(p_cap, 2))
        # println(u(p))
        # println(k, size(w_cap), size(u(p)))

        Q_k[k] = dot(w_cap, u(p))*detbk[k]
    end
    return sum(Q_k)
end

# Evaluation of a function
"""
    eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)

Evaluate a function at given points within an element.

# Arguments
- `u::Function`: The function to be evaluated.
- `points_elem::Matrix`: The points at which to evaluate the function.
- `mesh::Mesh`: The mesh object (ignored).
- `tri_idx::Integer`: The index of the current element (ignored).
- `quadrule::TriQuad`: The quadrature rule (ignored).

# Returns
- `u_evals::Matrix`: The evaluated function values at the given points.
"""
function eval_u(u::Function, points_elem::Matrix, mesh::Mesh, tri_idx::Integer, quadrule::TriQuad)
    ######################
    ### COMPLETARE QUI ###
    ######################
end

