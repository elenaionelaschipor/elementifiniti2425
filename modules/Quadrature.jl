# Author: Ivan Bioli (https://github.com/IvanBioli)

"""
    triarea(V1, V2, V3)

Calculate the area of a triangle given its vertices.

# Arguments
- `V1`: The first vertex of the triangle.
- `V2`: The second vertex of the triangle.
- `V3`: The third vertex of the triangle.

# Returns
- `area::Float64`: The area of the triangle.
"""
function triarea(V1, V2, V3)
    area = 0.5 * abs(V1[1] * (V2[2] - V3[2]) + V2[1] * (V3[2] - V1[2]) + V3[1] * (V1[2] - V2[2]))
    return area
end

"""
    b(T, p):

    returns array of baricenter points of the triangles
"""
function b(T_ind, T_points)
    N_tri = size(T_ind, 2)
    ba = zeros(N_tri, 2)
    for i in 1:N_tri
        i1, i2, i3 = T_ind[:, i]
        
        v1 = T_points[:, i1]
        v2 = T_points[:, i2]
        v3= T_points[:, i3]
        
        ba[i,:] = (v1+v2+v3)./3
    end
    return ba
end



"""
    Q0(u, T, p)

Perform numerical integration using the Q0 quadrature rule (i.e., baricenter formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""

function Q0(p, T, u)
    n_tri = size(T, 2) # Number of triangles
    I_approx::Float64 = 0
    for i = 1:n_tri
        # Get the vertices
        V1 = p[:, T[1, i]]
        V2 = p[:, T[2, i]]
        V3 = p[:, T[3, i]]
        # Barycenter
        G = (V1 + V2 + V3) ./ 3
        # Contribution of the triangle of vertices V1, V2, V3
        I_approx += u(G) * triarea(V1, V2, V3)
    end
    return I_approx
end

"""
    Q1(u, T, p)

Perform numerical integration using the Q1 quadrature rule (i.e., vertex formula) over a mesh.
This quadrature rule has order 1.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.

# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""

function Q1(p, T, u)
    n_tri = size(T, 2) # Number of triangles
    I_approx::Float64 = 0
    for i = 1:n_tri
        # Get the vertices
        V1 = p[:, T[1, i]]
        V2 = p[:, T[2, i]]
        V3 = p[:, T[3, i]]
        # Contribution of the triangle of vertices V1, V2, V3
        I_approx += (u(V1) + u(V2) + u(V3)) / 3 * triarea(V1, V2, V3)
    end
    return I_approx
end


"""
    Q2(p, T, u)

Perform numerical integration using the Q2 quadrature rule (i.e., midpoints rule) over a mesh.
This quadrature rule has order 2.

# Arguments
- `p::Matrix`: The coordinates of the mesh nodes.
- `T::Matrix`: The connectivity matrix of the mesh elements.
- `u::Function`: The function to be integrated.
# Returns
- `I_approx::Float64`: The approximate integral of the function over the mesh.
"""

function Q2(p, T, u)
    n_tri = size(T, 2) # Number of triangles
    I_approx::Float64 = 0
    for i = 1:n_tri
        # Get the vertices
        V1 = p[:, T[1, i]]
        V2 = p[:, T[2, i]]
        V3 = p[:, T[3, i]]
        # Mid-Points
        M1 = (V1 + V2) ./ 2
        M2 = (V1 + V3) ./ 2
        M3 = (V2 + V3) ./ 2
        # Contribution of the triangle of vertices V1, V2, V3
        I_approx += (u(M1) + u(M2) + u(M3)) / 3 * triarea(V1, V2, V3)
    end
    return I_approx
end

"""
quadratura(u, T, P, i)
performs either Q0 or Q1 or Q2 depending on the index i given
"""

function quadratura(u, T_ind, T_points, i)
    if i == 0
        return Q0(u, T_ind, T_points)
    elseif i == 1
        return Q1(u, T_ind, T_points)
    elseif i == 2
        return Q2(u, T_ind, T_points)
    end
end