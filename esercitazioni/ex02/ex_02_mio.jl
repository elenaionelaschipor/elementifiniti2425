import Pkg
Pkg.activate("C:/Users/elena_/Documents/GitHub/elementifiniti2425/elementifinitiunipv_pkg")
using LinearAlgebra
using SparseArrays
using LaTeXStrings
using Plots
using Revise
includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Meshing.jl")

# on meshes\
begin
        
    mesh_square(0.3, display = true)
    mesh_circle(0.3, display = true)
end

begin
    # cambiare eventualmente a square_003 o circle 003 per avere mesh più piccole (ci ho messo 5 minuti a fare i conti non rifacciamoli)
    T_circle_ind, T_circle_points= get_nodes_connectivity("tmp_circle.msh")

    T_square_ind, T_square_points = get_nodes_connectivity("tmp_square.msh")

    b_circle_ind, b_circle_p = get_boundary_nodes("tmp_circle.msh")
    b_square_ind, b_square_p = get_boundary_nodes("tmp_square.msh")


    function plot_mesh(T_ind, T_points, b_ind, b_p)
        p = plot(axis = false, ratio = 1)
        N_tri = size(T_ind, 2)
        N_bounds = size(b_ind,1)
        boundary_nodes = [b_p[:, j] for j in 1:N_bounds]
        for i in 1:N_tri
            i1, i2, i3 = T_ind[:, i]
            x1, y1 = T_points[:, i1]
            x2, y2 = T_points[:, i2]
            x3, y3 = T_points[:, i3]
            plot!([x1,x2,x3,x1], [y1, y2, y3, y1], label="", lw=2, color=:green)
            if [x1, y1] in boundary_nodes && [x2, y2] in boundary_nodes
                plot!([x1,x2], [y1, y2], label="", lw=2, color=:blue)
            end
            if [x2, y2] in boundary_nodes && [x3, y3] in boundary_nodes
                plot!([x2,x3], [y2, y3], label="", lw=2, color=:blue)
            end
            if [x3, y3] in boundary_nodes && [x1, y1] in boundary_nodes
                plot!([x3,x1], [y3,y1], label="", lw=2, color=:blue)
            end 
        end
        display(p)
    end
    plot_mesh(T_circle_ind, T_circle_points, b_circle_ind, b_circle_p)
    plot_mesh(T_square_ind, T_square_points, b_square_ind, b_square_p)
end

T_003, P_003 =get_nodes_connectivity("circle_003.msh") 
b_ind_003, b_p_003 = get_boundary_nodes("circle_003.msh")
plot_mesh(T_003, P_003, b_ind_003, b_p_003)

import Meshes
mesh = to_Meshes(T_circle_ind, T_circle_points)
Meshes.viz(mesh, showsegments = true)

##### QUADRATURE

begin
    function triarea(V1, V2, V3)
        area = 0.5 * abs(V1[1] * (V2[2] - V3[2]) + V2[1] * (V3[2] - V1[2]) + V3[1] * (V1[2] - V2[2]))
        return area
    end

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

    function Q0(u, T_ind, T_points)
        baricentro =  b(T_ind, T_points)
        N_tri = size(T_ind, 2)
        Q = 0
        for i in 1:N_tri
            i1, i2, i3 = T_ind[:, i]
            V1 = T_points[:, i1]
            V2 = T_points[:, i2]
            V3 = T_points[:, i3]
            area_T = triarea(V1, V2, V3) 
            Q += area_T * u(baricentro[i, :])
        end
        return Q
    end

    function Q1(u, T_ind, T_points)
        N_tri = size(T_ind, 2)
        Q = 0
        for i in 1:N_tri
            i1, i2, i3 = T_ind[:, i]
            V1 = T_points[:, i1]
            V2 = T_points[:, i2]
            V3 = T_points[:, i3]
            area_T = triarea(V1, V2, V3) 

            Q += area_T * ( u(V1) + u(V2) + u(V3))/3
        end
        return Q
    end

    function Q2(u, T_ind, T_points)
        N_tri = size(T_ind, 2)
        Q = 0
        for i in 1:N_tri
            i1, i2, i3 = T_ind[:, i]
            V1 = T_points[:, i1]
            V2 = T_points[:, i2]
            V3 = T_points[:, i3]
            m1 = (V1+V2)/2
            m2 = (V3+V2)/2
            m3 = (V1+V3)/2
            area_T = triarea(V1, V2, V3) 
            Q += area_T * ( u(m1) + u(m2) + u(m3))/3
        end
        return Q
    end

    function quadratura(u, T_ind, T_points, i)
        if i == 0
            return Q0(u, T_ind, T_points)
        elseif i == 1
            return Q1(u, T_ind, T_points)
        elseif i == 2
            return Q2(u, T_ind, T_points)
        end
    end
end


u1 = (x) -> exp(x[1]+x[2])
u2 = (x) ->  x[1]^2 + x[2]^2
u3 = (x) -> exp(x[1]^2+x[2]^2)
u4 = (x) -> 1

# test with specific case
i_es= [exp(2)-2*exp(1)+1, 2/3, pi *(exp(1)-1), pi] 

u_n= [u1, u2, u3, u4]

f = (x) -> 2*x

errors = zeros(10, 3, 4)
# H = LinRange(10^-2, 1, 15)
H = 10 .^ range(-2, 0, length=10)
for j in 1:10 
    h = H[j]
    mesh_square(h)
    mesh_circle(h)
    T_ind_S, T_points_S = get_nodes_connectivity("tmp_square.msh")
    T_ind_C, T_points_C = get_nodes_connectivity("tmp_circle.msh")
    
    for l in 1:4
        if l in 1:2
            z = l
            u = u_n[z]
            integrale_esatto = i_es[l]
            for i in 0:2
                q = quadratura(u, T_ind_S, T_points_S, i)
                errors[j, i+1, l] = abs(q-integrale_esatto)
            end
        else
            z = l-2
            u = u_n[l]
            integrale_esatto = i_es[l]
            for i in 0:2
                q = quadratura(u, T_ind_C, T_points_C, i)
                errors[j, i+1, l] = abs(q-integrale_esatto)
            end
        end
    end
    println(h, "-------------------------")
end

etichette = [L"Errore con $u(x,y)= e^{x+y}$ sul quadrato", L"Errore con $u(x,y)=x^2+y^2$ sul quadrato", 
L"Errore con $u(x,y)=e^{x^2+y^2}$ sul cerchio", L"Errore con $u(x,y) = 1$ sul cerchio"]
for z in 1:4
    p = plot(legend = :outerbottomright)
    for i in 0:2
        x = H
        y = errors[:, i+1, z]
        inds = y.>0
        plot!(x[inds], y[inds], xscale = :log10, yscale=:log10, label = L"errore di $Q_%$(i)$", ls=:solid, marker=:circle, lw = 2, ms = 5)
        title!(etichette[z])
    end
    plot!(H, H.^2,  xscale = :log10, yscale=:log10, label = L"riferimento $h^2$", lw = 2, ms = 2, ls = :dash)
    
    if z in 1:2
        plot!(H, exp(-5).*H.^3 , xscale = :log10, yscale=:log10, label = L"riferimento $h^3$", lw = 2, ms = 2, ls = :dash)
    end
    display(p)
end


