import Pkg
Pkg.activate("C:/Users/elena_/Documents/GitHub/elementifiniti2425/elementifinitiunipv_pkg")
using LinearAlgebra
using SparseArrays
using LaTeXStrings
using Plots
using Printf

function f(x)
    return (pi^2)*cos.(pi*x)
    end

function u(x)
    return cos.(pi*x)
    end


function poisson1d(N, f, ga, gb)
# build A
    h = 1/N
    # D = [0 => 1/(h^2)*2*ones(N-1), 1 => -1/(h^2)*ones(N-2), -1 => -1/(h^2)*ones(N-2) ]
    A = spdiagm(0 => 1/(h^2)*2*ones(N-1), 1 => -1/(h^2)*ones(N-2), -1 => -1/(h^2)*ones(N-2)
    )
    x = collect(range(0+h,1-h, step = h))
    b = f(x) + vcat([ga/h^2], vcat(zeros(N-3), [gb/h^2]))
    # println(b)
    # println(vcat(zeros(N-3), [4]))
    #return A, b
    uh = A\b
    return uh, x
    end



function compute_error(u, uh; N = 100)
    h = 1/N
    x = collect(range(0+h,1-h, step = h))
    u_esatta = u(x)

    error = maximum(abs.(u_esatta-uh))
    return error
    end


N = [10,20,40,80,160,320]
err = []
for n in N
push!(err, compute_error(u,poisson1d(n, f, 1, -1)[1];N=n))
end

begin
# plots
plot()
for n in N
uh, x = poisson1d(n, f, 1, -1)
scatter!(x, uh, ms = 4, ma = 0.5, label = L"$u_h$, N = %$(n)$")
end
h = 1/100
l = collect(range(0+h,1-h, step = h))
plot!(l, u(l), lw =2 , label = L"$u$")
xlims!(0, 1)
title!(L"Soluzione esatta $u(x) = cos(\pi x)$ vs soluzione approssimata $u(x)$", titlefontsize = 10)

end


plot(N, err, xscale=:log10, yscale=:log10, label = "errore", lw = 2, ms = 2)
scatter!(N, err,xscale=:log10, yscale=:log10, label = "", ms = 4)
plot!(N, 1 ./N, xscale=:log10, yscale=:log10, label =L"riferimento $1/N$")
plot!(N, 1 ./N.^2, xscale=:log10, yscale=:log10, label =L"riferimento $\frac{1}{N^2}$")

title!("Errore")
xlabel!(L"$N$")
ylabel!(L"$||u - u_h||_{\infty}$")

using Revise
includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Meshing.jl")


begin
    
mesh_square(0.3)
mesh_circle(0.3)


# cambiare eventualmente a square_003 o circle 003 per avere mesh più piccole (ci ho messo 5 minuti a fare i conti non rifacciamoli)
T_circle_ind, T_circle_points= get_nodes_connectivity("tmp_circle.msh")

T_square_ind, T_square_points = get_nodes_connectivity("tmp_square.msh")

b_circle_ind, b_circle_p = get_boundary_nodes("tmp_circle.msh")
b_square_ind, b_square_p = get_boundary_nodes("tmp_square.msh")


function plot_mesh(T_ind, T_points, b_ind, b_p)
    p = plot(axis = false, ratio = 1)
    N_tri = size(T_ind, 2)
    N_points = size(T_points, 2)
    N_bounds = size(b_ind,1)
    # list_points_x = []
    # list_points_y = []
    boundary_nodes = []
    for j in 1:N_bounds
        push!(boundary_nodes, b_p[:, j])
    end
    for i in 1:N_tri
        i1, i2, i3 = T_ind[:, i]
        
        x1, y1 = T_points[:, i1]
        x2, y2 = T_points[:, i2]
        x3, y3 = T_points[:, i3]
        # append!(list_points_x, [x1, x2, x3])
        # append!(list_points_y, [y1,y2,y3])
        plot!([x1,x2], [y1, y2], label="", lw=2, color=:green)
        plot!([x2, x3], [y2,y3], label="", lw=2, color=:green)
        plot!([x3,x1], [y3,y1], label="", lw=2, color=:green)

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
    # plot(list_points_x, list_points_y, label = "", color=:green)
    # return Nothing
end


plot_mesh(T_circle_ind, T_circle_points, b_circle_ind, b_circle_p)
plot_mesh(T_square_ind, T_square_points, b_square_ind, b_square_p)

end

import Meshes
mesh = to_Meshes(T_circle_ind, T_circle_points)
Meshes.viz(mesh, showsegments = true)