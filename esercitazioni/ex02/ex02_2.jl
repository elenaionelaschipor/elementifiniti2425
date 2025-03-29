# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg")
using Revise

# Load the necessary files
includet("../modules/Meshing.jl")
includet("../modules/Quadrature.jl")

using Plots
using LaTeXStrings

######################## u(x,y) = exp(x+y)  ########################
u(x) = exp(x[1] + x[2]);
u_str = L"e^{x+y}"
I_exact = (exp(1) - 1)^2

begin
    mesh_sizes = 10 .^ (range(0, -2, 10))
    Q0_approx = similar(mesh_sizes)
    Q1_approx = similar(mesh_sizes)
    Q2_approx = similar(mesh_sizes)
    for (i, h) in enumerate(mesh_sizes)
        # Build the mesh with mesh-size h    
        out_file = mesh_square(h)
        T, p = get_nodes_connectivity(out_file)
        # Compute the integral
        Q0_approx[i] = Q0(p, T, u)
        Q1_approx[i] = Q1(p, T, u)
        Q2_approx[i] = Q2(p, T, u)
    end
end

# Plotting
begin
    plt = plot(mesh_sizes, abs.(Q0_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_0$", marker=:circle)
    plot!(mesh_sizes, abs.(Q1_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_1$", marker=:rect)
    plot!(mesh_sizes, abs.(Q2_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_2$", marker=:diamond)
    plot!(mesh_sizes, mesh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")
    plot!(mesh_sizes, (mesh_sizes .^ 3) .* 1e-2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^3)$")
    plot!(plt, legend=:bottomright)
    title!("Andamento dell'errore con $u_str sul quadrato")
    xlabel!(L"$h$")
    ylabel!("Errore")
    savefig(plt, "./figures/ex02_2_1.pdf") # Save plot
    plot(plt) # Show plot
end

######################## u(x,y) = x^2 + y^2  ########################
u(x) = x[1]^2 + x[2]^2;
u_str = L"x^2 + y^2"
I_exact = 2 / 3

begin
    mesh_sizes = 10 .^ (range(0, -2, 10))
    Q0_approx = similar(mesh_sizes)
    Q1_approx = similar(mesh_sizes)
    Q2_approx = similar(mesh_sizes)
    for (i, h) in enumerate(mesh_sizes)
        # Build the mesh with mesh-size h    
        out_file = mesh_square(h)
        T, p = get_nodes_connectivity(out_file)
        # Compute the integral
        Q0_approx[i] = Q0(p, T, u)
        Q1_approx[i] = Q1(p, T, u)
        Q2_approx[i] = Q2(p, T, u)
    end
end

# Plotting
begin
    plt = plot(mesh_sizes, abs.(Q0_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_0$", marker=:circle)
    plot!(mesh_sizes, abs.(Q1_approx .- I_exact), xaxis=:log, yaxis=:log, label=L"$Q_1$", marker=:rect)
    plot!(mesh_sizes, abs.(Q2_approx .- I_exact) .+ eps(Float64), xaxis=:log, yaxis=:log, label=L"$Q_2$", marker=:diamond)
    plot!(mesh_sizes, mesh_sizes .^ 2, xaxis=:log, yaxis=:log, linestyle=:dash, label=L"$\mathcal{O}(h^2)$")
    plot!(plt, legend=:bottomright)
    title!("Andamento dell'errore con $u_str sul quadrato")
    xlabel!(L"$h$")
    ylabel!("Errore")
    savefig(plt, "./figures/ex02_2_2.pdf") # Save plot
    plot(plt) # Show plot
end
