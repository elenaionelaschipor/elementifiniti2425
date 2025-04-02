import Pkg
Pkg.activate("C:/Users/elena_/Documents/GitHub/elementifiniti2425/elementifinitiunipv_pkg")
using LinearAlgebra
using SparseArrays
using LaTeXStrings
using Plots
using Revise
includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Meshing.jl")
includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Quadrature_adv.jl")
# define functions and exact integrals
begin 
    u1 = (x) -> exp.(x[1, :]+x[2, :])
    u2 = (x) ->  x[1, :].^2 + x[2, :].^2
    u3 = (x) -> exp.(x[1, :].^2+x[2, :].^2)
    u4 = (x) -> ones(size(x,2))

    # test with specific case
    i_es= [exp(2)-2*exp(1)+1, 2/3, pi *(exp(1)-1), pi] 

    u_n= [u1, u2, u3, u4]
end 



# try with one value of h
methods = [Q0_ref, Q1_ref, Q2_ref]


errors = zeros( 3, 4)
# begin
    h = 0.7
    mesh_square(h)
    mesh_circle(h)
    T_sq, p_sq = get_nodes_connectivity("tmp_square.msh")
    T_circ, p_circ = get_nodes_connectivity("tmp_circle.msh")
    errors = zeros(3, 4)


    meshhhh = Mesh(T_sq, p_sq)
    Quadrature(u1, meshhhh, Q0_ref)
    Bk, ak = get_Bk!(meshhhh)
    
    for l in 1:4
        u = u_n[l]
        integrale_esatto = i_es[l]
        # println("cambio funzione -----------------")
        if l in 1:2  # square
            mesh_sq = Mesh(T_sq, p_sq)
            for i in 1:3
                # println("cambio quadratura ----------------- ")
                q = Quadrature(u, mesh_sq, methods[i])
                errors[i, l] = abs(q-integrale_esatto)
            end
        else        # circle
            mesh_circ = Mesh(T_circ,p_circ)
            for i in 1:3  # q0, q1, q2...
                # println("cambio quadratura ----------------- ")
                q = Quadrature(u, mesh_circ, methods[i])
                errors[i, l] = abs(q-integrale_esatto)
            end
        end
    end
# end


# multiple values of h 

H = 10 .^ range(-2, 0, length=5)

methods = [Q0_ref, Q1_ref, Q2_ref]

errors = zeros(5, 3, 4)
for j in 1:5 
    h = H[j]
    mesh_square(h)
    mesh_circle(h)
    T_ind_S, T_points_S = get_nodes_connectivity("tmp_square.msh")
    T_ind_C, T_points_C = get_nodes_connectivity("tmp_circle.msh")
    

    for l in 1:4
        u = u_n[l]
        integrale_esatto = i_es[l]
        if l in 1:2  # square
            mesh_sq = Mesh(T_ind_S, T_points_S)
            for i in 1:3
                q = Quadrature(u, mesh_sq, methods[i])
                errors[j, i, l] = abs(q-integrale_esatto)
            end
        else        # circle
            mesh_circ = Mesh(T_ind_C, T_points_C)
            for i in 1:3  # q0, q1, q2...
                q = Quadrature(u, mesh_circ, methods[i])
                errors[j, i, l] = abs(q-integrale_esatto)
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


