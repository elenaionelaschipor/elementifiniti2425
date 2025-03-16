import Pkg
Pkg.activate("elementifinitiunipv_pkg")

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