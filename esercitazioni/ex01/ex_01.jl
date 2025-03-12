begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
end

begin    # così runna tutte le linee di codice in mezzo
x = 1
y= 2
x+y   # quando fai ctrl + enter ti stampa a terminale 

end

begin
    function example(x, y=10; z=20)
        return x+y+z
    end
    x = 5
    example(4, 100, z=6)
end

begin
    using LinearAlgebra
    using SparseArrays
    using LaTeXStrings
    using Plots
    using 
end

begin 
    function f(x)
        return pi^2*cos(pi*x)
    end
    f(pi)
end
begin
    function u(x)
        return cos(pi*x)
    end
end

begin
    function poisson1d(N, f, ga, gb)
        # costruire matrice delle diff finite
        h = 1/N
        x = range(0+h, 1-h, N-1)
        f_ = f(x)

        A = spdiagsm(-1 => -1/(h^2)*ones(N-3), 0=> 2/h^2*ones(N-2), 1 => -1/(h^2)*ones(N-3))
        b = f_x + hcat([1/h^2 * ga], zeros(N-3), [1/h^2 * gb]) 
        
        u_h = A\b
        return u_h
    end



poisson1d(4, f, 1, 1)    
end


