begin
    import Pkg
    Pkg.activate("C:/Users/elena_/Documents/GitHub/elementifiniti2425/elementifinitiunipv_pkg")
    using Revise
    using Plots
    using LinearAlgebra

    # Load the necessary files
    includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Meshing.jl")
    includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Quadrature_adv.jl")
    includet("C:/Users/elena_/Documents/GitHub/elementifiniti2425/modules/Assembly.jl")
end


# build the mesh with size h
begin
    h = 0.5
    out_file = mesh_circle(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    boundary_tags, boundary_coords = get_boundary_nodes(out_file)
    set_dirichletdofs!(msh, boundary_tags)
end

# Assembly
begin
    f(x) = 1
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, f)
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
end

# ex 2: the ker is empty, ex 3: solve homogeneus problem
begin  
    D = boundary_tags
    n_points = size(p, 2)
    n_tri = size(T, 2)
    N = [i for i in 1:n_points]
    F = filter((x) -> !(x in D), N)
    ker = nullspace(Matrix(A)) # non è banale!!!
    A_cond = A[F, F]
    ker_cond = nullspace(Matrix(A_cond)) #  questo è banale

    b_cond = Vector(b[F])

    uF = A_cond\b_cond
    uD = zeros(size(D)[1])

    u_h = zeros(n_points)
    u_h[F] = uF
    u_h[D] = uD
end
# plots
begin
    plot_flat(msh, u_h)
    plot_surf(msh, u_h)
end

# controllo della norma infinito
begin
    u_ex = (x) -> 0.25*(1 .- x[1, :].^2 .- x[2, :].^2)
    u_ex_on_points = u_ex(p)   

    norm_inf = maximum(abs.(u_ex_on_points .- u_h))

end

# norma L2
begin
    quadratura = Q0_ref
    L2error(u_ex,u_h,msh,quadratura)
    
end

begin
    errs_norminf = []
    errs_norml2 = []
    # controllo norma infinito al variare di h
    
    H = 10 .^ range(-2, 0, length=10)
    for h in H
        out_file = mesh_circle(h)
        T, p = get_nodes_connectivity(out_file)
        msh = Mesh_constructor(T, p)

        boundary_tags, boundary_coords = get_boundary_nodes(out_file)
        set_dirichletdofs!(msh, boundary_tags)
        
        initialize_assembly!(msh)
        local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, f)
        # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
        A, b = assemble_global(msh, local_assembler)
        D = boundary_tags
        n_points = size(p, 2)
        n_tri = size(T, 2)
        N = [i for i in 1:n_points]
        F = filter((x) -> !(x in D), N)
        A_cond = A[F, F]
        b_cond = Vector(b[F])
        uF = A_cond\b_cond
        uD = zeros(size(D)[1])
        u_h = zeros(n_points)
        u_h[F] = uF
        u_h[D] = uD
        u_ex = (x) -> 0.25*(1 .- x[1, :].^2 .- x[2, :].^2)
        u_ex_on_points = u_ex(p)   

        append!(errs_norminf, maximum(abs.(u_ex_on_points .- u_h)))
        
        append!(errs_norml2, L2error(u_ex,u_h,msh,quadratura))
    end
end


# plotting

begin
    scatter(H, errs_norminf, label = "errore al variare di h " , yscale=:log10, xscale =:log10)
    plot!(H, H , xscale =:log10, yscale=:log10, label = "rif primo ordine")
    xaxis!("h")
    yaxis!("Errore")
    title!("Errore in norma inf")
end

begin
    scatter(H, errs_norml2, label = "errore al variare di h " , yscale=:log10, xscale =:log10)
    plot!(H, H.^2 , xscale =:log10, yscale=:log10, label = "rif secondo ordine")
    xaxis!("h")
    yaxis!("Errore")
    title!("Errore in norma L2")
end

    
# -------------------------------------- CONDIZIONI NON OMOGENEE ---------------------------------------------------
begin
    # test
    
    g = (x) -> - exp(-x[1]^2 + x[2]^2)
    func = (x) -> - exp(-x[1]^2 + x[2]^2)*(4*x[1]^2 + 4*x[2]^2)
    h = 0.03
    out_file = mesh_circle(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    boundary_tags, boundary_coords = get_boundary_nodes(out_file)
    set_dirichletdofs!(msh, boundary_tags)

    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = poisson_assemble_local!(Ke, fe, msh, cell_index, func)
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)


    A_cond, b_cond, u_h = impose_dirichlet(A, b, g, msh)


    u_ex = (x) -> - exp.(-x[1, :].^2 .+ x[2, :].^2)
    quadratura = Q2_ref
    println("errore:", L2error(u_ex,u_h,msh,quadratura))
end





