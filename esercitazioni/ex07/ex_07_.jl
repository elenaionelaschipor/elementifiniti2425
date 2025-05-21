
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


begin
    h = 0.05
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    F_tags = setdiff(1:get_ndofs(msh), D_tags)

    # assembly
    f(x) = 1
    K = (x) -> 0.0001
    beta = [1,0]
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta; stab = "SUPG", δ = 0.5 )
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    g = (x) -> 0
    
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g, D_tags)
    plot_surf(msh, u_h)
    
end