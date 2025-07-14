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


# aggiungiamo cond di neumann....
begin
    h = 0.1
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["lower"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    N_tags = setdiff(boundary_tags, D_tags)
    F_tags = setdiff(1:get_ndofs(msh), boundary_tags)

    # assembly
    f(x) = 1
    K = (x) -> 1
    beta = (x) -> [1., 0.]
    initialize_assembly!(msh)  
    local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta)
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    g = (x) -> 0 # x[1] + x[2]
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g, D_tags)
    plot_surf(msh, u_h)
    
end



# il problema diventa di traporto!!! (facciamo solo dirichlet)
begin
    h = 0.1
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    N_tags = setdiff(boundary_tags, D_tags)  # should be empty
    F_tags = setdiff(1:get_ndofs(msh), D_tags)

    # assembly
    f(x) = 1
    K = (x) -> 1
    
    beta = (x) -> [0., -15.]
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta)
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    g = (x) -> 0 #  0.5x[1] + exp(-x[2])
    
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g, D_tags)
    plot_surf(msh, u_h)
    
end


# aggiungiamo k 

begin
    h = 0.1
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["right"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    N_tags = setdiff(boundary_tags, D_tags)  # should be empty
    F_tags = setdiff(1:get_ndofs(msh), D_tags)

    # assembly
    f(x) = 1
    K = (x) -> 0.01
    
    beta = (x) -> [0., 10.]
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) =  transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta)
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    g = (x) -> 0
    
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g, D_tags)
    plot_surf(msh, u_h)
    
end