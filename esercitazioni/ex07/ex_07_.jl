
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
    import Meshes
end


begin
    h = 0.05
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["left", "right"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    F_tags = setdiff(1:get_ndofs(msh), D_tags)

    # assembly
    f(x) = 1
    K = (x) -> 0.001 #10^-4
    beta = (x) -> [1,0]
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta; δ = 0.5 )
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    g = (x) -> 0
    
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g, D_tags)
    plot_surf(msh, u_h)
end


# layer caratteristico

begin
    stab = nothing
    h = 0.05
    out_file = mesh_square(h)
    T, p = get_nodes_connectivity(out_file)
    msh = Mesh_constructor(T, p)

    D_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    boundary_tags, _ = get_boundary_nodes(out_file; labels = ["boundary"])
    F_tags = setdiff(1:get_ndofs(msh), D_tags)

    # assembly
    f(x) = 0
    K = (x) -> 0.0001 #10^-4
    beta = (x) -> [2,3]
    initialize_assembly!(msh)
    local_assembler(Ke, fe, msh, cell_index) = transport_assemble_local!(Ke, fe, msh, cell_index, f, K, beta; stab = stab, δ = 0.5 )
    # vuol dire che f è fissato, tutto il resto no quindi va in input delle altre cose
    A, b = assemble_global(msh, local_assembler)
    
    function g_(x)
        if x[2] <= 1-3*x[1]
            return 1
        else
            return 0
        end
    end
    
    A_cond, b_cond, u_h = impose_dirichlet_neumann(A,b,msh,g_, D_tags)
    plot_surf(msh, u_h)



end