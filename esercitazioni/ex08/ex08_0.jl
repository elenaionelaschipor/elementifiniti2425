# Author: Ivan Bioli (https://github.com/IvanBioli)
begin
    import Pkg
    Pkg.activate("elementifinitiunipv_pkg")
    using Revise
    
end

begin
    # Load the necessary files
    includet("../../modules/Meshing.jl")
    includet("../../modules/Quadrature_adv.jl")
    includet("../../modules/Assembly.jl")
    includet("../../modules/Assembly_mixed.jl")
    using LinearAlgebra
end

h = 0.1
out_file = mesh_circle(h)
edges2nodes, elems2edges, elems2orientation = get_edges_info(out_file)
msh = Mesh_constructor(get_nodes_connectivity(out_file)...)
set_edges_info!(msh, edges2nodes, elems2edges, elems2orientation)


quad = Q2_ref
quad.points
shapef_2D_RT0FE(quad)

divshapef_2D_RT0FE(quad)


Ae = zeros(3,3)
Be = zeros(1,3)
fe = zeros(1)
f = (x) -> 1*ones(axes(x))
mu = (x) -> 2*ones(2,2)
darcy_assemble_local_mixed!(Ae, Be, fe, msh, 4,f, mu )