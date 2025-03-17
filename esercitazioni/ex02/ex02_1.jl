import Pkg
Pkg.activate("elementifinitiunipv_pkg")
using Revise
includet("../../modules/Meshing.jl")

# Create the mesh
h = 0.05
out_file = mesh_circle(h; display=true)
T, p = get_nodes_connectivity(out_file)
T
p

i = 72
v1 = p[:, T[1, i]]

using Plots
x = [0, 1, 0, 0]
y = [0, 0, 1, 0]
plot!(x, y, label="", linewidth=1, color=:black)