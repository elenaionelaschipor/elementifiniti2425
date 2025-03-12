# Author: Ivan Bioli (https://github.com/IvanBioli)

import Pkg
Pkg.activate("elementifinitiunipv_pkg")
using Revise
includet("../../modules/Meshing.jl")
using Plots

# Function to plot the mesh
function plot_mesh(T, p)
    fig = plot()
    # Loop through each triangle in the mesh
    for i in 1:size(T, 2)  # Now iterating over columns (triangles)
        # Extract the coordinates of the vertices of the current triangle
        v1 = p[:, T[1, i]]  # Vertex 1 (using column index from the 1st row of T)
        v2 = p[:, T[2, i]]  # Vertex 2 (using column index from the 2nd row of T)
        v3 = p[:, T[3, i]]  # Vertex 3 (using column index from the 3rd row of T)
        
        # Create arrays of x and y coordinates for the triangle
        x = [v1[1], v2[1], v3[1], v1[1]]  # x-coordinates (closed loop)
        y = [v1[2], v2[2], v3[2], v1[2]]  # y-coordinates (closed loop)
        
        # Plot the triangle by connecting the vertices
        plot!(x, y, label="", linewidth=1, color=:black)
    end
    return fig 
end

# Create the mesh
h = 0.05
out_file = mesh_square(h)
T, p = get_nodes_connectivity(out_file)

# Plot the mesh
fig = plot_mesh(T, p)
plot!(fig, aspect_ratio=1)
savefig("figures/ex01_2.pdf")
display(fig)

# Compare with Meshes.jl
import Meshes
mesh = to_Meshes(T, p)
Meshes.viz(mesh, showsegments = true)
