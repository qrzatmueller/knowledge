#gmsh api is intalled in the (base) environment and the thesis environment

import gmsh 
geo = gmsh.model.geo
import numpy as np


gmsh.initialize()
gmsh.option.setNumber("Mesh.SaveAll", 1) #so not only physical entities are saved
gmsh.model.add("test")

##--- Mic Array ---
# -- Parameters --
ny, nz = 3, 3
Ly, Lz = 1.0, 1.0
mesh_size = 0.1

# --1. Node generation
node_tags = {}
for i in range(ny):
    for j in range(nz):
        y = (i) * (Ly /(ny - 1)) - Ly/2 
        z = (j) * (Lz /(nz - 1)) - Lz/2
        tag = geo.addPoint(0,y,z,mesh_size)
        node_tags[(i,j)] = tag 

#2 --Connect Horizontal Lines
for j in range(nz):
    for i in range(ny-1):
        p1 = node_tags[(i,j)]
        p2 = node_tags[(i+1, j)]
        geo.addLine(p1,p2)

#3 --Connect Vertical Lines
for i in range(ny):
    for j in range(nz-1):
        p1 = node_tags[(i,j)]
        p2 = node_tags[(i,j+1)]
        geo.addLine(p1,p2)

##--SpeakerFrame
#--parameters
x_S = -1
Ly_S = 1
Lz_S = 0.5
ny_S = 2
nz_S = 2
# --1. Node generation
node_tags_S = {}
for i in range(ny_S):
    for j in range(nz_S):
        y = (i) * (Ly_S /(ny_S-1)) -Ly_S/2
        z = (j) * (Lz_S /(ny_S-1) )- Lz_S
        tag = geo.addPoint(x_S,y,z,mesh_size)
        node_tags_S[(i,j)] = tag 

#2 --Connect Horizontal Lines
for j in range(nz_S):
    for i in range(ny_S-1):
        p1 = node_tags_S[(i,j)]
        p2 = node_tags_S[(i+1, j)]
        geo.addLine(p1,p2)

#3 --Connect Vertical Lines
for i in range(ny_S):
    for j in range(nz_S-1):
        p1 = node_tags_S[(i,j)]
        p2 = node_tags_S[(i,j+1)]
        geo.addLine(p1,p2)

#supports
Ly_sup = Ly + 0.2
level=1

p1 = geo.addPoint(0, -Ly_sup/2, 0, mesh_size)
p2 = node_tags[(0,level)]
p3 = node_tags_S[(0,1)]
p4 = node_tags_S[(0,0)]
geo.addLine(p1,p2)
geo.addLine(p1,p3)
geo.addLine(p1,p4)


p5 = geo.addPoint(0, Ly_sup/2, 0, mesh_size)
p6 = node_tags[(2,level)]
geo.addLine(p5,p6)





geo.synchronize()

# gmsh.model.addPhysicalGroup(1, [1], 5)


# gmsh.model.mesh.generate(1) #for meshing

gmsh.write("micarray.geo_unrolled")
gmsh.finalize()
print("Gmsh ran successfully and wrote a .msh / .geo file")
