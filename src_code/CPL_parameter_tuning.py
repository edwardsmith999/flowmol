#!/usr/bin/env python
import math

# "Inputs"
density = 0.8
rcutoff = 2.2
dr_nbr  = 0.29
d_xL    = 41.1
d_zL    = 41.1 
d_yL_md = 31.3
d_yL_cfd = 41.0 
wall_layers = 2
cellsize_ratio = 2.0 # must be a float representation of an integer
eps = 0.01

# Calculate number of initial units required to get close to desired yL_md
a = 1.0 / (math.pow((density/4.0),(1.0/3.0)))    # lattice parameter a
d_yL_md_pluswall = (d_yL_md + 0.25*a + 
                    0.5*(wall_layers - 1))       # account for wall

nunits_x_md = int(math.floor(d_xL/a))            # best number of units
nunits_y_md = int(math.floor(d_yL_md_pluswall/a))# best number of units
nunits_z_md = int(math.floor(d_zL/a))            # best number of units
xL_md = nunits_x_md*a                            # new xL_md
yL_md = nunits_y_md*a                            # new yL_md
zL_md = nunits_z_md*a                            # new zL_md

teth_dist = a * (0.25 + 0.5*(wall_layers - 1))
teth_dist = round(teth_dist,2) + eps

# Calculate computational cell size from rcutoff and yL_md (floor to be
# conservative)
ncells_x_md = int(round(xL_md / (rcutoff + dr_nbr)))
ncells_y_md = int(round(yL_md / (rcutoff + dr_nbr)))
ncells_z_md = int(round(zL_md / (rcutoff + dr_nbr)))
dx_md = xL_md / float(ncells_x_md)
dy_md = yL_md / float(ncells_y_md)
dz_md = zL_md / float(ncells_z_md)

# Calculate cfd cell size and closest possible yL_cfd  

dx_cfd = dx_md*cellsize_ratio
dy_cfd = dy_md*cellsize_ratio
dz_cfd = dz_md*cellsize_ratio

xL_cfd = xL_md  # Constraint
zL_cfd = zL_md  # Constraint

ncells_x_cfd = int(round(xL_cfd/dx_cfd))
ncells_y_cfd = int(round(d_yL_cfd/dy_cfd))
ncells_z_cfd = int(round(zL_cfd/dz_cfd))

yL_cfd = ncells_y_cfd*dy_cfd

message = (
           "Inputs: \n" + 
           "\tdensity:           \t" + str(density)        + "\n" + 
           "\trcutoff:           \t" + str(rcutoff)        + "\n" + 
           "\tdr_nbr:            \t" + str(dr_nbr)         + "\n" + 
           "\tDesired xL:        \t" + str(d_xL)           + "\n" + 
           "\tDesired zL:        \t" + str(d_zL)           + "\n" + 
           "\tDesired MD yL:     \t" + str(d_yL_md)        + "\n" + 
           "\tHCP wall layers:   \t" + str(wall_layers)    + "\n" + 
           "\tDesired CFD yL:    \t" + str(d_yL_cfd)       + "\n" + 
           "\tCellsize ratio:    \t" + str(cellsize_ratio) + "\n" + 

           "Results: \n" + 
           "\tMD FCC units(x):   \t" + str(nunits_x_md)    + "\n" + 
           "\tMD FCC units(y):   \t" + str(nunits_y_md)    + "\n" + 
           "\tMD FCC units(z):   \t" + str(nunits_z_md)    + "\n" + 
           "\tFCC lattice param: \t" + str(a)              + "\n" + 
           "\tMD tether dist:    \t" + str(teth_dist)      + "\n" + 
           "\tMD ncells(x):      \t" + str(ncells_x_md)    + "\n" +  
           "\tMD ncells(y):      \t" + str(ncells_y_md)    + "\n" +  
           "\tMD ncells(z):      \t" + str(ncells_z_md)    + "\n\n" +  
  
           "\tMD  xL:            \t" + str(xL_md)          + "\n" + 
           "\tCFD xL:            \t" + str(xL_cfd)         + "\n" + 
           "\tMD  yL:            \t" + str(yL_md)          + "\n" + 
           "\tCFD yL:            \t" + str(yL_cfd)         + "\n" + 
           "\tMD  zL:            \t" + str(zL_md)          + "\n" + 
           "\tCFD zL:            \t" + str(zL_cfd)         + "\n\n" + 
  
           "\tCFD ncells(x):     \t" + str(ncells_x_cfd)   + "\n" +  
           "\tCFD ncells(y):     \t" + str(ncells_y_cfd)   + "\n" +  
           "\tCFD ncells(z):     \t" + str(ncells_z_cfd)   + "\n" +  
		   "\tMD cell length(x): \t" + str(dx_md)          + "\n" +  
           "\tMD cell length(y): \t" + str(dy_md)          + "\n" +  
           "\tMD cell length(z): \t" + str(dz_md)          + "\n" +  
           "\tCFD cell length(x):\t" + str(dx_cfd)         + "\n" +  
           "\tCFD cell length(y):\t" + str(dy_cfd)         + "\n" +  
           "\tCFD cell length(z):\t" + str(dz_cfd)
          )
print(message)
