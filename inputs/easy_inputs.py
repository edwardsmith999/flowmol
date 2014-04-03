#! /usr/bin/env python
import math

# Desired "Inputs" 
density = 1.0
rcutoff = 2.0**(1.0/6.0)
dr_nbr  = 0.3
D_xL    = 1562                                 # Prefix "D_" for "desired"
D_zL    = 1071 
D_yL_md = 568
D_yL_cfd = 568
wall_layers = 0
eps = 0.01
# npx_md = 12
# npy_md = 6
# npz_md = 4
# cfd_cells_per_md_proc_x = 7
# cfd_cells_per_md_proc_y = 33
# cfd_cells_per_md_proc_z = 17
npx_md = 12
npy_md = 6
npz_md = 10
cfd_cells_per_md_proc_x = 7
cfd_cells_per_md_proc_y = 33
cfd_cells_per_md_proc_z = 5

ncy_olap = 0


# Calculate lattice properties and size 
a = 1.0 / (math.pow((density/4.0),(1.0/3.0)))    # lattice parameter a
D_yL_md_pluswall = (D_yL_md + 0.25*a + 
                    0.5*(wall_layers - 1))       # account for wall
nunits_x_md = int(math.floor(D_xL/a))            # best number of units
nunits_y_md = int(math.floor(D_yL_md_pluswall/a))
nunits_z_md = int(math.floor(D_zL/a))
teth_dist = a * (0.25 + 0.5*(wall_layers - 1))
teth_dist = round(teth_dist,2) + eps

# Calculate MD domain, proc domain, ncells and cellsize
xL_md = nunits_x_md*a
yL_md = nunits_y_md*a
zL_md = nunits_z_md*a
xL_md_l = xL_md / npx_md
yL_md_l = yL_md / npy_md
zL_md_l = zL_md / npz_md
ncx_md_l = int(math.floor(xL_md_l / (rcutoff + dr_nbr)))
ncy_md_l = int(math.floor(yL_md_l / (rcutoff + dr_nbr)))
ncz_md_l = int(math.floor(zL_md_l / (rcutoff + dr_nbr)))
ncx_md = npx_md * ncx_md_l
ncy_md = npy_md * ncy_md_l
ncz_md = npz_md * ncz_md_l
dx_md = xL_md / float(ncx_md)
dy_md = yL_md / float(ncy_md)
dz_md = zL_md / float(ncz_md)

# Calculate CFD cellsize, ncells, and domain size 
dx_cfd = xL_md_l / cfd_cells_per_md_proc_x
dy_cfd = yL_md_l / cfd_cells_per_md_proc_y
dz_cfd = zL_md_l / cfd_cells_per_md_proc_z
D_xL_cfd = xL_md                                # Constraint
D_yL_cfd = D_yL_cfd
D_zL_cfd = zL_md                                # Constraint
ncx_cfd = int(round(D_xL_cfd/dx_cfd))
ncy_cfd = int(round(D_yL_cfd/dy_cfd))
ncz_cfd = int(round(D_zL_cfd/dz_cfd))
xL_cfd = ncx_cfd*dx_cfd
yL_cfd = ncy_cfd*dy_cfd
zL_cfd = ncz_cfd*dz_cfd

# Other interesting information
globalnp = 4*nunits_x_md*nunits_y_md*nunits_z_md
yL_olap = ncy_olap*dy_cfd

message = (
           "CPL_parameter_tuning inputs: \n" + 
           "\tdensity:           \t" + str(density)        + "\n" + 
           "\trcutoff:           \t" + str(rcutoff)        + "\n" + 
           "\tdr_nbr:            \t" + str(dr_nbr)         + "\n" + 
           "\tDesired xL:        \t" + str(D_xL)           + "\n" + 
           "\tDesired zL:        \t" + str(D_zL)           + "\n" + 
           "\tDesired MD yL:     \t" + str(D_yL_md)        + "\n" + 
           "\tHCP wall layers:   \t" + str(wall_layers)    + "\n" + 
           "\tDesired CFD yL:    \t" + str(D_yL_cfd)       + "\n" + 
           "\tnproc MD (x,y,z):  \t" + str(npx_md)         + "\t" + 
                                       str(npy_md)         + "\t" +
                                       str(npz_md)         + "\n" +
           "\tcfdcells/MDproc:   \t" + str(cfd_cells_per_md_proc_x) + "\t" + 
                                       str(cfd_cells_per_md_proc_y) + "\t" +
                                       str(cfd_cells_per_md_proc_z) + "\n" +
           "\toverlap cells (y): \t" + str(ncy_olap)       + "\n\n" +

           "---------- TUNING CALCULATIONS RESULTS ---------- \n\n" +

           "Results: \n" + 
           "\tMD FCC units(x):   \t" + str(nunits_x_md)    + "\n" + 
           "\tMD FCC units(y):   \t" + str(nunits_y_md)    + "\n" + 
           "\tMD FCC units(z):   \t" + str(nunits_z_md)    + "\n" + 
           "\tFCC lattice param: \t" + str(a)              + "\n" + 
           "\tMD tether dist:    \t" + str(teth_dist)      + "\n" + 
           "\tMD local cells(x): \t" + str(ncx_md_l)       + "\n" +  
           "\tMD local cells(y): \t" + str(ncy_md_l)       + "\n" +  
           "\tMD local cells(z): \t" + str(ncz_md_l)       + "\n" +  
           "\tMD ncells(x):      \t" + str(ncx_md)         + "\n" +  
           "\tMD ncells(y):      \t" + str(ncy_md)         + "\n" +  
           "\tMD ncells(z):      \t" + str(ncz_md)         + "\n\n" +  
  
           "\tMD  xL:            \t" + str(xL_md)          + "\n" + 
           "\tCFD xL:            \t" + str(xL_cfd)         + "\n" + 
           "\tMD  yL:            \t" + str(yL_md)          + "\n" + 
           "\tCFD yL:            \t" + str(yL_cfd)         + "\n" + 
           "\tMD  zL:            \t" + str(zL_md)          + "\n" + 
           "\tCFD zL:            \t" + str(zL_cfd)         + "\n\n" + 
  
           "\tCFD ncells(x):     \t" + str(ncx_cfd)        + "\n" +  
           "\tCFD ncells(y):     \t" + str(ncy_cfd)        + "\n" +  
           "\tCFD ncells(z):     \t" + str(ncz_cfd)        + "\n" +  
		   "\tMD cell length(x): \t" + str(dx_md)          + "\n" +  
           "\tMD cell length(y): \t" + str(dy_md)          + "\n" +  
           "\tMD cell length(z): \t" + str(dz_md)          + "\n" +  
           "\tCFD cell length(x):\t" + str(dx_cfd)         + "\n" +  
           "\tCFD cell length(y):\t" + str(dy_cfd)         + "\n" +  
           "\tCFD cell length(z):\t" + str(dz_cfd)         + "\n\n" + 

           "\tcell2bin_ratio (x):\t" + str(dx_cfd/dx_md)   + "\n" +
           "\tcell2bin_ratio (y):\t" + str(dy_cfd/dy_md)   + "\n" +
           "\tcell2bin_ratio (z):\t" + str(dz_cfd/dz_md)   + "\n" +
           "\tbin2cell_ratio (x):\t" + str(dx_md/dx_cfd)   + "\n" +
           "\tbin2cell_ratio (y):\t" + str(dy_md/dy_cfd)   + "\n" +
           "\tbin2cell_ratio (z):\t" + str(dz_md/dz_cfd)   + "\n" +
           "\tMD proc domain(y): \t" + str(yL_md_l)        + "\n" +           
           "\tO'lap size (y):    \t" + str(yL_olap)        + "\n\n\n" + 


           "----------- SUGGESTED INPUT PARAMETERS ---------- \n\n" +
           "MD.in: \n\n" +
           "\tDENSITY:           \t" + str(density) + "\n" +
		   "\tINITIALNUNITS(x):  \t" + str(nunits_x_md) + "\n" +
           "\tINITIALNUNITS(y):  \t" + str(nunits_y_md) + "\n" +
           "\tINITIALNUNITS(z):  \t" + str(nunits_z_md) + "\n" +
           "\tDELTA_RNEIGHBR:    \t" + str(dr_nbr) + "\n" +
           "\tTETHERDISTBOTTOM:  \t" + str(teth_dist) + "\n" +
           "\tTHERMSTATBOTTOM:   \t" + "don't forget it!" + "\n" +
           "\tPROCESSORS:        \t" + "don't forget it!" + "\n" +
           "\tVELOCITY_OUTFLAG:  \t" + "match timing with CFD save_dt!" + "\n" +
           "\tMASS_OUTFLAG:      \t" + "match timing with CFD save_dt!" + "\n" +
           "\tPRESSURE_OUTFLAG:  \t" + "match timing with CFD save_dt!" + "\n\n" +
           
          "param.inc: \n\n" +
          "\tngx:                \t" + str(ncx_cfd) + "+1\n" +
          "\tngy:                \t" + str(ncy_cfd) + "+1\n" +
          "\tngz:                \t" + str(ncz_cfd) + "+1\n" +
          "\tnp(xyz):            \t" + "must be a factor of nc(xyz)_cfd" + "\n\n" +

          "input: \n\n" 
          "\tRe:                 \t" + "don't forget it!" + "\n" +
          "\tnsteps:             \t" + "don't forget it!" + "\n" +
          "\tsave dt:            \t" + "match timing with MD outputs" + "\n\n" +

          "input.file: \n\n" +
          "\tLx:                 \t" + str(xL_cfd) + "\n" +
          "\tLy:                 \t" + str(yL_cfd) + "\n" +
          "\tngx:                \t" + str(ncx_cfd+1) + "\n" +
          "\tngy:                \t" + str(ncy_cfd+1) + "\n\n" +

          "input.setup: \n\n" +
          "\tnix:                \t" + str(ncx_cfd+1) + "\n" +
          "\tniy:                \t" + str(ncy_cfd+1) + "\n" +
          "\tniz:                \t" + str(ncz_cfd+1) + "\n" +
          "\txL:                 \t" + str(xL_cfd) + "\n" +
          "\tyL:                 \t" + str(yL_cfd) + "\n" +
          "\tzL:                 \t" + str(zL_cfd) + "\n\n" +

          "COUPLER.in: \n\n" +
          "\tDENSITY_CFD:        \t" + str(density) + "\n" +
          "\tTIMESTEP_RATIO:     \t" + "match with MD and CFD dts" + "\n" +
          "\tCONSTRAINT_INFO:    \t" + "don't forget to set the cells correctly!" + "\n" +
          "\tOVERLAP_EXTENTS:    \t" + "don't forget to set the cells correctly! \n\n\n"

          )

print(message)
