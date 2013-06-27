/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/disp,PPPMDisp)

#else

#ifndef LMP_PPPM_DISP_H
#define LMP_PPPM_DISP_H

#include "lmptype.h"
#include "mpi.h"

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_FLOAT
#else
typedef double FFT_SCALAR;
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include "kspace.h"

namespace LAMMPS_NS {


#define EWALD_MAXORDER	6
#define EWALD_FUNCS	3

class PPPMDisp : public KSpace {
 public:
  PPPMDisp(class LAMMPS *, int, char **);
  virtual ~PPPMDisp();
  virtual void init();
  virtual void setup();
  void setup_grid();
  virtual void compute(int, int);
  virtual int timing_1d(int, double &);
  virtual int timing_3d(int, double &);
  virtual double memory_usage();

 protected:

/* ----------------------------------------------------------------------
Variables needed for calculating the 1/r and 1/r^6 potential
------------------------------------------------------------------------- */

  int function[EWALD_FUNCS];

  int me,nprocs;
  int nfactors;
  int *factors;
  double qsum,qsqsum;
  double csumij;
  double csum;
  double *csumi;  //needed as correction term for per atom calculations!
  double *cii;
  int csumflag;
  double cutoff, cutoff_lj;
  double volume;
  double *B;
  double virial_1[6], virial_6[6];
  double sf_coeff[6], sf_coeff_6[6];
  int peratom_allocate_flag;


  double delxinv,delyinv,delzinv,delvolinv;
  double delxinv_6,delyinv_6,delzinv_6,delvolinv_6;
    
  double shift,shiftone;
  int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
  int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
  int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
  int nlower,nupper;
  int ngrid,nfft,nfft_both;

  double shift_6,shiftone_6;
  int nxlo_in_6,nylo_in_6,nzlo_in_6,nxhi_in_6,nyhi_in_6,nzhi_in_6;
  int nxlo_out_6,nylo_out_6,nzlo_out_6,nxhi_out_6,nyhi_out_6,nzhi_out_6;
  int nxlo_fft_6,nylo_fft_6,nzlo_fft_6,nxhi_fft_6,nyhi_fft_6,nzhi_fft_6;
  int nlower_6,nupper_6;
  int ngrid_6,nfft_6,nfft_both_6;

  //// the following variables are needed for every structure factor
  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
  FFT_SCALAR *density_fft;
  FFT_SCALAR ***u_brick;
  FFT_SCALAR ***v0_brick,***v1_brick,***v2_brick,***v3_brick,***v4_brick,***v5_brick;

  FFT_SCALAR ***density_brick_g;
  FFT_SCALAR ***vdx_brick_g,***vdy_brick_g,***vdz_brick_g;
  FFT_SCALAR *density_fft_g;
  FFT_SCALAR ***u_brick_g;
  FFT_SCALAR ***v0_brick_g,***v1_brick_g,***v2_brick_g,***v3_brick_g,***v4_brick_g,***v5_brick_g;

  FFT_SCALAR ***density_brick_a0;
  FFT_SCALAR ***vdx_brick_a0,***vdy_brick_a0,***vdz_brick_a0;
  FFT_SCALAR *density_fft_a0;
  FFT_SCALAR ***u_brick_a0;
  FFT_SCALAR ***v0_brick_a0,***v1_brick_a0,***v2_brick_a0,***v3_brick_a0,***v4_brick_a0,***v5_brick_a0;

  FFT_SCALAR ***density_brick_a1;
  FFT_SCALAR ***vdx_brick_a1,***vdy_brick_a1,***vdz_brick_a1;
  FFT_SCALAR *density_fft_a1;
  FFT_SCALAR ***u_brick_a1;
  FFT_SCALAR ***v0_brick_a1,***v1_brick_a1,***v2_brick_a1,***v3_brick_a1,***v4_brick_a1,***v5_brick_a1;

  FFT_SCALAR ***density_brick_a2;
  FFT_SCALAR ***vdx_brick_a2,***vdy_brick_a2,***vdz_brick_a2;
  FFT_SCALAR *density_fft_a2;
  FFT_SCALAR ***u_brick_a2;
  FFT_SCALAR ***v0_brick_a2,***v1_brick_a2,***v2_brick_a2,***v3_brick_a2,***v4_brick_a2,***v5_brick_a2;

  FFT_SCALAR ***density_brick_a3;
  FFT_SCALAR ***vdx_brick_a3,***vdy_brick_a3,***vdz_brick_a3;
  FFT_SCALAR *density_fft_a3;
  FFT_SCALAR ***u_brick_a3;
  FFT_SCALAR ***v0_brick_a3,***v1_brick_a3,***v2_brick_a3,***v3_brick_a3,***v4_brick_a3,***v5_brick_a3;

  FFT_SCALAR ***density_brick_a4;
  FFT_SCALAR ***vdx_brick_a4,***vdy_brick_a4,***vdz_brick_a4;
  FFT_SCALAR *density_fft_a4;
  FFT_SCALAR ***u_brick_a4;
  FFT_SCALAR ***v0_brick_a4,***v1_brick_a4,***v2_brick_a4,***v3_brick_a4,***v4_brick_a4,***v5_brick_a4;

  FFT_SCALAR ***density_brick_a5;
  FFT_SCALAR ***vdx_brick_a5,***vdy_brick_a5,***vdz_brick_a5;
  FFT_SCALAR *density_fft_a5;
  FFT_SCALAR ***u_brick_a5;
  FFT_SCALAR ***v0_brick_a5,***v1_brick_a5,***v2_brick_a5,***v3_brick_a5,***v4_brick_a5,***v5_brick_a5;

  FFT_SCALAR ***density_brick_a6;
  FFT_SCALAR ***vdx_brick_a6,***vdy_brick_a6,***vdz_brick_a6;
  FFT_SCALAR *density_fft_a6;
  FFT_SCALAR ***u_brick_a6;
  FFT_SCALAR ***v0_brick_a6,***v1_brick_a6,***v2_brick_a6,***v3_brick_a6,***v4_brick_a6,***v5_brick_a6;

  //// needed for each interaction type
  double *greensfn;
  double **vg;
  double **vg2;

  double *greensfn_6;
  double **vg_6;
  double **vg2_6;

  double *fkx,*fky,*fkz;
  double *fkx2, *fky2, *fkz2;
  double *fkx_6, *fky_6, *fkz_6;
  double *fkx2_6, *fky2_6, *fkz2_6;
  double *gf_b;
  double *gf_b_6;

  double *sf_precoeff1, *sf_precoeff2, *sf_precoeff3, *sf_precoeff4, 
    *sf_precoeff5, *sf_precoeff6;
  double *sf_precoeff1_6, *sf_precoeff2_6, *sf_precoeff3_6, 
    *sf_precoeff4_6, *sf_precoeff5_6, *sf_precoeff6_6;
  FFT_SCALAR **rho1d,**rho_coeff;
  FFT_SCALAR **drho1d, **drho_coeff;
  FFT_SCALAR **rho1d_6, **rho_coeff_6;
  FFT_SCALAR **drho1d_6, **drho_coeff_6;
  FFT_SCALAR *work1,*work2;
  FFT_SCALAR *work1_6, *work2_6;


  class FFT3d *fft1,*fft2 ;
  class FFT3d *fft1_6, *fft2_6;
  class Remap *remap;
  class Remap *remap_6;
  class CommGrid *cg;
  class CommGrid *cg_peratom;
  class CommGrid *cg_6;
  class CommGrid *cg_peratom_6;

  int **part2grid;             // storage for particle -> grid mapping
  int **part2grid_6;
  int nmax;

  int triclinic;               // domain settings, orthog or triclinic
  double *boxlo;
                               // TIP4P settings
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  double qdist;                // distance from O site to negative charge
  double alpha;                // geometric factor

  void init_coeffs();	

  void set_grid();
  void set_grid_6();
  void set_init_g6();
  void set_fft_parameters(int&, int&, int&, int&, int&,int&,
                          int&, int&,int&, int&, int&,int&,
                          int&, int&,int&, int&, int&,int&,
                          int&, int&,int&, int&, int&,
			  int&, int&, int&,
		          double&, double&, int&);
  void set_n_pppm_6();
  void adjust_gewald();
  void adjust_gewald_6();
  double f();
  double derivf();
  double f_6();
  double derivf_6();
  double final_accuracy();
  double final_accuracy_6();
  double lj_rspace_error();
  double compute_qopt();
  double compute_qopt_6();
  double compute_qopt_ik();
  double compute_qopt_ad();
  double compute_qopt_6_ik();
  double compute_qopt_6_ad();

  void calc_csum();

  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  double rms(double, double, bigint, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom(double*, int);
  double gf_denom(double, double, double, double*, int);
  

  void compute_sf_precoeff(int, int, int, int, 
                           int, int, int,
                           int, int, int,
                           double*, double*, double*,
                           double*, double*, double*);
  void compute_gf();
  void compute_sf_coeff();
  void compute_gf_6();
  void compute_sf_coeff_6();


  virtual void particle_map(double, double, double,
                             double, int **, int, int,
                             int, int, int,
                             int, int, int);
  virtual void particle_map_c(double, double, double,
			      double, int **, int, int,
                              int, int, int,
                              int, int, int );
  virtual void make_rho_c();
  virtual void make_rho_g();
  virtual void make_rho_a();

  virtual void brick2fft(int, int, int, int, int, int,
			 FFT_SCALAR ***, FFT_SCALAR *, FFT_SCALAR *,
                         LAMMPS_NS::Remap *);
  virtual void brick2fft_a();

  virtual void poisson_ik(FFT_SCALAR *, FFT_SCALAR *,
		          FFT_SCALAR *, LAMMPS_NS::FFT3d *,LAMMPS_NS::FFT3d *, 
                          int, int, int, int, int, int, int,
		          int, int, int, int, int, int,
                          int, int, int, double&, double *,
                          double *, double *, double *,
                          double *, double *, double *,
		          FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, double *, double **, double **,
                          FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                          FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_ad(FFT_SCALAR*, FFT_SCALAR*,
                          FFT_SCALAR*, LAMMPS_NS::FFT3d*,LAMMPS_NS::FFT3d*, 
                          int, int, int, int,
                          int, int, int, int, int, int,
                          int, int, int, int, int, int,
                          double&, double*,
                          double*, double**, double**,
                          FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***,
                          FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***);

  virtual void poisson_peratom(FFT_SCALAR*, FFT_SCALAR*, LAMMPS_NS::FFT3d*, 
                               double**, double**, int,
                               int, int, int, int, int, int,
                               FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***,
                               FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***);
  virtual void poisson_2s_ik(FFT_SCALAR *, FFT_SCALAR *,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
			     FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
			     FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);
  virtual void poisson_2s_ad(FFT_SCALAR *, FFT_SCALAR *,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
			     FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
			     FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_2s_peratom(FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***,
				  FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***,
                                  FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***,
                                  FFT_SCALAR***, FFT_SCALAR***, FFT_SCALAR***);


  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  virtual void fieldforce_c_peratom();
  virtual void fieldforce_g_ik();
  virtual void fieldforce_g_ad();
  virtual void fieldforce_g_peratom();
  virtual void fieldforce_a_ik();
  virtual void fieldforce_a_ad();
  virtual void fieldforce_a_peratom();
  void procs2grid2d(int,int,int,int *, int*);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &, 
		     const FFT_SCALAR &, int, FFT_SCALAR **, FFT_SCALAR **);
  void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &, 
		      const FFT_SCALAR &, int, FFT_SCALAR **, FFT_SCALAR **);
  void compute_rho_coeff(FFT_SCALAR **,FFT_SCALAR **, int);
  void slabcorr(int);

  // grid communication

  void pack_forward(int, FFT_SCALAR *, int, int *);
  void unpack_forward(int, FFT_SCALAR *, int, int *);
  void pack_reverse(int, FFT_SCALAR *, int, int *);
  void unpack_reverse(int, FFT_SCALAR *, int, int *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot (yet) use PPPMDisp with triclinic box

UNDOCUMENTED

E: Cannot use PPPMDisp with 2d simulation

UNDOCUMENTED

E: Cannot use nonperiodic boundaries with PPPMDisp

UNDOCUMENTED

E: Incorrect boundaries with slab PPPMDisp

UNDOCUMENTED

E: PPPMDisp coulomb order cannot be greater than %d

UNDOCUMENTED

E: KSpace style is incompatible with Pair style

Setting a kspace style requires that a pair style with a long-range
Coulombic and Dispersion component be selected.

E: Unsupported mixing rule in kspace_style pppm/disp for pair_style %s

UNDOCUMENTED

E: Unsupported order in kspace_style pppm/disp pair_style %s

UNDOCUMENTED

W: Charges are set, but coulombic solver is not used

UNDOCUMENTED

E: Kspace style with selected options requires atom attribute q

The atom style defined does not have these attributes.
Change the atom style or switch of the coulomb solver.

E: Cannot use kspace solver with selected options on system with no charge

No atoms in system have a non-zero charge. Change charges or change 
options of the kspace solver/pair style.

W: System is not charge neutral, net charge = %g

The total charge on all atoms on the system is not 0.0, which
is not valid for Ewald or PPPM coulombic solvers.

E: Bond and angle potentials must be defined for TIP4P

Cannot use TIP4P pair potential unless bond and angle potentials
are defined.

E: Bad TIP4P angle type for PPPMDisp/TIP4P

UNDOCUMENTED

E: Bad TIP4P bond type for PPPMDisp/TIP4P

UNDOCUMENTED

W: Reducing PPPMDisp Coulomb order b/c stencil extends beyond neighbor processor.

UNDOCUMENTED

E: PPPMDisp Coulomb grid is too large

UNDOCUMENTED

E: Coulomb PPPMDisp order has been reduced below minorder

UNDOCUMENTED

W: Reducing PPPMDisp Dispersion order b/c stencil extends beyond neighbor processor

UNDOCUMENTED

E: PPPMDisp Dispersion grid is too large

UNDOCUMENTED

E: Dispersion PPPMDisp order has been reduced below minorder

UNDOCUMENTED

E: PPPM grid stencil extends beyond nearest neighbor processor

UNDOCUMENTED

E: epsilon or sigma reference not set by pair style in PPPMDisp

UNDOCUMENTED

E: KSpace accuracy too large to estimate G vector

UNDOCUMENTED

E: Could not compute grid size for Coulomb interaction!

UNDOCUMENTED

E: Could not compute g_ewald

UNDOCUMENTED

E: Could not adjust g_ewald_6

UNDOCUMENTED

E: Cannot compute initial g_ewald_disp

LAMMPS failed to compute an initial guess for the PPPM_disp g_ewald_6
factor that partitions the computation between real space and k-space
for Disptersion interactions.

E: Could not compute grid size for Dispersion!

UNDOCUMENTED

E: Out of range atoms - cannot compute PPPMDisp

UNDOCUMENTED

U: Cannot (yet) use PPPM_disp with triclinic box

This feature is not yet supported.

U: Cannot use PPPM_disp with 2d simulation

The kspace style pppm_disp cannot be used in 2d simulations.  You can use
2d PPPM_disp in a 3d simulation; see the kspace_modify command.

U: Cannot use nonperiodic boundaries with PPPM_disp

For kspace style pppm_disp, all 3 dimensions must have periodic boundaries
unless you use the kspace_modify command to define a 2d slab with a
non-periodic z dimension.

U: Incorrect boundaries with slab PPPM_disp

Must have periodic x,y dimensions and non-periodic z dimension to use
2d slab option with PPPM_disp.

U: PPPM_disp coulomb order cannot be greater than %d

Self-explanatory.

U: PPPM_disp dispersion order cannot be greater than %d

Self-explanatory.

U: Unsupported mixing rule in kspace_style pppm_disp for pair_style %s

PPPM_disp requires arithemtic or geometric mixing rules.

U: Unsupported order in kspace_style pppm_disp pair_style %s

PPPM_disp only works for 1/r and 1/r^6 potentials

U: Charges are set, but coulombic long-range solver is not used.

Charges have been specified, however, calculations are performed
as if they were zero.

U: Bad TIP4P angle type for PPPM_disp/TIP4P

Specified angle type is not valid.

U: Bad TIP4P bond type for PPPM_disp/TIP4P

Specified bond type is not valid.

U: Reducing PPPM_disp Coulomb order b/c stencil extends beyond neighbor processor

LAMMPS is attempting this in order to allow the simulation
to run.  It should not effect the PPPM_disp accuracy.

U: Reducing PPPM_disp Dispersion order b/c stencil extends beyond neighbor processor

LAMMPS is attempting this in order to allow the simulation
to run.  It should not effect the PPPM_disp accuracy.

U: PPPM_disp Coulomb grid is too large

The global PPPM_disp grid for Coulomb interactions is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 16384.  You likely need to decrease the
requested precision.

U: PPPM_grid Dispersion grid is too large

One of the PPPM_disp grids for Dispersion interactions is larger than OFFSET in one or more dimensions.
OFFSET is currently set to 16384.  You likely need to decrease the
requested precision.

U: Coulomb PPPM_disp order has been reduced to 0

LAMMPS has attempted to reduce the PPPM_disp coulomb order to enable the simulation
to run, but can reduce the order no further.  Try increasing the
accuracy of PPPM_disp coulomb by reducing the tolerance size, thus inducing a 
larger PPPM_disp coulomb grid.

U: Dispersion PPPM_disp order has been reduced to 0

LAMMPS has attempted to reduce the PPPM_disp dispersion order to enable the simulation
to run, but can reduce the order no further.  Try increasing the
accuracy of PPPM_disp dispersion by reducing the tolerance size, thus inducing a 
larger PPPM_disp dispersion grid.

U: Cannot compute PPPM_disp g_ewald

LAMMPS failed to compute a valid approximation for the PPPM_disp g_ewald
factor that partitions the computation between real space and k-space
for Coulomb interactions.

U: Cannot compute final g_ewald_disp

LAMMPS failed to compute a final value for the PPPM_disp g_ewald_6
factor that partitions the computation between real space and k-space
for Disptersion interactions.

U: Out of range atoms - cannot compute PPPM_disp

One or more atoms are attempting to map their charge to a PPPM_disp grid
point that is not owned by a processor.  This is likely for one of two
reasons, both of them bad.  First, it may mean that an atom near the
boundary of a processor's sub-domain has moved more than 1/2 the
"neighbor skin distance"_neighbor.html without neighbor lists being
rebuilt and atoms being migrated to new processors.  This also means
you may be missing pairwise interactions that need to be computed.
The solution is to change the re-neighboring criteria via the
"neigh_modify"_neigh_modify command.  The safest settings are "delay 0
every 1 check yes".  Second, it may mean that an atom has moved far
outside a processor's sub-domain or even the entire simulation box.
This indicates bad physics, e.g. due to highly overlapping atoms, too
large a timestep, etc.

*/
