/* -*- c++ -*- ----------------------------------------------------------
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

KSpaceStyle(pppm/old,PPPMOld)

#else

#ifndef LMP_PPPM_OLD_H
#define LMP_PPPM_OLD_H

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

class PPPMOld : public KSpace {
 public:
  PPPMOld(class LAMMPS *, int, char **);
  virtual ~PPPMOld();
  virtual void init();
  virtual void setup();
  virtual void compute(int, int);
  virtual int timing_1d(int, double &);
  virtual int timing_3d(int, double &);
  virtual double memory_usage();

  virtual void compute_group_group(int, int, int);

 protected:
  int me,nprocs;
  int nfactors;
  int *factors;
  double qsum,qsqsum;
  double cutoff;
  double volume;
  double delxinv,delyinv,delzinv,delvolinv;
  double shift,shiftone;
  int peratom_allocate_flag;

  int nxlo_in,nylo_in,nzlo_in,nxhi_in,nyhi_in,nzhi_in;
  int nxlo_out,nylo_out,nzlo_out,nxhi_out,nyhi_out,nzhi_out;
  int nxlo_ghost,nxhi_ghost,nylo_ghost,nyhi_ghost,nzlo_ghost,nzhi_ghost;
  int nxlo_fft,nylo_fft,nzlo_fft,nxhi_fft,nyhi_fft,nzhi_fft;
  int nlower,nupper;
  int ngrid,nfft,nfft_both;
  int nbuf,nbuf_peratom;

  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick,***vdy_brick,***vdz_brick;
  FFT_SCALAR ***u_brick;
  FFT_SCALAR ***v0_brick,***v1_brick,***v2_brick;
  FFT_SCALAR ***v3_brick,***v4_brick,***v5_brick;
  double *greensfn;
  double **vg;
  double *fkx,*fky,*fkz;
  FFT_SCALAR *density_fft;
  FFT_SCALAR *work1,*work2;
  FFT_SCALAR *buf1,*buf2,*buf3,*buf4;

  double *gf_b;
  FFT_SCALAR **rho1d,**rho_coeff;

  // group-group interactions

  int group_allocate_flag;
  FFT_SCALAR ***density_A_brick,***density_B_brick;
  FFT_SCALAR *density_A_fft,*density_B_fft;


  class FFT3d *fft1,*fft2;
  class Remap *remap;

  int **part2grid;             // storage for particle -> grid mapping
  int nmax;

  int triclinic;               // domain settings, orthog or triclinic
  double *boxlo;
                               // TIP4P settings
  int typeH,typeO;             // atom types of TIP4P water H and O atoms
  double qdist;                // distance from O site to negative charge
  double alpha;                // geometric factor

  void set_grid();
  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  double rms(double, double, bigint, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom();

  virtual void particle_map();
  virtual void make_rho();
  virtual void brick2fft();
  virtual void fillbrick();
  virtual void fillbrick_peratom();
  virtual void poisson();
  virtual void poisson_peratom();
  virtual void fieldforce();
  virtual void fieldforce_peratom();
  void procs2grid2d(int,int,int,int *, int*);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &,
                     const FFT_SCALAR &);
  void compute_rho_coeff();
  void slabcorr();

  // group-group interactions

  virtual void allocate_groups();
  virtual void deallocate_groups();
  virtual void make_rho_groups(int, int, int);
  virtual void poisson_groups(int);

/* ----------------------------------------------------------------------
   denominator for Hockney-Eastwood Green's function
     of x,y,z = sin(kx*deltax/2), etc

            inf                 n-1
   S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
           j=-inf               l=0

          = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
   gf_b = denominator expansion coeffs
------------------------------------------------------------------------- */

  inline double gf_denom(const double &x, const double &y,
                         const double &z) const {
    double sx,sy,sz;
    sz = sy = sx = 0.0;
    for (int l = order-1; l >= 0; l--) {
      sx = gf_b[l] + sx*x;
      sy = gf_b[l] + sy*y;
      sz = gf_b[l] + sz*z;
    }
    double s = sx*sy*sz;
    return s*s;
  };
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
