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

#ifdef PAIR_CLASS

PairStyle(eam/fs/gpu,PairEAMFSGPU)

#else

#ifndef LMP_PAIR_EAM_FS_GPU_H
#define LMP_PAIR_EAM_FS_GPU_H

#include "pair_eam_gpu.h"

namespace LAMMPS_NS {

class PairEAMFSGPU : public PairEAMGPU {
public:
  PairEAMFSGPU(class LAMMPS *);
  virtual ~PairEAMFSGPU() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: No matching element in EAM potential file

The EAM potential file does not contain elements that match the
requested elements.

E: Cannot open EAM potential file %s

The specified EAM potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect element names in EAM potential file

The element names in the EAM file do not match those requested.

*/
