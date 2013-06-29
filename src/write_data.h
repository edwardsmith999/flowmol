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

#ifdef COMMAND_CLASS

CommandStyle(write_data,WriteData)

#else

#ifndef LMP_WRITE_DATA_H
#define LMP_WRITE_DATA_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class WriteData : protected Pointers {
 public:
  WriteData(class LAMMPS *);
  void command(int, char **);
  void write(char *);

 private:
  int me,nprocs;
  int pairflag;
  FILE *fp;
  bigint nbonds_local,nbonds;
  bigint nangles_local,nangles;

  void header();
  void type_arrays();
  void force_fields();
  void atoms();
  void velocities();
  void bonds();
  void angles();
  void dihedrals();
  void impropers();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Write_data command before simulation box is defined

UNDOCUMENTED

E: Illegal ... command

UNDOCUMENTED

E: Atom count is inconsistent, cannot write data file

UNDOCUMENTED

E: Cannot open data file %s

UNDOCUMENTED

*/
