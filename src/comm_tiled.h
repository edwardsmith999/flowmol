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

#ifndef LMP_COMM_TILED_H
#define LMP_COMM_TILED_H

#include "comm.h"

namespace LAMMPS_NS {

class CommTiled : public Comm {
 public:
  CommTiled(class LAMMPS *);
  CommTiled(class LAMMPS *, class Comm *);
  virtual ~CommTiled();

  void init();
  void setup();                        // setup comm pattern
  void forward_comm(int dummy = 0);    // forward comm of atom coords
  void reverse_comm();                 // reverse comm of forces
  void exchange();                     // move atoms to new procs
  void borders();                      // setup list of atoms to comm

  void forward_comm_pair(class Pair *);    // forward comm from a Pair
  void reverse_comm_pair(class Pair *);    // reverse comm from a Pair
  void forward_comm_fix(class Fix *);      // forward comm from a Fix
  void reverse_comm_fix(class Fix *);      // reverse comm from a Fix
  void forward_comm_variable_fix(class Fix *); // variable-size variant
  void reverse_comm_variable_fix(class Fix *); // variable-size variant
  void forward_comm_compute(class Compute *);  // forward from a Compute
  void reverse_comm_compute(class Compute *);  // reverse from a Compute
  void forward_comm_dump(class Dump *);    // forward comm from a Dump
  void reverse_comm_dump(class Dump *);    // reverse comm from a Dump

  void forward_comm_array(int, double **);          // forward comm of array
  int exchange_variable(int, double *, double *&);  // exchange on neigh stencil
  bigint memory_usage();

 private:
  int triclinic;                    // 0 if domain is orthog, 1 if triclinic
  int map_style;                    // non-0 if global->local mapping is done
  int size_forward;                 // # of per-atom datums in forward comm
  int size_reverse;                 // # of datums in reverse comm
  int size_border;                  // # of datums in forward border comm

  int nswap;                    // # of swaps to perform = 2*dim
  int *nsendproc,*nrecvproc;    // # of procs to send/recv to/from in each swap
  int *sendother;               // 1 if send to any other proc in each swap
  int *sendself;                // 1 if send to self in each swap
  int *nprocmax;                // current max # of send procs for each swap
  int **sendproc,**recvproc;    // proc to send/recv to/from per swap/proc
  int **sendnum,**recvnum;      // # of atoms to send/recv per swap/proc
  int **size_forward_recv;      // # of values to recv in each forward swap/proc
  int **firstrecv;              // where to put 1st recv atom per swap/proc
  int **size_reverse_send;      // # to send in each reverse comm per swap/proc
  int **size_reverse_recv;      // # to recv in each reverse comm per swap/proc
  int **forward_recv_offset;  // forward comm offsets in buf_recv per swap/proc
  int **reverse_recv_offset;  // reverse comm offsets in buf_recv per swap/proc

  int ***sendlist;              // list of atoms to send per swap/proc
  int **maxsendlist;            // max size of send list per swap/proc
  int **pbc_flag;               // general flag for sending atoms thru PBC
  int ***pbc;                   // dimension flags for PBC adjustments

  double ***sendbox;            // bounding box of atoms to send per swap/proc

  double *buf_send;             // send buffer for all comm
  double *buf_recv;             // recv buffer for all comm
  int maxsend,maxrecv;          // current size of send/recv buffer
  int maxforward,maxreverse;    // max # of datums in forward/reverse comm

  int maxexchange;              // max # of datums/atom in exchange comm 
  int bufextra;                 // extra space beyond maxsend in send buffer

  int maxreqstat;               // max size of Request and Status vectors
  MPI_Request *requests;
  MPI_Status *statuses;

  int comm_x_only,comm_f_only;      // 1 if only exchange x,f in for/rev comm

  struct RCBinfo {
    double mysplit[3][2];      // fractional RCB bounding box for one proc
    double cut;        	       // position of cut this proc owns
    int dim;	               // dimension = 0/1/2 of cut
  };

  int *overlap;
  RCBinfo *rcbinfo;       // list of RCB info for all procs

  void init_buffers();

  // box drop and other functions

  typedef void (CommTiled::*BoxDropPtr)(int, double *, double *, int &, int &);
  BoxDropPtr box_drop;
  void box_drop_brick(int, double *, double *, int &, int &);
  void box_drop_tiled(int, double *, double *, int &, int &);
  void box_drop_tiled_recurse(double *, double *, int, int, int &, int &);

  typedef void (CommTiled::*BoxOtherPtr)(int, int, int, double *, double *);
  BoxOtherPtr box_other;
  void box_other_brick(int, int, int, double *, double *);
  void box_other_tiled(int, int, int, double *, double *);

  void grow_send(int, int);            // reallocate send buffer
  void grow_recv(int);                 // free/allocate recv buffer
  void grow_list(int, int, int);       // reallocate sendlist for one swap/proc
  void allocate_swap(int);             // allocate swap arrays
  void grow_swap_send(int, int, int);  // grow swap arrays for send and recv
  void grow_swap_recv(int, int);
  void deallocate_swap(int);           // deallocate swap arrays
};

}

#endif

/* ERROR/WARNING messages:

*/
