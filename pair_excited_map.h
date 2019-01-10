/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
//Written by Steven E Strong, based on pair_coul_cut.cpp
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(excited/map,PairExcitedMap)

#else

#ifndef LMP_PAIR_EXCITED_MAP_H
#define LMP_PAIR_EXCITED_MAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairExcitedMap : public Pair {
 public:
  PairExcitedMap(class LAMMPS *);
  virtual ~PairExcitedMap() {};
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void *extract(const char *, int &);

 protected:
  double cut_global,cutsq;
  tagint tagO,tagH,tagH0;       //ids of excited chromophore; H0=non-excited H
  double mapA,mapB;             //coefficients of map energy (a*E + b*E^2)
};

}

#endif
#endif
