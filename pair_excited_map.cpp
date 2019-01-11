/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
//Written by Steven E Strong, based on pair_coul_cut.cpp

//input file syntax is:
//pair_style excited/map cutoff idO mapA mapB
//no pair_coeff command should be used
//water molecules must be in order OHH OHH ...
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_excited_map.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairExcitedMap::PairExcitedMap(LAMMPS *lmp) : Pair(lmp) {
  single_enable=0;
  restartinfo  =0;
  reinitflag   =0;
  one_coeff    =1;

  ewaldflag = pppmflag = msmflag = dipoleflag = 1;
  if (!force->newton)
    error->all(FLERR,"Newton's 3rd law must be enabled to add forces to ghost atoms");
}

PairExcitedMap::~PairExcitedMap() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

void PairExcitedMap::compute(int eflag, int vflag)
{
  int j,jj,jnum,itype,jtype,ih;
  double qtmp,delx,dely,delz,epair,fpair;
  double rsq,r2inv,rinv,factor_coul;
  int *jlist,*numneigh,**firstneigh;

  epair = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal+nghost;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  //get ids of excited molecule on this proc
  int idO =atom->map(tagO);
  int idH =atom->map(tagH);
  int idH0=atom->map(tagH0);
  int typeO=type[idO];
  itype = type[idH];

  //because neighbor list includes everything within cutoff,
  //don't need to worry about communicating E field

  //TODO: this assumes that idO,idH,and idH0 are always on this proc (id>-1)
  //might be too expensive to test this every time
  //should verify that it is ensured
  if (idO==-1)
    return;  //no relevant atoms on this proc, only in a very big simulation

  double *xO=x[idO];
  double *xH=x[idH];
  double oh[3];
  oh[0] = xH[0] - xO[0];
  oh[1] = xH[1] - xO[1];
  oh[2] = xH[2] - xO[2];
  double rOH=oh[0]*oh[0] + oh[1]*oh[1] + oh[2]*oh[2];
  rOH = sqrt(rOH);
  double rOHinv= 1.0/rOH;
  oh[0] *= rOHinv; //normalize
  oh[1] *= rOHinv;
  oh[2] *= rOHinv;

  //neighbor list info
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  jlist = firstneigh[idH];
  jnum = numneigh[idH];

  //compute electric field at H, and calculate derivative components
  double eHvec[3] = {0.0,0.0,0.0};
  double fO[3]    = {0.0,0.0,0.0};
  double fH[3]    = {0.0,0.0,0.0};
  double uhj[3];   //unit vector in j->H direction
  //TODO: could compute one of these from the sum of the others
  double tmpdot,qfact,eHtmp;
  double **fI;
  //TODO: can this proc add forces to ghost atoms?
  memory->create(fI,ntotal,3,"pairExcitedMap:forces");
  size_t nbytes = sizeof(double)*ntotal*3;
  //TODO: not necessary, but could help to identify atoms with force
  if (nbytes) memset(&fI[0][0],0.0,nbytes);

  //loop through neighs of H
  for (jj = 0; jj < jnum; jj++) {
    j = jlist[jj];

    //exclude same molecule from contributing to eH
    factor_coul = special_coul[sbmask(j)];
    if (factor_coul==0.0)
      continue;

    //filters out bits that encode special neighbors
    j &= NEIGHMASK;

    if (type[j]!=typeO)  //only test cutoff wrt O
      continue;

    delx = xH[0] - x[j][0];
    dely = xH[1] - x[j][1];
    delz = xH[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];  //TODO: don't need to care about type?

    if (rsq < cut2) {
      qtmp=q[j];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);

      //normalize hj vector
      uhj[0]=delx*rinv;
      uhj[1]=dely*rinv;
      uhj[2]=delz*rinv;

      eHtmp = qqrd2e * qtmp * r2inv;

      //accumulate eH
      eHvec[0] += uhj[0]*eHtmp;  //has lammps units of force/charge
      eHvec[1] += uhj[1]*eHtmp;
      eHvec[2] += uhj[2]*eHtmp;

      //TODO: do we need any unit conversions (qqrd2e, etc) here
      //get force vectors
      qfact  = qtmp * r2inv;
      tmpdot = uhj[0]*oh[0] + uhj[1]*oh[1] + uhj[2]*oh[2];
      fO[0] += qfact * ( oh[0]*tmpdot - uhj[0] ) * rOHinv;
      fO[1] += qfact * ( oh[1]*tmpdot - uhj[1] ) * rOHinv;
      fO[2] += qfact * ( oh[2]*tmpdot - uhj[2] ) * rOHinv;

      fI[j][0] = qfact * rinv * (3*uhj[0]*tmpdot - oh[0]);
      fI[j][1] = qfact * rinv * (3*uhj[1]*tmpdot - oh[1]);
      fI[j][2] = qfact * rinv * (3*uhj[2]*tmpdot - oh[2]);

      fH[0] +=  qfact * ( oh[0]*rinv + uhj[0]*rOHinv -
			  (oh[0]*rOHinv + 3*uhj[0]*rinv)*tmpdot );
      fH[1] +=  qfact * ( oh[1]*rinv + uhj[1]*rOHinv -
			  (oh[1]*rOHinv + 3*uhj[1]*rinv)*tmpdot );
      fH[2] +=  qfact * ( oh[2]*rinv + uhj[2]*rOHinv -
			  (oh[2]*rOHinv + 3*uhj[2]*rinv)*tmpdot );

      //now compute field from hydrogens on that oxygen
      //TODO: probably faster to unroll this loop
      for (ih=1; ih<3; ih++) { //TODO: are locals IDs in order OHH too?
	qtmp = q[j];
	delx = xH[0] - x[j+ih][0];
	dely = xH[1] - x[j+ih][1];
	delz = xH[2] - x[j+ih][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];

	r2inv = 1.0/rsq;
	rinv = sqrt(r2inv);

	//normalize hj vector
	uhj[0]=delx*rinv;
	uhj[1]=dely*rinv;
	uhj[2]=delz*rinv;

	eHtmp = qqrd2e * qtmp * r2inv;

	//accumulate eH
	eHvec[0] += uhj[0]*eHtmp;  //has lammps units of force/charge
	eHvec[1] += uhj[1]*eHtmp;
	eHvec[2] += uhj[2]*eHtmp;

	//get force vectors
	//TODO: it seems like there is an extra factor of L^-3 here
	//maybe this is accounted for in mapA and mapB
	qfact  = qtmp * r2inv;
	tmpdot = uhj[0]*oh[0] + uhj[1]*oh[1] + uhj[2]*oh[2];
	fO[0] += qfact * ( oh[0]*tmpdot - uhj[0] ) * rOHinv;
	fO[1] += qfact * ( oh[1]*tmpdot - uhj[1] ) * rOHinv;
	fO[2] += qfact * ( oh[2]*tmpdot - uhj[2] ) * rOHinv;

	fI[j+ih][0] = qfact * rinv * (3*uhj[0]*tmpdot - oh[0]);
	fI[j+ih][1] = qfact * rinv * (3*uhj[1]*tmpdot - oh[1]);
	fI[j+ih][2] = qfact * rinv * (3*uhj[2]*tmpdot - oh[2]);

	fH[0] +=  qfact * ( oh[0]*rinv + uhj[0]*rOHinv -
			    (oh[0]*rOHinv + 3*uhj[0]*rinv)*tmpdot );
	fH[1] +=  qfact * ( oh[1]*rinv + uhj[1]*rOHinv -
			    (oh[1]*rOHinv + 3*uhj[1]*rinv)*tmpdot );
	fH[2] +=  qfact * ( oh[2]*rinv + uhj[2]*rOHinv -
			    (oh[2]*rOHinv + 3*uhj[2]*rinv)*tmpdot );
      }
    }
  }
  double eH = eHvec[0]*oh[0] + eHvec[1]*oh[1] + eHvec[2]*oh[2];
  double mapBE = mapB * eH;

  //consolidate fO and fH into fI
  for (int ii=0; ii<3; ii++) {
    if ( fI[idO][ii] || fI[idH][ii] || fI[idH0][ii] )
      error->one(FLERR,"forces are being added to the excited molecule");
    fI[idO][ii]=fO[ii];
    fI[idH][ii]=fH[ii];
  }

  // loop over neighbors of excited H and compute force due to eH
  // this should do all neighbors, but not itself
  for (jj = 0; jj < jnum; jj++) {
    j = jlist[jj];

    //factor_coul = special_coul[sbmask(j)];
    j &= NEIGHMASK;

    f[j][0] += mapA*fI[j][0] + mapBE*fI[j][0];
    f[j][1] += mapA*fI[j][1] + mapBE*fI[j][1];
    f[j][2] += mapA*fI[j][2] + mapBE*fI[j][2];


  }

  //add force on excited H
  j=idH;
  f[j][0] += fI[j][0];
  f[j][1] += fI[j][0];
  f[j][2] += fI[j][0];

  if (eflag)
    epair = - mapA*eH - mapBE*eH/2;  //TODO: energy units?
  //TODO: need to figure out virial if want to use pressure
  //if (evflag) ev_tally(i,j,nlocal,newton_pair,
  //			 0.0,epair,fpair,delx,dely,delz);

  memory->destroy(fI);

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
//read when pair_style command is issued
------------------------------------------------------------------------- */

void PairExcitedMap::settings(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  cut2 = cut_global*cut_global;

  tagO = force->inumeric(FLERR,arg[1]);
  tagH = tagO+1; //this requires that hydrogen atoms always follow O
  tagH0= tagH+1;

  //these should be input in lammps units of energy/efield (charge * length)
  //we incorporate the negative sign and constants here
  mapA = - force->numeric(FLERR,arg[2]);
  mapB = -2 * force->numeric(FLERR,arg[3]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairExcitedMap::coeff(int narg, char **arg)
{
  if (narg != 2)
    error->all(FLERR,"Incorrect args for pair coefficients");

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  //allocate arrays that pair class expects
  if (!allocated) allocate();
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairExcitedMap::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style excited/map requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairExcitedMap::init_one(int i, int j)
{
  return cut_global;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairExcitedMap::write_restart(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&cut2,sizeof(double),1,fp);
  fwrite(&mapA,sizeof(double),1,fp);
  fwrite(&mapB,sizeof(double),1,fp);
  fwrite(&tagO,sizeof(tagint),1,fp);
  fwrite(&tagH,sizeof(tagint),1,fp);
  fwrite(&tagH0,sizeof(tagint),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairExcitedMap::read_restart(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&cut2,sizeof(double),1,fp);
    fread(&mapA,sizeof(double),1,fp);
    fread(&mapB,sizeof(double),1,fp);

    fread(&tagO,sizeof(tagint),1,fp);
    fread(&tagH,sizeof(tagint),1,fp);
    fread(&tagH0,sizeof(tagint),1,fp);
  }

  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mapA,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mapB,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&tagO,1,MPI_LMP_TAGINT,0,world);
  MPI_Bcast(&tagH,1,MPI_LMP_TAGINT,0,world);
  MPI_Bcast(&tagH0,1,MPI_LMP_TAGINT,0,world);

  allocate();
}

void *PairExcitedMap::extract(const char *str, int &dim) { return NULL; }

void PairExcitedMap::allocate() {
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 1;
      cutsq[i][j] = cut2;
    }

  allocated=1;
}
