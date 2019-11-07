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
//pair_style excited/map cutoff typeO idO mapA mapB
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
#include "domain.h"
#include "utils.h"

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

  memory->create(fI,neighbor->oneatom,3,"pairExcitedMap:forces");
  memory->create(indF,neighbor->oneatom,"pairExcitedMap:index");
}

PairExcitedMap::~PairExcitedMap() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
  memory->destroy(fI);
  memory->destroy(indF);
}

void PairExcitedMap::compute(int eflag, int vflag)
{
  int j,jj,ih,iih;
  double qtmp,delx,dely,delz;
  double rsq,r2inv,rinv,factor_coul;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal+nghost;

  //get ids of excited molecule on this proc
  int idH =atom->map(tagH);
  int idO =atom->map(tagO);
  idO = domain->closest_image(idH,idO);
  int idH0=atom->map(tagH0); //TODO: prob don't need to track idH0
  idH0 = domain->closest_image(idO,idH0);
  if (type[idO]!=typeO)
    error->one(FLERR,"excited O atom has the wrong type");

  //because neighbor list includes everything within cutoff,
  //don't need to worry about communicating E field
  //only compute forces if excited H is a local atom, not a ghost
  if (idH >= nlocal)
    return;

  //TODO: this assumes that idO,idH,and idH0 are always on the same proc
  //should verify that it is ensured
  if (idH==-1)
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

  //compute electric field at H, and calculate derivative components
  double eHvec[3] = {0.0,0.0,0.0};
  double fO[3]    = {0.0,0.0,0.0};
  double fH[3]    = {0.0,0.0,0.0};
  double uhj[3];   //unit vector in j->H direction
  double tmpdot,qfact,eHtmp;

  //loop through neighs of H
  int count=0;
  for (j=0; j<ntotal; j++) {

    if (type[j]!=typeO)  //only test cutoff wrt O
      continue;
    if (j==idO) //skip same molec
      continue;

    delx = xH[0] - x[j][0];
    dely = xH[1] - x[j][1];
    delz = xH[2] - x[j][2];
    //PBCs are accounted for in ghost communication
    rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < cut2) {
      qtmp=q[j];
      r2inv = 1.0/rsq;
      rinv = sqrt(r2inv);

      //normalize hj vector
      uhj[0]=delx*rinv;
      uhj[1]=dely*rinv;
      uhj[2]=delz*rinv;

      eHtmp = qtmp * r2inv;  //will convert to energy after loop

      //accumulate eH
      eHvec[0] += uhj[0]*eHtmp;
      eHvec[1] += uhj[1]*eHtmp;
      eHvec[2] += uhj[2]*eHtmp;

      //get force vectors
      qfact  = qtmp * r2inv;  //will convert after loop for efficiency
      tmpdot = uhj[0]*oh[0] + uhj[1]*oh[1] + uhj[2]*oh[2];
      fO[0] += qfact * ( oh[0]*tmpdot - uhj[0] ) * rOHinv; //Q/L^3
      fO[1] += qfact * ( oh[1]*tmpdot - uhj[1] ) * rOHinv;
      fO[2] += qfact * ( oh[2]*tmpdot - uhj[2] ) * rOHinv;

      fI[count][0] = qfact * rinv * (3*uhj[0]*tmpdot - oh[0]); //Q/L^3
      fI[count][1] = qfact * rinv * (3*uhj[1]*tmpdot - oh[1]);
      fI[count][2] = qfact * rinv * (3*uhj[2]*tmpdot - oh[2]);

      fH[0] +=  qfact * ( oh[0]*rinv + uhj[0]*rOHinv -
			  (oh[0]*rOHinv + 3*uhj[0]*rinv)*tmpdot ); //Q/L^3
      fH[1] +=  qfact * ( oh[1]*rinv + uhj[1]*rOHinv -
			  (oh[1]*rOHinv + 3*uhj[1]*rinv)*tmpdot );
      fH[2] +=  qfact * ( oh[2]*rinv + uhj[2]*rOHinv -
			  (oh[2]*rOHinv + 3*uhj[2]*rinv)*tmpdot );

      indF[count++]=j;

      //now compute field from hydrogens on that oxygen
      //TODO: probably faster to unroll this loop
      for (iih=1; iih<3; iih++) {
	ih = atom->map(tag[j] + iih);  //get local id of next H atom
	if (ih==-1)
          error->one(FLERR,"hydrogen is missing");
	ih = domain->closest_image(j,ih);

	qtmp = q[ih];
	delx = xH[0] - x[ih][0];
	dely = xH[1] - x[ih][1];
	delz = xH[2] - x[ih][2];
	rsq = delx*delx + dely*dely + delz*delz;
	r2inv = 1.0/rsq;
	rinv = sqrt(r2inv);

	//normalize hj vector
	uhj[0]=delx*rinv;
	uhj[1]=dely*rinv;
	uhj[2]=delz*rinv;

	eHtmp = qtmp * r2inv;  //will convert to energy after loop

	//accumulate eH
	eHvec[0] += uhj[0]*eHtmp;
	eHvec[1] += uhj[1]*eHtmp;
	eHvec[2] += uhj[2]*eHtmp;

	//get force vectors
	qfact  = qtmp * r2inv;  //will convert after loop for efficiency
	tmpdot = uhj[0]*oh[0] + uhj[1]*oh[1] + uhj[2]*oh[2];
	fO[0] += qfact * ( oh[0]*tmpdot - uhj[0] ) * rOHinv;   //Q/L^3
	fO[1] += qfact * ( oh[1]*tmpdot - uhj[1] ) * rOHinv;
	fO[2] += qfact * ( oh[2]*tmpdot - uhj[2] ) * rOHinv;

	fI[count][0] = qfact * rinv * (3*uhj[0]*tmpdot - oh[0]);  //Q/L^3
	fI[count][1] = qfact * rinv * (3*uhj[1]*tmpdot - oh[1]);
	fI[count][2] = qfact * rinv * (3*uhj[2]*tmpdot - oh[2]);

	fH[0] +=  qfact * ( oh[0]*rinv + uhj[0]*rOHinv -
			    (oh[0]*rOHinv + 3*uhj[0]*rinv)*tmpdot );  //Q/L^3
	fH[1] +=  qfact * ( oh[1]*rinv + uhj[1]*rOHinv -
			    (oh[1]*rOHinv + 3*uhj[1]*rinv)*tmpdot );
	fH[2] +=  qfact * ( oh[2]*rinv + uhj[2]*rOHinv -
			    (oh[2]*rOHinv + 3*uhj[2]*rinv)*tmpdot );

	indF[count++]=ih;
      }
    }
  }
  double eH = eHvec[0]*oh[0] + eHvec[1]*oh[1] + eHvec[2]*oh[2];
  eH *= force->qqrd2e; //energy/charge*length
  double mapC = mapA + mapB*eH;  //charge*length
  mapC *= force->qqrd2e; //convert Q/L^3 to E/QL^2 and to F

  //test that no forces already added to excited molecule
  for (int ii=0; ii<count; ii++)
    if (indF[ii]==idH || indF[ii]==idH0 || indF[ii]==idO )
      error->one(FLERR,"forces are being added to the excited molecule");

  //consolidate fO and fH into fI
  indF[count]=idH;
  indF[count+1]=idO;
  for (int ii=0; ii<3; ii++) {
    fI[count][ii]=fH[ii];
    fI[count+1][ii]=fO[ii];
  }
  count+=2;

  // loop over neighbors of excited H and compute force due to eH
  // this should do all neighbors, but not itself
  double fThis[3];
  double fTot[3] = {0.0,0.0,0.0};
  for (jj=0; jj<count; jj++) {
    j=indF[jj];

    //fI = force/charge*length
    fThis[0] = mapC*fI[jj][0];
    fThis[1] = mapC*fI[jj][1];
    fThis[2] = mapC*fI[jj][2];

    f[j][0] += fThis[0];
    f[j][1] += fThis[1];
    f[j][2] += fThis[2];

    fTot[0] += fThis[0];
    fTot[1] += fThis[1];
    fTot[2] += fThis[2];
  }

  //TODO: fTot should be zero
  if (fabs(fTot[0]) > 1e-10 ||
      fabs(fTot[1]) > 1e-10 ||
      fabs(fTot[2]) > 1e-10   ) {
    fprintf(screen,"%e %e %e\n",fTot[0],fTot[1],fTot[2]);
    error->one(FLERR,"total force is non-zero");
  }

  //this neglects the constant term in the map energy
  if (eflag_global)
    eng_vdwl +=  -mapA*eH - mapB*eH*eH/2;

  //TODO: is there a way to compute the contribution to the pressure
  //without fdotr?
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
//read when pair_style command is issued
------------------------------------------------------------------------- */

void PairExcitedMap::settings(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  cut2 = cut_global*cut_global;

  typeO= force->inumeric(FLERR,arg[1]);
  tagO = force->inumeric(FLERR,arg[2]);
  tagH = tagO+1; //this requires that hydrogen atoms always follow O
  tagH0= tagH+1;

  //these should be input in lammps units
  //we incorporate the negative sign and constants here
  mapA = - force->numeric(FLERR,arg[3]);  //charge*length
  mapB = -2 * force->numeric(FLERR,arg[4]);  //(charge*length)^2/energy
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
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut2,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&mapA,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&mapB,sizeof(double),1,fp,NULL,error);

    utils::sfread(FLERR,&tagO,sizeof(tagint),1,fp,NULL,error);
    utils::sfread(FLERR,&tagH,sizeof(tagint),1,fp,NULL,error);
    utils::sfread(FLERR,&tagH0,sizeof(tagint),1,fp,NULL,error);
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
