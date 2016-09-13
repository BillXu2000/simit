/*
  This is a Version 2.0 MPI + OpenMP implementation of LULESH

                 Copyright (c) 2010-2013.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 2.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

//////////////
DIFFERENCES BETWEEN THIS VERSION (2.x) AND EARLIER VERSIONS:
* Addition of regions to make work more representative of multi-material codes
* Default size of each domain is 30^3 (27000 elem) instead of 45^3. This is
  more representative of our actual working set sizes
* Single source distribution supports pure serial, pure OpenMP, MPI-only, 
  and MPI+OpenMP
* Addition of ability to visualize the mesh using VisIt 
  https://wci.llnl.gov/codes/visit/download.html
* Various command line options (see ./lulesh2.0 -h)
 -q              : quiet mode - suppress stdout
 -i <iterations> : number of cycles to run
 -s <size>       : length of cube mesh along side
 -r <numregions> : Number of distinct regions (def: 11)
 -b <balance>    : Load balance between regions of a domain (def: 1)
 -c <cost>       : Extra cost of more expensive regions (def: 1)
 -f <filepieces> : Number of file parts for viz output (def: np/9)
 -p              : Print out progress
 -v              : Output viz file (requires compiling with -DVIZ_MESH
 -h              : This message

 printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");

*Notable changes in LULESH 2.0

* Split functionality into different files
lulesh.cc - where most (all?) of the timed functionality lies
lulesh-comm.cc - MPI functionality
lulesh-init.cc - Setup code
lulesh-viz.cc  - Support for visualization option
lulesh-util.cc - Non-timed functions
*
* The concept of "regions" was added, although every region is the same ideal
*    gas material, and the same sedov blast wave problem is still the only
*    problem its hardcoded to solve.
* Regions allow two things important to making this proxy app more representative:
*   Four of the LULESH routines are now performed on a region-by-region basis,
*     making the memory access patterns non-unit stride
*   Artificial load imbalances can be easily introduced that could impact
*     parallelization strategies.  
* The load balance flag changes region assignment.  Region number is raised to
*   the power entered for assignment probability.  Most likely regions changes
*   with MPI process id.
* The cost flag raises the cost of ~45% of the regions to evaluate EOS by the
*   entered multiple. The cost of 5% is 10x the entered multiple.
* MPI and OpenMP were added, and coalesced into a single version of the source
*   that can support serial builds, MPI-only, OpenMP-only, and MPI+OpenMP
* Added support to write plot files using "poor mans parallel I/O" when linked
*   with the silo library, which in turn can be read by VisIt.
* Enabled variable timestep calculation by default (courant condition), which
*   results in an additional reduction.
* Default domain (mesh) size reduced from 45^3 to 30^3
* Command line options to allow numerous test cases without needing to recompile
* Performance optimizations and code cleanup beyond LULESH 1.0
* Added a "Figure of Merit" calculation (elements solved per microsecond) and
*   output in support of using LULESH 2.0 for the 2017 CORAL procurement
*
* Possible Differences in Final Release (other changes possible)
*
* High Level mesh structure to allow data structure transformations
* Different default parameters
* Minor code performance changes and cleanup

TODO in future versions
* Add reader for (truly) unstructured meshes, probably serial only
* CMake based build system

//////////////

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/

#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <unistd.h>

#include "lulesh.h"

// DLU : using simit namespace
#include "graph.h"
#include "program.h"
#include "mesh.h"
#include <cmath>
using namespace simit;


/******************************************/

int main(int argc, char *argv[])
{
  Domain *locDom ;
   Int_t numRanks ;
   Int_t myRank ;
   struct cmdLineOpts opts;

   numRanks = 1;
   myRank = 0;

   /* Set defaults that can be overridden by command line opts */
   opts.its = 9999999;
   opts.nx  = 30;
   opts.numReg = 11;
   opts.numFiles = (int)(numRanks+10)/9;
   opts.showProg = 0;
   opts.quiet = 0;
   opts.viz = 0;
   opts.balance = 1;
   opts.cost = 1;

   ParseCommandLineOptions(argc, argv, myRank, &opts);

   if ((myRank == 0) && (opts.quiet == 0)) {
      printf("Running problem size %d^3 per domain until completion\n", opts.nx);
      printf("Num processors: %d\n", numRanks);
      printf("Total number of elements: %lld\n\n", (long long int)(numRanks*opts.nx*opts.nx*opts.nx));
      printf("To run other sizes, use -s <integer>.\n");
      printf("To run a fixed number of iterations, use -i <integer>.\n");
      printf("Not Available in the SIMIT version : To run a more or less balanced region set, use -b <integer>.\n");
      printf("To change the relative costs of regions, use -c <integer>.\n");
      printf("To print out progress, use -p\n");
      printf("To write an output file for VisIt, use -v\n");
      printf("See help (-h) for more options\n\n");
   }

   // DLU - Initialize Simit
   simit::init("cpu", sizeof(double));

   // Set up the mesh and decompose. Assumes regular cubes for now
   Int_t col, row, plane, side;
   InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);

   // Build the main data structure and initialize it
   locDom = new Domain(numRanks, col, row, plane, opts.nx,
                       side, opts.numReg, opts.balance, opts.cost) ;

   // DLU - Create a graph and initialize it with Domain data
   Set nodes;
   Set elems(nodes, nodes, nodes, nodes, nodes, nodes, nodes, nodes);
   Set connects(elems,elems,elems,elems,elems,elems,elems);
   // The fields of the nodes set
   FieldRef<double,3> coord 	= nodes.addField<double,3>("coord");	// Nodal coordinates
   FieldRef<double,3> coord_local 	= nodes.addField<double,3>("coord_local");	// Nodal coordinates
   FieldRef<double,3> vel     	= nodes.addField<double,3>("vel");		// Nodal velocities
   FieldRef<double,3> a     	= nodes.addField<double,3>("a");		// Nodal accelerations
   FieldRef<double,3> f     	= nodes.addField<double,3>("f");		// Nodal forces
   FieldRef<double> nodalMass 	= nodes.addField<double>("nodalMass");	// Nodal mass
   FieldRef<int,3> symm     	= nodes.addField<int,3>("symm");		// Nodes on symmetry planes

   std::vector<ElementRef> nodeRefs;
   // initialize nodes values with Domain
   Index_t edgeElems = opts.nx;
   Index_t edgeNodes = edgeElems+1;  // In Lulesh same number of elems and nodes in each dimension
   Index_t nidx = 0 ;
   std::vector<int> symmX(edgeNodes*edgeNodes*edgeNodes,0) ;
   std::vector<int> symmY(edgeNodes*edgeNodes*edgeNodes,0) ;
   std::vector<int> symmZ(edgeNodes*edgeNodes*edgeNodes,0) ;
   for (int bc =0; bc<edgeNodes*edgeNodes; ++bc){
	   if (!locDom->symmXempty() != 0) {
		   symmX[locDom->symmX(bc)]=1;
	   }
	   if (!locDom->symmYempty() != 0) {
		   symmY[locDom->symmY(bc)]=1;
	   }
	   if (!locDom->symmZempty() != 0) {
		   symmZ[locDom->symmZ(bc)]=1;
	   }
   }
   for (Index_t plane=0; plane<edgeNodes; ++plane) {
     for (Index_t row=0; row<edgeNodes; ++row) {
       for (Index_t col=0; col<edgeNodes; ++col) {
    	     ElementRef node = nodes.add();
    	     nodeRefs.push_back(node);
    	     coord.set(node, {locDom->x(nidx),locDom->y(nidx),locDom->z(nidx)});
    	     vel.set(node, {locDom->xd(nidx),locDom->yd(nidx),locDom->zd(nidx)});
    	     a.set(node, {locDom->xdd(nidx),locDom->ydd(nidx),locDom->zdd(nidx)});
    	     f.set(node, {locDom->fx(nidx),locDom->fy(nidx),locDom->fz(nidx)});
    	     nodalMass.set(node, locDom->nodalMass(nidx));
    	     symm.set(node, {symmX[nidx],symmY[nidx],symmZ[nidx]});
    	     ++nidx ;
       }
     }
   }

   // The fields of the elems set
   FieldRef<int,3,6> elemBC 	= elems.addField<int,3,6>("elemBC");	/* symmetry/free-surface flags for each elem face */
   FieldRef<double,3> dxyz     	= elems.addField<double,3>("dxyz");		/* principal strains -- temporary */
   FieldRef<double,3> delvel  	= elems.addField<double,3>("delvel");	/* velocity gradient -- temporary */
   FieldRef<double,3> delx  	= elems.addField<double,3>("delx");		/* coordinate gradient -- temporary */
   FieldRef<double> e     		= elems.addField<double>("e");   		/* energy */
   FieldRef<double> p     		= elems.addField<double>("p");   		/* pressure */
   FieldRef<double> q     		= elems.addField<double>("q");   		/* q */
   FieldRef<double> ql    		= elems.addField<double>("ql");  		/* linear term for q */
   FieldRef<double> qq    		= elems.addField<double>("qq");   		/* quadratic term for q */
   FieldRef<double> v     		= elems.addField<double>("v");   		/* relative volume */
   FieldRef<double> volo     	= elems.addField<double>("volo");   	/* reference volume */
   FieldRef<double> vnew   		= elems.addField<double>("vnew");   	/* new relative volume -- temporary */
   FieldRef<double> delv     	= elems.addField<double>("delv");   	/* m_vnew - m_v */
   FieldRef<double> vdov     	= elems.addField<double>("vdov");    	/* volume derivative over volume */
   FieldRef<double> arealg     	= elems.addField<double>("arealg");   	/* characteristic length of an element */
   FieldRef<double> ss     		= elems.addField<double>("ss");        	/* "sound speed" */
   FieldRef<double> elemMass    = elems.addField<double>("elemMass");   /* mass */

   // initialize elems values with Domain
   std::vector<ElementRef> elemRefs;
   nidx = 0 ;
   for (Index_t plane=0; plane<edgeElems; ++plane) {
     for (Index_t row=0; row<edgeElems; ++row) {
       for (Index_t col=0; col<edgeElems; ++col) {
    	     ElementRef elem = elems.add(nodeRefs[locDom->nodelist(nidx)[0]], nodeRefs[locDom->nodelist(nidx)[1]],
										 nodeRefs[locDom->nodelist(nidx)[2]], nodeRefs[locDom->nodelist(nidx)[3]],
										 nodeRefs[locDom->nodelist(nidx)[4]], nodeRefs[locDom->nodelist(nidx)[5]],
										 nodeRefs[locDom->nodelist(nidx)[6]], nodeRefs[locDom->nodelist(nidx)[7]]);
    	     elemRefs.push_back(elem);
    	     //dxyz.set(elem, {locDom->dxx(nidx),locDom->dyy(nidx),locDom->dzz(nidx)});
    	     //delvel.set(elem, {locDom->delv_xi(nidx),locDom->delv_eta(nidx),locDom->delv_zeta(nidx)});
    	     //delx.set(elem, {locDom->delx_xi(nidx),locDom->delx_eta(nidx),locDom->delx_zeta(nidx)});
    	     e.set(elem,{locDom->e(nidx)});
    	     p.set(elem,{locDom->p(nidx)});
    	     q.set(elem,{locDom->q(nidx)});
    	     ql.set(elem,{locDom->ql(nidx)});
    	     qq.set(elem,{locDom->qq(nidx)});
    	     v.set(elem,{locDom->v(nidx)});
    	     volo.set(elem,{locDom->volo(nidx)});
    	     delv.set(elem,{locDom->delv(nidx)});
    	     vdov.set(elem,{locDom->vdov(nidx)});
    	     arealg.set(elem,{locDom->arealg(nidx)});
    	     ss.set(elem,{locDom->ss(nidx)});
    	     elemMass.set(elem,{locDom->elemMass(nidx)});
    	     Int_t bcMask = locDom->elemBC(nidx) ;
    	     std::vector<int> maskxim(3,0);
    	     switch (bcMask & XI_M) {
    	     	 case XI_M_COMM: maskxim[0]=0;maskxim[1]=0;maskxim[2]=1; break;
    	         case 0:         break;
    	         case XI_M_SYMM: maskxim[0]=0;maskxim[1]=1;maskxim[2]=0; break;
    	         case XI_M_FREE: maskxim[0]=1;maskxim[1]=0;maskxim[2]=0; break ;
    	     }
    	     std::vector<int> maskxip(3,0);
    	      switch (bcMask & XI_P) {
    	         case XI_P_COMM: maskxip[0]=0;maskxip[1]=0;maskxip[2]=1; break;
    	         case 0:         break ;
    	         case XI_P_SYMM: maskxip[0]=0;maskxip[1]=1;maskxip[2]=0; break ;
    	         case XI_P_FREE: maskxip[0]=1;maskxip[1]=0;maskxip[2]=0; break ;
    	      }
     	     std::vector<int> masketam(3,0);
    	      switch (bcMask & ETA_M) {
    	         case ETA_M_COMM: masketam[0]=0;masketam[1]=0;masketam[2]=1; break;
    	         case 0:          break ;
    	         case ETA_M_SYMM: masketam[0]=0;masketam[1]=1;masketam[2]=0; break ;
    	         case ETA_M_FREE: masketam[0]=1;masketam[1]=0;masketam[2]=0; break ;
    	      }
     	     std::vector<int> masketap(3,0);
    	      switch (bcMask & ETA_P) {
    	         case ETA_P_COMM: masketap[0]=0;masketap[1]=0;masketap[2]=1; break;
    	         case 0:          break ;
    	         case ETA_P_SYMM: masketap[0]=0;masketap[1]=1;masketap[2]=0; break ;
    	         case ETA_P_FREE: masketap[0]=1;masketap[1]=0;masketap[2]=0; break ;
    	      }
     	     std::vector<int> maskzetam(3,0);
    	      switch (bcMask & ZETA_M) {
    	         case ZETA_M_COMM: maskzetam[0]=0;maskzetam[1]=0;maskzetam[2]=1; break;
    	         case 0:           break ;
    	         case ZETA_M_SYMM: maskzetam[0]=0;maskzetam[1]=1;maskzetam[2]=0; break ;
    	         case ZETA_M_FREE: maskzetam[0]=1;maskzetam[1]=0;maskzetam[2]=0; break ;
    	      }
     	     std::vector<int> maskzetap(3,0);
    	      switch (bcMask & ZETA_P) {
    	         case ZETA_P_COMM: maskzetap[0]=0;maskzetap[1]=0;maskzetap[2]=1; break;
    	         case 0:           break ;
    	         case ZETA_P_SYMM: maskzetap[0]=0;maskzetap[1]=1;maskzetap[2]=0; break ;
    	         case ZETA_P_FREE: maskzetap[0]=1;maskzetap[1]=0;maskzetap[2]=0; break ;
    	      }
//    	      elemBC.set(elem,{0,0,0,0,0,0, 0,0,0,0,0,0, 0,1,0,0,0,0});
    	      elemBC.set(elem,{maskxim[2],maskxip[2],masketam[2],masketap[2],maskzetam[2],maskzetap[2],
    	    		  	  	   maskxim[1],maskxip[1],masketam[1],masketap[1],maskzetam[1],maskzetap[1],
							   maskxim[0],maskxip[0],masketam[0],masketap[0],maskzetam[0],maskzetap[0]});
    	      ++nidx ;
       }
     }
   }

   /* element connectivity across each face */
   nidx=0;
   std::vector<ElementRef> connectRefs;
   for (Index_t plane=0; plane<edgeElems; ++plane) {
     for (Index_t row=0; row<edgeElems; ++row) {
       for (Index_t col=0; col<edgeElems; ++col) {
    	   ElementRef connect = connects.add(elemRefs[nidx],
    			   	   	   	    elemRefs[locDom->lxim(nidx)],elemRefs[locDom->lxip(nidx)],
								elemRefs[locDom->letam(nidx)],elemRefs[locDom->letap(nidx)],
								elemRefs[locDom->lzetam(nidx)],elemRefs[locDom->lzetap(nidx)]);
    	   connectRefs.push_back(connect);
    	   ++nidx;
       }
     }
   }
	simit::Tensor<int,2> cycle;
	cycle(0)= locDom->cycle();
	simit::Tensor<double,2> time, deltatime, stoptime, dtfixed,
							dtcourant, dthydro, deltatimemultlb,
							deltatimemultub, dtmax,
							hgcoef, u_cut, qstop, monoq_limiter_mult,
							monoq_max_slope, qlc_monoq, qqc_monoq,
							eosvmin, eosvmax, v_cut, qqc, dvovmax,
							rho0, ss4o3, e_cut, p_cut, q_cut, emin, pmin;
	time(0)=locDom->time();
	deltatime(0)=locDom->deltatime();
	stoptime(0)=locDom->stoptime();
	dtfixed(0)=locDom->dtfixed();
	dtcourant(0)=locDom->dtcourant();
	dthydro(0)=locDom->dthydro();
	deltatimemultlb(0)=locDom->deltatimemultlb();
	deltatimemultub(0)=locDom->deltatimemultub();
	dtmax(0)=locDom->dtmax();
	hgcoef(0)=locDom->hgcoef();
	u_cut(0)=locDom->u_cut();
	qstop(0)=locDom->qstop();
	monoq_limiter_mult(0)=locDom->monoq_limiter_mult();
	monoq_max_slope(0)=locDom->monoq_max_slope();
	qlc_monoq(0)=locDom->qlc_monoq();
	qqc_monoq(0)=locDom->qqc_monoq();
	eosvmin(0)=locDom->eosvmin();
	eosvmax(0)=locDom->eosvmax();
	v_cut(0)=locDom->v_cut();
	qqc(0)=locDom->qqc();
	dvovmax(0)=locDom->dvovmax();
	rho0(0)=locDom->refdens();
	ss4o3(0)=locDom->ss4o3();
	e_cut(0)=locDom->e_cut();
	p_cut(0)=locDom->p_cut();
	q_cut(0)=locDom->q_cut();
	emin(0)=locDom->emin();
	pmin(0)=locDom->pmin();

	std::string codefile = "../lulesh.sim";
   // Compile program and bind arguments
   Program program;
   program.loadFile(codefile);

   Function TimeIncrement_sim = program.compile("TimeIncrement_sim");
   TimeIncrement_sim.bind("nodes",  &nodes);
   TimeIncrement_sim.bind("elems", &elems);
   TimeIncrement_sim.bind("connects", &connects);
   TimeIncrement_sim.bind("cycle", &cycle);
   TimeIncrement_sim.bind("time", &time);
   TimeIncrement_sim.bind("deltatime", &deltatime);
   TimeIncrement_sim.bind("stoptime", &stoptime);
   TimeIncrement_sim.bind("dtfixed", &dtfixed);
   TimeIncrement_sim.bind("dtcourant", &dtcourant);
   TimeIncrement_sim.bind("dthydro", &dthydro);
   TimeIncrement_sim.bind("deltatimemultlb", &deltatimemultlb);
   TimeIncrement_sim.bind("deltatimemultub", &deltatimemultub);
   TimeIncrement_sim.bind("dtmax", &dtmax);
   TimeIncrement_sim.bind("hgcoef", &hgcoef);
   TimeIncrement_sim.bind("u_cut", &u_cut);
   TimeIncrement_sim.bind("qstop", &qstop);
   TimeIncrement_sim.bind("monoq_limiter_mult", &monoq_limiter_mult);
   TimeIncrement_sim.bind("monoq_max_slope", &monoq_max_slope);
   TimeIncrement_sim.bind("qlc_monoq", &qlc_monoq);
   TimeIncrement_sim.bind("qqc_monoq", &qqc_monoq);
   TimeIncrement_sim.bind("eosvmin", &eosvmin);
   TimeIncrement_sim.bind("eosvmax", &eosvmax);
   TimeIncrement_sim.bind("v_cut", &v_cut);
   TimeIncrement_sim.bind("qqc", &qqc);
   TimeIncrement_sim.bind("dvovmax", &dvovmax);
   TimeIncrement_sim.bind("rho0", &rho0);
   TimeIncrement_sim.bind("ss4o3", &ss4o3);
   TimeIncrement_sim.bind("e_cut", &e_cut);
   TimeIncrement_sim.bind("p_cut", &p_cut);
   TimeIncrement_sim.bind("q_cut", &q_cut);
   TimeIncrement_sim.bind("emin", &emin);
   TimeIncrement_sim.bind("pmin", &pmin);

   Function LagrangeLeapFrog_sim = program.compile("LagrangeLeapFrog_sim");
   LagrangeLeapFrog_sim.bind("nodes",  &nodes);
   LagrangeLeapFrog_sim.bind("elems", &elems);
   LagrangeLeapFrog_sim.bind("connects", &connects);
   LagrangeLeapFrog_sim.bind("cycle", &cycle);
   LagrangeLeapFrog_sim.bind("time", &time);
   LagrangeLeapFrog_sim.bind("deltatime", &deltatime);
   LagrangeLeapFrog_sim.bind("stoptime", &stoptime);
   LagrangeLeapFrog_sim.bind("dtfixed", &dtfixed);
   LagrangeLeapFrog_sim.bind("dtcourant", &dtcourant);
   LagrangeLeapFrog_sim.bind("dthydro", &dthydro);
   LagrangeLeapFrog_sim.bind("deltatimemultlb", &deltatimemultlb);
   LagrangeLeapFrog_sim.bind("deltatimemultub", &deltatimemultub);
   LagrangeLeapFrog_sim.bind("dtmax", &dtmax);
   LagrangeLeapFrog_sim.bind("hgcoef", &hgcoef);
   LagrangeLeapFrog_sim.bind("u_cut", &u_cut);
   LagrangeLeapFrog_sim.bind("qstop", &qstop);
   LagrangeLeapFrog_sim.bind("monoq_limiter_mult", &monoq_limiter_mult);
   LagrangeLeapFrog_sim.bind("monoq_max_slope", &monoq_max_slope);
   LagrangeLeapFrog_sim.bind("qlc_monoq", &qlc_monoq);
   LagrangeLeapFrog_sim.bind("qqc_monoq", &qqc_monoq);
   LagrangeLeapFrog_sim.bind("eosvmin", &eosvmin);
   LagrangeLeapFrog_sim.bind("eosvmax", &eosvmax);
   LagrangeLeapFrog_sim.bind("v_cut", &v_cut);
   LagrangeLeapFrog_sim.bind("qqc", &qqc);
   LagrangeLeapFrog_sim.bind("dvovmax", &dvovmax);
   LagrangeLeapFrog_sim.bind("rho0", &rho0);
   LagrangeLeapFrog_sim.bind("ss4o3", &ss4o3);
   LagrangeLeapFrog_sim.bind("e_cut", &e_cut);
   LagrangeLeapFrog_sim.bind("p_cut", &p_cut);
   LagrangeLeapFrog_sim.bind("q_cut", &q_cut);
   LagrangeLeapFrog_sim.bind("emin", &emin);
   LagrangeLeapFrog_sim.bind("pmin", &pmin);

   TimeIncrement_sim.init();
   LagrangeLeapFrog_sim.init();

//   std::cout << locDom->x(locDom->nodelist(0)[7]) << std::endl;
//   std::cout << locDom->y(locDom->nodelist(0)[7]) << std::endl;
//   std::cout << locDom->z(locDom->nodelist(0)[7]) << std::endl;

   // BEGIN timestep to solution */
   timeval start;
   gettimeofday(&start, NULL) ;

   while((locDom->time() < locDom->stoptime()) && (locDom->cycle() < 4)) {

//	   std::cout << " set coord " << coord << std::endl;
//	   std::cout << " End of CPP " << std::endl;
	   TimeIncrement_sim.run() ;
      locDom->cycle() = cycle(0);
      locDom->time() = time(0);
      locDom->deltatime() = deltatime(0);
//  	stoptime(0)=locDom->stoptime();
//  	dtfixed(0)=locDom->dtfixed();
//  	deltatimemultlb(0)=locDom->deltatimemultlb();
//  	deltatimemultub(0)=locDom->deltatimemultub();
//  	dtmax(0)=locDom->dtmax();
//  	hgcoef(0)=locDom->hgcoef();
//  	u_cut(0)=locDom->u_cut();
//  	qstop(0)=locDom->qstop();
//  	monoq_limiter_mult(0)=locDom->monoq_limiter_mult();
//  	monoq_max_slope(0)=locDom->monoq_max_slope();
//  	qlc_monoq(0)=locDom->qlc_monoq();
//  	qqc_monoq(0)=locDom->qqc_monoq();
//  	eosvmin(0)=locDom->eosvmin();
//  	eosvmax(0)=locDom->eosvmax();
//  	v_cut(0)=locDom->v_cut();
//  	qqc(0)=locDom->qqc();
//  	dvovmax(0)=locDom->dvovmax();
      //LagrangeLeapFrog_sim.bind("deltatime", Tensor<double>(double(locDom->deltatime())));
      LagrangeLeapFrog_sim.run();
      locDom->dtcourant() = dtcourant(0);
      locDom->dthydro() = dthydro(0);

      if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank == 0)) {
         printf("cycle = %d, time = %e, dt=%e\n",
                locDom->cycle(), double(locDom->time()), double(locDom->deltatime()) ) ;
      }
   }

   // Use reduced max elapsed time
   double elapsed_time;
   timeval end;
   gettimeofday(&end, NULL) ;
   elapsed_time = (double)(end.tv_sec - start.tv_sec) + ((double)(end.tv_usec - start.tv_usec))/1000000 ;
   double elapsed_timeG;
   elapsed_timeG = elapsed_time;

   // Write out final viz file */
   //if (opts.viz) {
   //   DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   //}
   
   if ((myRank == 0) && (opts.quiet == 0)) {
      VerifyAndWriteFinalOutput(elapsed_timeG, *locDom, opts.nx, numRanks);
   }

   return 0 ;
}
