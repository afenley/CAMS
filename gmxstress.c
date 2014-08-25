/*  
    Calculator of Atomistic Mechanical Stress (CAMS) is an open-source 
    software package for computing atomic stresses from molecular dynamics 
    simulations.

    Copyright (C) 2014  Andrew T. Fenley (afenley@ucsd.edu), Hari M. Muddana 
    (hmuddana@gmail.com), and Michael K. Gilson (mgilson@ucsd.edu).

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<sys/time.h>
#include"include/xdrfile.h"
#include"include/xdrfile_trr.h"
#include"include/xdrfile_xtc.h"
#include"include/structdefs.h"
#include"include/gmxlib.h"
#include"include/energy.h"
#include"include/stress.h"
#define MAX_FILENAME 200

char   ndxfilename[MAX_FILENAME];
char   tprfilename[MAX_FILENAME];
char   trrfilename[MAX_FILENAME];
char   outfilename[MAX_FILENAME];
int    debugflag=0;
int    gbflag=0;
float  gbsalt=0.145;
int    bondiflag=0;
int    kineticflag=0;
int    splitflag=0;

int main(int argc, char *argv[])
{
	int i,j,natoms;
	char command[500];
	struct topology molecule;
	int numframes, count,atomindex,resindex,prevresindex;
	XDRFILE *infile;
	timeval start,stop;
	gettimeofday(&start,0);
	
	help();                                  //print help menu
	commandlineargs(argc, argv);             //read command line arguments
	readindexgroup(ndxfilename,&molecule);   //read index group
	readtpr(tprfilename,&molecule);	         //read topology file
	molecule.setchargeij();
	if(bondiflag || molecule.missingradii())
	{
		printf("Atomic radii missing or user has set the bondi flag. Now, using Bondi Radii\n");
		molecule.setradii();                     //set atomic radii to Bondi Radii
	}


	numframes = 0;
	read_trr_natoms(trrfilename,&natoms);
	if(natoms != molecule.numatoms)
	{
		printf("ERROR: the number of atoms in the trajectory file and the tpr file do not match. exiting now.\n"); exit(1);
	}
	infile = xdrfile_open(trrfilename,"r");
	if(debugflag)
	{
		printf("%12s %12s %12s %12s %12s %12s %12s %12s %12s\n","BOND","ANGLE","PDIH","IDIH","RBDIH","PAIR","NBOND","KINETIC","GB(OBC)");
	}
	while(readframe(infile,&molecule))
	{
		numframes = numframes + 1;
		printf("Frame No: %d\n",numframes);
		getnblist_grid(&molecule);
		if(debugflag)
		{
			printf("%12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f\n",bondenergy(&molecule),angleenergy(&molecule),pdihedralenergy(&molecule),idihedralenergy(&molecule),
									rbdihedralenergy(&molecule),pairenergy(&molecule),nbenergy(&molecule),kineticenergy(&molecule),gbenergy(&molecule));
		}
		initializestress(&molecule);
		bondstress(&molecule);
		anglestress(&molecule);
		pdihedralstress(&molecule);
		idihedralstress(&molecule);
		rbdihedralstress(&molecule);	
		pairstress(&molecule);
		nbstress(&molecule);
		if(gbflag)
		{
			gbstress(&molecule);
		}
		if(kineticflag)
		{
			kineticstress(&molecule);
		}
		normalization(&molecule);
		writeout(outfilename,&molecule);
	}
	xdrfile_close(infile);
	writestat(outfilename,&molecule,numframes);
	gettimeofday(&stop,0);
	printf("calculation time: %d seconds and %d milliseconds\n",stop.tv_sec - start.tv_sec, (stop.tv_usec - start.tv_usec)/1000);
}
