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
#include "../include/xdrfile.h"
#include "../include/xdrfile_trr.h"
#include "../include/xdrfile_xtc.h"
#include"../include/gmxlib.h"
#include"../include/structdefs.h"

#define NBCUTOFF 10.0
#define I_NBCUTOFF 0.10

extern int  splitflag;

/* this function reads the groups in the .ndx file, and asks the user to select a group.
*/
void readindexgroup(char ndxfilename[],struct topology *mol)
{
	FILE* infile;
	int grpcount,atomcount;
	char buffer[1000];
	int i;
	int indexgroupnr;

	/* output all the index groups and the user selects the groups
	 */
	grpcount = 0;
	infile = fopen(ndxfilename,"r");
	if(infile == NULL) { printf("%s file not found\n",ndxfilename); exit(1);}
	while((fgets(buffer,1000,infile)) != NULL)
	{
		if (buffer[0] == '[')
		{
			printf("Atom group %4d: %s",grpcount,buffer);
			grpcount++;
		}	
	}
	fclose(infile);
	printf("Select the group for output (0 - %d):",grpcount-1);
	scanf("%d",&indexgroupnr);



	/* memory allocation for the index group array is done dynamically.
	 * so we first count the number of atoms in the group before reading
	 * the atom numbers in the group
	 */
	infile = fopen(ndxfilename,"r");
	grpcount = -1;
	while((fgets(buffer,1000,infile)) != NULL)
	{
		if (buffer[0] == '[')
		{
			grpcount++;
		}
		if(grpcount == indexgroupnr) break;
	}
	atomcount = 0;
	while(1)
	{
		if(!fgets(buffer,1000,infile)) break;
		else if(buffer[0] == '[') break;
		else if(buffer[0] == ';' || buffer[0] == '@' || buffer[0] == '\0' || buffer[0] == '\n') continue;
		else { atomcount = atomcount + wordcount(buffer);}
	}
	mol->indexgroupcount = atomcount;
	fclose(infile);



	/* now with the atomcount and the index group at hand, we read the atom
	 * number in the index group
	 * atom numbers in the index file start from 1, whereas the indices in the tpr file
	 * and in this program start from 0. so, we are shifting the indices to start from
	 * zero.
	 */
	mol->initializegroup();
	if(splitflag)
	{
		mol->initializesplit();
        }
	infile = fopen(ndxfilename,"r");
	grpcount = -1;
	while((fgets(buffer,1000,infile)) != NULL)
	{
		if (buffer[0] == '[')
		{
			grpcount++;
		}
		if(grpcount == indexgroupnr) break;
	}
	for(i=0;i<mol->indexgroupcount;i++)
	{
		fscanf(infile,"%d",&mol->indexgroup[i]);
		mol->indexgroup[i] = mol->indexgroup[i] - 1;
	}
	fclose(infile);
}
