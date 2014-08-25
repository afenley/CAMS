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
#include "../include/xdrfile.h"
#include "../include/xdrfile_trr.h"
#include "../include/xdrfile_xtc.h"
#include"../include/structdefs.h"

extern char trrfilename[];
/* this function reads one frame from the trr file
 */
int readframe(XDRFILE *infile,struct topology *mol)
{
	int numatoms;
	int i,j;
	char buffer[1000];
	char temp[100];
	int tempint;
	float tempfloat;
	int found;
	fpos_t fposition;

	int step;
	float time, lambda, box[3][3];

	read_trr_natoms(trrfilename,&numatoms);
	found = read_trr(infile,numatoms,&step,&time,&lambda,box,mol->crd,mol->velocity,NULL);	
	if(found != exdrOK)
	{
		return 0;
	}
	mol->xbox = 10.0 * box[0][0];
	mol->ybox = 10.0 * box[1][1];
	mol->zbox = 10.0 * box[2][2];
	mol->volume = (mol->xbox)*(mol->ybox)*(mol->zbox);
	mol->volume = (mol->volume)/(numatoms);
	
	if(mol->volume == 0.0)
	{
		mol->volume = 20.58;   //volume of a carbon atom
	}	
	
	for(i=0;i<numatoms;i++)
	{
		mol->crd[i][0] = 10.0 * mol->crd[i][0];
		mol->crd[i][1] = 10.0 * mol->crd[i][1];
		mol->crd[i][2] = 10.0 * mol->crd[i][2];
	}
	return 1;	
}
