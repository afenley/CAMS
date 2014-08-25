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
#include<math.h>
#include"../include/mdmath.h"
#include"../include/structdefs.h"

/* smoothing the stresses
 */
void normalization(struct topology *mol)
{
	float normalization1,normalization2;
	int   i,j,atomi,numatoms;
	float *bondstress,*anglestress,*dihedralstress,*vdwstress,*qqstress,*solventstress,*kineticstress;

	bondstress      = mol->bondstress;
	anglestress     = mol->anglestress;
	dihedralstress  = mol->dihedralstress;
	vdwstress       = mol->vdwstress;
	qqstress        = mol->qqstress;
	solventstress   = mol->solventstress;
	kineticstress   = mol->kineticstress;
        normalization1  = 11.579500245/(mol->volume);  /* hydrostatic pressure (kbar): 69.47700/(6.0 * mol->volume) */
        normalization2  = 23.159000470/(mol->volume);  /* hydrostatic pressure (kbar): 69.47700/(3.0 * mol->volume) */
/*	normalization1  = 16.6667/(mol->volume);*/  /* normalization1 = 100.0/(6.0 * mol->volume) */
/*	normalization2  = 33.3333/(mol->volume);*/  /* normalization2 = 100.0/(3.0 * mol->volume) */

	numatoms = mol->numatoms;
	for(atomi=0;atomi<numatoms;atomi++)
	{
		bondstress[atomi]     = (bondstress[atomi])*normalization1;
		anglestress[atomi]    = (anglestress[atomi])*normalization1;
		dihedralstress[atomi] = (dihedralstress[atomi])*normalization1;
		vdwstress[atomi]      = (vdwstress[atomi])*normalization1;
		qqstress[atomi]       = (qqstress[atomi])*normalization1;
		solventstress[atomi]  = (solventstress[atomi])*normalization1;
		kineticstress[atomi]  = (kineticstress[atomi])*normalization2;
	}
}
