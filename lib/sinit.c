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
#include"../include/structdefs.h"


/* initialize the stresses to zeros
 */
void initializestress(struct topology *mol)
{
	int i;
	int numatoms;
	
	numatoms = mol->numatoms;
	memset(mol->bondstress,0,numatoms*sizeof(float));
	memset(mol->anglestress,0,numatoms*sizeof(float));
	memset(mol->dihedralstress,0,numatoms*sizeof(float));
	memset(mol->vdwstress,0,numatoms*sizeof(float));
	memset(mol->qqstress,0,numatoms*sizeof(float));
	memset(mol->solventstress,0,numatoms*sizeof(float));
	memset(mol->kineticstress,0,numatoms*sizeof(float));
}
