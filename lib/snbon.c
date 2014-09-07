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

/* this function computes the lj and coulombic interaction energy
 */
void nbstress(struct topology *mol)
{
	float stressij,rij,rij2,rij6,c6,c12;
	int *nbcount,*ljtypes,*nblist2;
	int **nblist;
	float **forcefield;
	float *vdwstress,*qqstress,**chargeij;
	int i,j,ljtype,atomi, atomj,currentcount,numatoms,numljtypes;

	numatoms   = mol->numatoms;
	nbcount    = mol->nbcount;
	ljtypes    = mol->ljtypes;
	numljtypes = mol->numljtypes;
	nblist     = mol->nblist;
	forcefield = mol->forcefield;
	vdwstress  = mol->vdwstress;
	qqstress   = mol->qqstress;
	chargeij   = mol->chargeij;

	for(i=0;i<numatoms;i++)
	{
		atomi = i;
		currentcount = nbcount[i];
		nblist2 = nblist[i];

		for(j=0;j<currentcount;j++)
		{
			atomj = nblist2[j];

			ljtype = ljtypes[atomi]*numljtypes + ljtypes[atomj];
			rij2 = 1.0/distance2(mol,atomi,atomj);
			rij6 = rij2*rij2*rij2;
			
			stressij         = -(forcefield[ljtype][1])*rij6*rij6 + (forcefield[ljtype][0])*rij6;
			vdwstress[atomi] = vdwstress[atomi] - stressij;// negative sign to represent hydrostatic pressure

			stressij        = -332.06*chargeij[atomi][atomj]*sqrt(rij2);
			qqstress[atomi] = qqstress[atomi] - stressij;// negative sign to represent hydrostatic pressure
			
		}

	}
}
