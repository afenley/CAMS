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


/* NOTE: energy calculations of bonded interactions rely on the assumption
 * that the molecule is whole (pbc condition). Periodic images are not taken 
 * into account in computing the energies.
 */


/* this function computes the lj and coulombic interaction energy
 */

float nbenergy(struct topology *mol)
{
	int i,j,ljtype1,ljtype2,ljtype,atomi,atomj;
	float energy,rij,rij2,rij6,rij12,c6,c12;
	
	energy = 0.0;
	for(i=0;i<mol->numatoms;i++)
	{
		for(j=0;j<mol->nbcount[i];j++)
		{
			atomi = i;
			atomj = mol->nblist[i][j];
			
			ljtype1 = mol->ljtypes[atomi];
			ljtype2 = mol->ljtypes[atomj];
			ljtype = ljtype1*mol->numljtypes + ljtype2;
			c6 = (mol->forcefield[ljtype][0])/6.0;
			c12 = (mol->forcefield[ljtype][1])/12.0;

			rij = distance(mol,atomi,atomj);
			rij2 = rij*rij;
			rij6 = rij2*rij2*rij2;
			rij12 = rij6*rij6;
			energy = energy + (c12/rij12) - (c6/rij6);
			energy = energy + 332.06*(mol->charge[atomi])*(mol->charge[atomj])/rij;
		}
	}
	return energy/2.0;
}
