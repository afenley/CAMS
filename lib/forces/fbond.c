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

/* The following function computes the bond stresses
 */
void bondforce(struct topology *mol)
{
	int i;
	float forceij,rij,dx,dy,dz,forceij_x,forceij_y,forceij_z;
	int numbonds;
	struct bond *bonds;
	struct bond *bonds2;
	float **force;
	int atomi,atomj;

	numbonds     = mol->numbonds;
	bonds        = mol->bonds;
	force        = mol->force;

	for(i=0;i<numbonds;i++)
	{
		/* bond force is computed as 2*k*(rij - ro)
		 */	
		bonds2 = &bonds[i];
		atomi = bonds2->atomi;
		atomj = bonds2->atomj;

		dx = mol->crd[atomj][0] - mol->crd[atomi][0];
		dy = mol->crd[atomj][1] - mol->crd[atomi][1];
		dz = mol->crd[atomj][2] - mol->crd[atomi][2];

		rij        = sqrt(dx*dx + dy*dy + dz*dz);
		forceij    = 2.0*(bonds2->fc)*(rij - bonds2->length);
		forceij_x  = forceij*dx/rij;
		forceij_y  = forceij*dy/rij;
		forceij_z  = forceij*dz/rij;

		force[atomi][0] = force[atomi][0] + forceij_x;
		force[atomi][1] = force[atomi][1] + forceij_y;
		force[atomi][2] = force[atomi][2] + forceij_z;

		force[atomj][0] = force[atomj][0] - forceij_x;
		force[atomj][1] = force[atomj][1] - forceij_y;
		force[atomj][2] = force[atomj][2] - forceij_z;
	}
}
