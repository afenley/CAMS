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


/* the following function computes the improper dihedral energy
 */
float idihedralenergy(struct topology *mol)
{
	int i,atomi,atomj,atomk,atoml;
	float ij_x, ij_y, ij_z;
	float ik_x, ik_y, ik_z;
	float il_x, il_y, il_z;
	float jk_x, jk_y, jk_z;
	float jl_x, jl_y, jl_z;
	float kl_x, kl_y, kl_z;
	float n1_x, n1_y, n1_z;
	float n2_x, n2_y, n2_z;
	float sign, n1, n2, n12;
	float energy,theta;

	energy = 0.0;
	for(i=0;i<mol->numidihedrals;i++)
	{
		atomi = mol->idihedrals[i].atomi;
		atomj = mol->idihedrals[i].atomj;
		atomk = mol->idihedrals[i].atomk;
		atoml = mol->idihedrals[i].atoml;

		ij_x = mol->crd[atomj][0] - mol->crd[atomi][0]; ik_x = mol->crd[atomk][0] - mol->crd[atomi][0]; il_x = mol->crd[atoml][0] - mol->crd[atomi][0];
		ij_y = mol->crd[atomj][1] - mol->crd[atomi][1]; ik_y = mol->crd[atomk][1] - mol->crd[atomi][1]; il_y = mol->crd[atoml][1] - mol->crd[atomi][1];
		ij_z = mol->crd[atomj][2] - mol->crd[atomi][2]; ik_z = mol->crd[atomk][2] - mol->crd[atomi][2]; il_z = mol->crd[atoml][2] - mol->crd[atomi][2];

		jk_x = mol->crd[atomk][0] - mol->crd[atomj][0]; jl_x = mol->crd[atoml][0] - mol->crd[atomj][0];
		jk_y = mol->crd[atomk][1] - mol->crd[atomj][1]; jl_y = mol->crd[atoml][1] - mol->crd[atomj][1];
		jk_z = mol->crd[atomk][2] - mol->crd[atomj][2]; jl_z = mol->crd[atoml][2] - mol->crd[atomj][2];

		kl_x = mol->crd[atoml][0] - mol->crd[atomk][0];
		kl_y = mol->crd[atoml][1] - mol->crd[atomk][1];
		kl_z = mol->crd[atoml][2] - mol->crd[atomk][2];

		n1_x = ij_y*jk_z - ij_z*jk_y;
		n1_y = ij_z*jk_x - ij_x*jk_z;
		n1_z = ij_x*jk_y - ij_y*jk_x;

		n2_x = jk_y*kl_z - jk_z*kl_y;
		n2_y = jk_z*kl_x - jk_x*kl_z;
		n2_z = jk_x*kl_y - jk_y*kl_x;

		sign = (n1_y*n2_z - n1_z*n2_y)*jk_x + (n1_z*n2_x - n1_x*n2_z)*jk_y + (n1_x*n2_y - n1_y*n2_x)*jk_z;
		if (sign < 0)
		{
			sign = -1.0;
		}
		else
		{
			sign = 1.0;
		}

		n1 = sqrt(n1_x*n1_x + n1_y*n1_y + n1_z*n1_z);
		n2 = sqrt(n2_x*n2_x + n2_y*n2_y + n2_z*n2_z);
		n12 = n1_x*n2_x + n1_y*n2_y + n1_z*n2_z;
		if (n12/(n1*n2) >= 1.0)
		{
			theta = 0.0;
		}
		else
		{
			if (n12/(n1*n2) <= -1.0)
			{
				theta = M_PI;
			}
			else
			{
				theta = acos(n12/(n1*n2));
			}
		}

		energy = energy + mol->idihedrals[i].fc*(1.0+cos(2.0*theta*sign - mol->idihedrals[i].phi));
	}
	return energy;
}
