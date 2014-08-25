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

/* the following function computes the improper dihedral stress
 */
void idihedralstress(struct topology *mol)
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
	float stressij,stressik,stressil,stressjk,stressjl,stresskl,theta;
	float rij,rik,ril,rjk,rjl,rkl,temp1,temp2,temp3,prefactor;

	int numidihedrals;
	struct idihedral *idihedrals;
	float (*crd)[3];
	float *stress;

	numidihedrals = mol->numidihedrals;
	idihedrals    = mol->idihedrals;
	crd           = mol->crd;
	stress        = mol->dihedralstress;

	for(i=0;i<numidihedrals;i++)
	{
		atomi = idihedrals[i].atomi;
		atomj = idihedrals[i].atomj;
		atomk = idihedrals[i].atomk;
		atoml = idihedrals[i].atoml;
		ij_x = crd[atomj][0] - crd[atomi][0]; ik_x = crd[atomk][0] - crd[atomi][0]; il_x = crd[atoml][0] - crd[atomi][0];
		ij_y = crd[atomj][1] - crd[atomi][1]; ik_y = crd[atomk][1] - crd[atomi][1]; il_y = crd[atoml][1] - crd[atomi][1];
		ij_z = crd[atomj][2] - crd[atomi][2]; ik_z = crd[atomk][2] - crd[atomi][2]; il_z = crd[atoml][2] - crd[atomi][2];

		jk_x = crd[atomk][0] - crd[atomj][0]; jl_x = crd[atoml][0] - crd[atomj][0];
		jk_y = crd[atomk][1] - crd[atomj][1]; jl_y = crd[atoml][1] - crd[atomj][1];
		jk_z = crd[atomk][2] - crd[atomj][2]; jl_z = crd[atoml][2] - crd[atomj][2];

		kl_x = crd[atoml][0] - crd[atomk][0];
		kl_y = crd[atoml][1] - crd[atomk][1];
		kl_z = crd[atoml][2] - crd[atomk][2];

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
		
		rij = distance(mol,atomi,atomj);
		rik = distance(mol,atomi,atomk);
		rjk = distance(mol,atomj,atomk);
		ril = distance(mol,atomi,atoml);
		rjl = distance(mol,atomj,atoml);
		rkl = distance(mol,atomk,atoml);
			
                if (idihedrals[i].phi > 3.14) // phi = PI, switch sign of sin
                {
			prefactor = sign*-4.0*(idihedrals[i].fc)*cos(theta*sign); // sin(2*theta) = sin(theta)*2*cos(theta)
                        //printf ("PHI: %f\n",idihedrals[i].phi);
	        }
                else // phi = 0
                {
                        prefactor = sign*4.0*(idihedrals[i].fc)*cos(theta*sign); // sin(2*theta) = sin(theta)*2*cos(theta)
                        //printf ("PHI: %f\n",idihedrals[i].phi);
                }
 	  	stressij = -rij*rij/n1;
		temp1 = (jk_x*kl_x + jk_y*kl_y + jk_z*kl_z)/n2;
		temp2 = cos(theta)*(jk_x*ik_x + jk_y*ik_y + jk_z*ik_z)/n1;
		stressij = stressij*(temp1 + temp2);

		stressik = rik*rik/n1;
		temp1 = (jk_x*jl_x + jk_y*jl_y + jk_z*jl_z)/n2;
		temp2 = cos(theta)*(ij_x*jk_x + ij_y*jk_y + ij_z*jk_z)/n1;
		stressik = stressik*(temp1+temp2);

		stressil = (-rjk*rjk*ril*ril)/(n1*n2);

		stressjk = -rjk*rjk;
		temp1 = ((ik_x*jl_x + ik_y*jl_y + ik_z*jl_z) + (ij_x*kl_x + ij_y*kl_y + ij_z*kl_z))/(n1*n2);
		temp2 = cos(theta)*(ij_x*ik_x + ij_y*ik_y + ij_z*ik_z)/(n1*n1);
		temp3 = cos(theta)*(kl_x*jl_x + kl_y*jl_y + kl_z*jl_z)/(n2*n2);
		stressjk = stressjk*(temp1+temp2+temp3);

		stressjl = rjl*rjl/n2;
		temp1 = (jk_x*ik_x + jk_y*ik_y + jk_z*ik_z)/n1;
		temp2 = cos(theta)*(jk_x*kl_x + jk_y*kl_y + jk_z*kl_z)/n2;
		stressjl = stressjl*(temp1+temp2);

		stresskl = -rkl*rkl/n2;
		temp1 = (ij_x*jk_x + ij_y*jk_y + ij_z*jk_z)/n1; 
		temp2 = cos(theta)*(jk_x*jl_x + jk_y*jl_y + jk_z*jl_z)/n2; 
		stresskl = stresskl*(temp1+temp2);

                // negative sign to represent hydrostatic pressure
		stress[atomi] = stress[atomi] - prefactor*(stressij + stressik + stressil);
		stress[atomj] = stress[atomj] - prefactor*(stressij + stressjk + stressjl);
		stress[atomk] = stress[atomk] - prefactor*(stressik + stressjk + stresskl);
		stress[atoml] = stress[atoml] - prefactor*(stressil + stressjl + stresskl);

	}
	
}
