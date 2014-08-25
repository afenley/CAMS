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


/* The following function computes the GB(OBC) energy
 */
float gbenergy(struct topology *mol)
{
	int i,j,n;
	float offset=0.09;
	float alpha=1.0;
	float beta=0.8;
	float gamma=4.85;
	int   numatoms;
	float Ri, Rj, rij, I, energy, Ukk, Lkk;
	float psi, fgb;
	float *bradii;
	int   *ncount;
	int   **nlist;
	

	numatoms = mol->numatoms;
	ncount   = mol->ncount;
	nlist    = mol->nlist;

	bradii = (float*) malloc(numatoms*sizeof(float));


        /*for(i=0;i<numatoms;i++)
        {
                I   = 0.0;
                Ri  = mol->radius[i] - offset;
                for(n=0;n<ncount[i];n++)
                {
                        j = nlist[i][n];
                        Rj  = (mol->radius[j] - offset) * (mol->s_hct[j]);
                        rij = distance(mol,i,j);
                        if(i!=j)
                        {
                                if((rij + Rj) <= Ri) Lkk = 1.0;
                                if(((rij - Rj) <= Ri) && (Ri < (rij + Rj))) Lkk = Ri;
                                if(Ri <= (rij - Rj)) Lkk = rij - Rj;
                                if((rij + Rj) <= Ri) Ukk = 1.0;
                                if(Ri < (rij + Rj)) Ukk = rij + Rj;

                                I = I + 1.0/Lkk;
                                I = I - 1.0/Ukk;
                                I = I + (rij/4.0)*((1.0/Ukk)*(1.0/Ukk) - (1.0/Lkk)*(1.0/Lkk));
                                I = I + (0.5/rij)*log(Lkk/Ukk);
                                I = I + 0.25*(Rj*Rj/rij)*((1.0/Lkk)*(1.0/Lkk) - (1.0/Ukk)*(1.0/Ukk));
                        }
                }
                I         = 0.5*I;
                psi       = I*(Ri);
                bradii[i] = 1.0/(Ri) - (1.0/(Ri+offset))*tanh(alpha*psi - beta*psi*psi + gamma*psi*psi*psi);
                bradii[i] = 1.0/bradii[i];
        }


        energy = 0.0;
        for(i=0;i<numatoms;i++)
        {
                for(n=0;n<ncount[i];n++)
                {
                        j   = nlist[i][n];
                        rij = distance(mol,i,j);
                        fgb = sqrt(rij*rij + bradii[i]*bradii[j]*exp(-rij*rij/(4.0*bradii[i]*bradii[j])));
                        energy = energy - (166.0 * (1.0 - 1.0/80.0) * (mol->charge[i]) * (mol->charge[j]))/fgb;
                }
        }*/
	
	
	for(i=0;i<numatoms;i++)
	{
		I   = 0.0;
		Ri  = (mol->radius[i]) - offset;
		for(j=0;j<numatoms;j++)
		{
			Rj  = (mol->radius[j] - offset) * (mol->s_hct[j]);
			rij = distance(mol,i,j);
			if(i!=j)
			{
				if((rij + Rj) <= Ri) Lkk = 1.0;
				if(((rij - Rj) <= Ri) && (Ri < (rij + Rj))) Lkk = Ri;
				if(Ri <= (rij - Rj)) Lkk = rij - Rj;
				if((rij + Rj) <= Ri) Ukk = 1.0;
				if(Ri < (rij + Rj)) Ukk = rij + Rj;
				
				I = I + 1.0/Lkk;
				I = I - 1.0/Ukk;
				I = I + (rij/4.0)*((1.0/Ukk)*(1.0/Ukk) - (1.0/Lkk)*(1.0/Lkk));
				I = I + (0.5/rij)*log(Lkk/Ukk);
				I = I + 0.25*(Rj*Rj/rij)*((1.0/Lkk)*(1.0/Lkk) - (1.0/Ukk)*(1.0/Ukk));
			}
		}
		I         = 0.5*I;
		psi       = I*(Ri); 
		bradii[i] = 1.0/(Ri) - (1.0/(Ri+offset))*tanh(alpha*psi - beta*psi*psi + gamma*psi*psi*psi);
		bradii[i] = 1.0/bradii[i];
	}

	
	energy = 0.0;
	for(i=0;i<numatoms;i++)
	{
		for(j=0;j<numatoms;j++)
		{
			rij = distance(mol,i,j);
			fgb = sqrt(rij*rij + bradii[i]*bradii[j]*exp(-rij*rij/(4.0*bradii[i]*bradii[j])));
			energy = energy - (166.0 * (1.0 - 1.0/80.0) * (mol->charge[i]) * (mol->charge[j]))/fgb;
		}
	}
	return energy;
}
