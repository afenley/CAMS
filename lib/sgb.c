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


/* The following function computes the GB(OBC) stress
 */
float gbstress(struct topology *mol)
{
	int i,j,n,numatoms,mycount;
	float offset=0.09,alpha=1.0,beta=0.8,gamma=4.85,Ri, Ri_i, Rj, rij, rij2, rij_Rj, rij_Rj2,br2, fgb1, fgb2, stressij, I, IUkk, ILkk;
	float psi, fgb,UL,irij;
	int   *ncount,*mylist;
	int   **nlist;
	float *bradii,*solventstress,*radius,*s_hct,**chargeij;

	numatoms      = mol->numatoms;
	ncount        = mol->ncount;
	nlist         = mol->nlist;
	solventstress = mol->solventstress;
	radius        = mol->radius;
	s_hct         = mol->s_hct;
	chargeij      = mol->chargeij;

	bradii = (float*) malloc(numatoms*sizeof(float));
        for(i=0;i<numatoms;i++)
        {
                I       = 0.0;
                Ri      = radius[i] - offset;
		Ri_i    = 1.0/Ri;
		mycount = ncount[i];
		mylist  = nlist[i];
                for(n=0;n<mycount;n++)
                {
                        j = mylist[n];
                        Rj  = (radius[j] - offset) * (s_hct[j]);
                        rij = distance(mol,i,j);
                        if(i!=j)
                        {
				rij_Rj  = rij + Rj;
				rij_Rj2 = rij - Rj;

                                if(rij_Rj <= Ri) ILkk = 1.0;
                                if((rij_Rj2 <= Ri) && (Ri < rij_Rj)) ILkk = Ri_i;
                                if(Ri <= rij_Rj2) ILkk = 1.0/rij_Rj2;
                                if(rij_Rj <= Ri) IUkk = 1.0;
                                if(Ri < rij_Rj) IUkk = 1.0/rij_Rj;

				UL = (IUkk*IUkk - ILkk*ILkk)*0.25;
				irij = 1.0/rij;
                                I = I + ILkk - IUkk + rij*UL + irij*(0.5*log(IUkk/ILkk) - Rj*Rj*UL);
                        }
                }
                psi       = 0.5*I*Ri;      //I = 0.5*I;
                bradii[i] = 1.0/(Ri_i - (1.0/radius[i])*tanh(psi*(alpha - psi*(beta - gamma*psi))));
        }
	
	for(i=0;i<numatoms;i++)
	{
		mycount = ncount[i];
		mylist  = nlist[i];

		for(n=0;n<mycount;n++)
		{
			j = mylist[n];
			if(i!=j)
			{
				rij2 = distance2(mol,i,j);
				br2  = bradii[i]*bradii[j];
				fgb1 = exp(-rij2/(4.0*br2));
				fgb2 = rij2 + br2*fgb1;
		
				stressij = 163.96*chargeij[i][j]*rij2*(1.0-fgb1/4.0)/sqrt(fgb2*fgb2*fgb2);
				solventstress[i] = solventstress[i] - stressij;// negative sign to represent hydrostatic pressure
				solventstress[j] = solventstress[j] - stressij;
			}
		}
	}


}
