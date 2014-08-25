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

/* the following function computes the pair stresses
 */
void pairstress(struct topology *mol)
{
	int i,atomi,atomj;
	float stressij, rij,rij2,rij6,rij12;
	float *vdwstress;
	float *qqstress;
	struct pair *pairs;
	int numpairs;
	float **chargeij;
	float fudgeQQ;

	numpairs  = mol->numpairs;
	vdwstress = mol->vdwstress;
	qqstress  = mol->qqstress;
	pairs     = mol->pairs;
	chargeij  = mol->chargeij;
	fudgeQQ   = (mol->fudgeQQ)*332.06;

	for(i=0;i<numpairs;i++)
	{
		atomi = pairs[i].atomi;
		atomj = pairs[i].atomj;

		/* pair lj and qq stresses are computed as:
		 * LJ: stressij = -12.0*c12/rij12 + 6.0*c6/rij6
		 * QQ: stressij = -fudgeQQ*332.06*qi*qj/rij
		 */
		
		rij2 = 1.0/distance2(mol,atomi,atomj);
		rij6 = rij2*rij2*rij2;

		stressij         = -(pairs[i].c12)*rij6*rij6 + (pairs[i].c6)*rij6;
		vdwstress[atomi] = vdwstress[atomi] - stressij;// negative sign to represent hydrostatic pressure
		vdwstress[atomj] = vdwstress[atomj] - stressij;
		
		stressij        = -fudgeQQ*chargeij[atomi][atomj]*sqrt(rij2);
		qqstress[atomi] = qqstress[atomi] - stressij;// negative sign to represent hydrostatic pressure
		qqstress[atomj] = qqstress[atomj] - stressij;
	}
}
