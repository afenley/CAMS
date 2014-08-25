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

/* the following function computes the angle stresses
 */
void angleforce(struct topology *mol)
{
	int i;
	int atomi,atomj,atomk;
	float theta,stressij,stressik,stressjk,prefactor,rij,rjk,rik;
	float forceij,forceik,forcejk;
	int numangles;
	struct angle *angles;
	struct angle *angles2;
	float **force;
	float sintheta,costheta;

	numangles = mol->numangles;
	angles    = mol->angles;
	force     = mol->force;
	
	for(i=0;i<numangles;i++)
	{
		angles2 = &angles[i];
		atomi = angles2->atomi;
		atomj = angles2->atomj;
		atomk = angles2->atomk;

		/* angle stress is computed as:
		 * stressij = [2*k*(theta - theta0)/sin(theta)]*[cos(theta) - rij/rjk)]
		 * stressik = [2*k*(theta - theta0)/sin(theta)]*[rik^2/(rij*rjk)]
		 * stressjk = [2*k*(theta - theta0)/sin(theta)]*[cos(theta) - rjk/rij]
		 */
	
		costheta = cosangle(mol,atomi,atomj,atomk);
		theta    = acos(costheta);
		sintheta = sqrt(1-(costheta*costheta));

		rij = distance(mol,atomi,atomj);
		rik = distance(mol,atomi,atomk);
		rjk = distance(mol,atomj,atomk);

		prefactor = (angles2->fc)*(theta - angles2->angle)/sintheta;
		prefactor = prefactor + prefactor;

		forceij = prefactor * (costheta - rij/rjk)/rij;
		forceik = prefactor * (rik*rik)/(rij*rjk)/rik;
		forcejk = prefactor * (costheta - rjk/rij)/rjk;
	
		force[atomi][0] = force[atomi][0] + forceij* (mol->crd[atomj][0] - mol->crd[atomi][0])/rij;
		force[atomi][1] = force[atomi][1] + forceij* (mol->crd[atomj][1] - mol->crd[atomi][1])/rij;
		force[atomi][2] = force[atomi][2] + forceij* (mol->crd[atomj][2] - mol->crd[atomi][2])/rij;
		force[atomj][0] = force[atomj][0] - forceij* (mol->crd[atomj][0] - mol->crd[atomi][0])/rij;
		force[atomj][1] = force[atomj][1] - forceij* (mol->crd[atomj][1] - mol->crd[atomi][1])/rij;
		force[atomj][2] = force[atomj][2] - forceij* (mol->crd[atomj][2] - mol->crd[atomi][2])/rij;
	
		force[atomi][0] = force[atomi][0] + forceik* (mol->crd[atomk][0] - mol->crd[atomi][0])/rik;
		force[atomi][1] = force[atomi][1] + forceik* (mol->crd[atomk][1] - mol->crd[atomi][1])/rik;
		force[atomi][2] = force[atomi][2] + forceik* (mol->crd[atomk][2] - mol->crd[atomi][2])/rik;
		force[atomk][0] = force[atomk][0] - forceik* (mol->crd[atomk][0] - mol->crd[atomi][0])/rik;
		force[atomk][1] = force[atomk][1] - forceik* (mol->crd[atomk][1] - mol->crd[atomi][1])/rik;
		force[atomk][2] = force[atomk][2] - forceik* (mol->crd[atomk][2] - mol->crd[atomi][2])/rik;

		force[atomj][0] = force[atomj][0] + forcejk* (mol->crd[atomk][0] - mol->crd[atomj][0])/rjk;
		force[atomj][1] = force[atomj][1] + forcejk* (mol->crd[atomk][1] - mol->crd[atomj][1])/rjk;
		force[atomj][2] = force[atomj][2] + forcejk* (mol->crd[atomk][2] - mol->crd[atomj][2])/rjk;
		force[atomk][0] = force[atomk][0] - forcejk* (mol->crd[atomk][0] - mol->crd[atomj][0])/rjk;
		force[atomk][1] = force[atomk][1] - forcejk* (mol->crd[atomk][1] - mol->crd[atomj][1])/rjk;
		force[atomk][2] = force[atomk][2] - forcejk* (mol->crd[atomk][2] - mol->crd[atomj][2])/rjk;
	}
}
