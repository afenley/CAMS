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
void anglestress(struct topology *mol)
{
	int i;
	int atomi,atomj,atomk;
	float theta,stressij,stressik,stressjk,prefactor,rij,rjk,rik;
	int numangles;
	struct angle *angles;
	struct angle *angles2;
	float *stress;
	float sintheta,costheta;

	numangles = mol->numangles;
	angles    = mol->angles;
	stress    = mol->anglestress;
	
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

		stressij = prefactor * (costheta - rij/rjk);
		stressik = prefactor * (rik*rik)/(rij*rjk);
		stressjk = prefactor * (costheta - rjk/rij);
		
		stress[atomi] = stress[atomi] - stressij - stressik;// negative sign to represent hydrostatic pressure
		stress[atomj] = stress[atomj] - stressij - stressjk;
		stress[atomk] = stress[atomk] - stressik - stressjk;
	}
}
