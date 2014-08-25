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
#include"../include/mdmath.h"
#include"../include/structdefs.h"


/* NOTE: energy calculations of bonded interactions rely on the assumption
 * that the molecule is whole (pbc condition). Periodic images are not taken 
 * into account in computing the energies.
 */


/* the following function computes the pair energy
 */
float pairenergy(struct topology *mol)
{
	int i,atomi,atomj;
	float energy,rij,rij2,rij6,rij12;
	float c6, c12;

	energy = 0.0;
	for(i=0;i<(mol->numpairs);i++)
	{
		atomi = mol->pairs[i].atomi;
		atomj = mol->pairs[i].atomj;
		/* pair lj and qq energy is computed using LJ and Columbic potential
		 * C6 and C12 coefficient of the pairs already include the fudgeLJ
		 * FudgeQQ is used for the electrostatic calculations
		 */
		rij = distance(mol,atomi,atomj);
		rij2 = rij*rij;
		rij6 = rij2*rij2*rij2;
		rij12 = rij6*rij6;
		c6 = (mol->pairs[i].c6)/6.0;
		c12 = (mol->pairs[i].c12)/12.0;
		energy = energy + (c12/rij12) - (c6/rij6);
	        energy = energy + (mol->fudgeQQ)*332.06*(mol->charge[atomi])*(mol->charge[atomj])/rij;
	}
	return energy;
}
