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


/* the following function computes the angle energy
 */
float angleenergy(struct topology *mol)
{
	int i;
	float theta, theta2, energy;

	energy = 0.0;
	for(i=0;i<(mol->numangles);i++)
	{
		/* angle energy is computed as k* (theta-thetao)^2
		 */
		theta  = mol->angles[i].angle - angle(mol,mol->angles[i].atomi,mol->angles[i].atomj,mol->angles[i].atomk);
		theta2 = theta*theta;
		energy = energy + mol->angles[i].fc * theta2;
	}
	return energy;
}
