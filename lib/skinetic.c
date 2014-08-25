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
#define I_JTOKCAL 0.239006

/* The following function computes the kinetic stress
 */
void kineticstress(struct topology *mol)
{
	int   i;
	int   numatoms;
	float *kineticstress;
	float (*velocity)[3];
	float *mass;

	numatoms        = mol->numatoms;
	kineticstress   = mol->kineticstress;
	velocity        = mol->velocity;
	mass	        = mol->mass;

	for(i=0;i<numatoms;i++)
	{
		/* kinetic stress is computed as m*v*v
		 */	
		kineticstress[i] = -1.0 * I_JTOKCAL*(mass[i] * ((velocity[i][0])*(velocity[i][0]) + (velocity[i][1])*(velocity[i][1]) + (velocity[i][2])*(velocity[i][2]))); // negative sign to represent hydrostatic pressure 
	}
}
