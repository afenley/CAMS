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
#include"../include/structdefs.h"

/* computes distance squared between two atoms
 */
float distance2(struct topology *mol,int i,int j)
{
	float dx, dy, dz;
	float *crd1, *crd2;

	crd1 = mol->crd[i];
	crd2 = mol->crd[j];

	dx = crd1[0] - crd2[0];
	dy = crd1[1] - crd2[1];
	dz = crd1[2] - crd2[2];

 	return (dx*dx + dy*dy + dz*dz);
}
