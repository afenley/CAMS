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
#include"../include/structdefs.h"

/* computes angle between three points
 */
float cosangle(struct topology *mol, int i, int j, int k)
{
	float xij, yij, zij, xkj, ykj, zkj;
	float magij, magkj;
	float dotproduct;
	float (*crd)[3];

	crd = mol->crd;

	xij = crd[i][0] - crd[j][0];
	yij = crd[i][1] - crd[j][1];
	zij = crd[i][2] - crd[j][2];

	xkj = crd[k][0] - crd[j][0];
	ykj = crd[k][1] - crd[j][1];
	zkj = crd[k][2] - crd[j][2];
	
	magij = sqrt(xij*xij + yij*yij + zij*zij);
	magkj = sqrt(xkj*xkj + ykj*ykj + zkj*zkj);
	
	dotproduct = xij*xkj + yij*ykj + zij*zkj;
	return (dotproduct/(magij*magkj));
}
