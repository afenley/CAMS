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

/* The following function computes the bond stresses
 */
void bondstress(struct topology *mol)
{
	int i;
	float stressij,rij;
	int numbonds;
	struct bond *bonds;
	struct bond *bonds2;
	float *bondstress;
	int atomi,atomj;

	numbonds     = mol->numbonds;
	bonds        = mol->bonds;
	bondstress   = mol->bondstress;

	for(i=0;i<numbonds;i++)
	{
		/* bond stress is computed as 2*k*rij*(rij - ro)
		 */	
		bonds2 = &bonds[i];
		atomi = bonds2->atomi;
		atomj = bonds2->atomj;

		rij      = distance(mol,atomi,atomj);
		stressij = (bonds2->fc)*rij*(rij - bonds2->length);
		stressij = stressij + stressij;                           // replacing the 2.0* by a double sum

		bondstress[atomi] = bondstress[atomi] - stressij; // negative sign to represent hydrostatic pressure
		bondstress[atomj] = bondstress[atomj] - stressij;
	}
}
