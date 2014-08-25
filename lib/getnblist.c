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

#define NBCUTOFF 10.0

/* generate the non-bonded list of the atoms using the brute-force approach
 * this also ensures that exclusions were taken into account
 */
void getnblist(struct topology *mol)
{	
	int i,j,k,found;
	float rij;
	
	/* This is the ugly way of looking at all pairs if they are within the cutoff
	 * and if they are not exclusions
	 */
	for(i=0;i<mol->numatoms;i++)
	{
		mol->nbcount[i] = 0;
		mol->ncount[i]  = 0;
		for(j=0;j<mol->numatoms;j++)
		{
			rij = distance(mol,i,j);
			if(rij <= NBCUTOFF)
			{
				found = 0;
				for(k=0;k<mol->exclusioncount[i];k++)
				{
					if(mol->exclusions[i][k] == j)
					{
						found = 1;
						break;
					}
				}
				if(found == 0)
				{
					mol->nblist[i][mol->nbcount[i]] = j;
					mol->nbcount[i] = mol->nbcount[i] + 1;
				}
				mol->nlist[i][mol->ncount[i]] = j;
				mol->ncount[i] = mol->ncount[i] + 1;
			}
		}
	}
}
