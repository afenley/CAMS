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

#define NBCUTOFF 10.0
#define I_NBCUTOFF 0.10

/* generate the non-bonded list of the atoms using the gromacs grid-based method
 * this also ensures that exclusions were taken into account
 */
void getnblist_grid(struct topology *mol)
{	
	int i,j,k,p,m,n;
	int count;
	float x,y,z,xmin,xmax,ymin,ymax,zmin,zmax, xsize,ysize,zsize,dx,dy,dz;
	int numxbins, numybins, numzbins, numcells, xbin, ybin, zbin, gridindex;
	int **grid;
	int *gridcount;
	int *tempgrid1;
	int *tempgrid2;
	int cell,neighborcell;

	int numatoms;
	int *nbcount;
	int *ncount;
	float (*crd)[3];
	int **nblist;
	int **nlist;
	float r_cutoff;
	int *exclusioncount;
	int **exclusions;
	int currentcount;
	int found;
	float rij2;
	int atomi,atomj;
	int count1, count2;
	int newxbin,newybin,newzbin;

	numatoms       = mol->numatoms;
	nbcount        = mol->nbcount;
	ncount         = mol->ncount;
	crd            = mol->crd;
	nblist         = mol->nblist;
	nlist	       = mol->nlist;
	exclusioncount = mol->exclusioncount;
	exclusions     = mol->exclusions;
	r_cutoff       = NBCUTOFF*NBCUTOFF;
	
	/* all the counts and nblists were set to zeros
	 */
	for(i=0;i<numatoms;i++)
	{
		nbcount[i] = 0;
		ncount[i]  = 0;
	}
	
	/* compute the size of the box
	 */
	xmin = 10000.0; xmax = -10000.0; ymin = 10000.0; ymax = -10000.0; zmin = 10000.0; zmax = -10000.0;
	for(i=0;i<numatoms;i++)
	{
		x = crd[i][0];
		y = crd[i][1];
		z = crd[i][2];

		if(x < xmin) { xmin = x;}
		else if(x > xmax) { xmax = x;}
		if(y < ymin) { ymin = y;}
		else if(y > ymax) { ymax = y;}
		if(z < zmin) { zmin = z;}
		else if(z > zmax) { zmax = z;}
	}
	xmin = xmin - 5.0; xmax = xmax + 5.0;
	ymin = ymin - 5.0; ymax = ymax + 5.0;
	zmin = zmin - 5.0; zmax = zmax + 5.0;

	xsize = xmax-xmin;
	ysize = ymax-ymin;
	zsize = zmax-zmin;

	numxbins = int(floor((xsize+xsize)*I_NBCUTOFF));
	numybins = int(floor((ysize+ysize)*I_NBCUTOFF));
	numzbins = int(floor((zsize+zsize)*I_NBCUTOFF));

	if(numxbins == 0) numxbins = 1;
	if(numybins == 0) numybins = 1;
	if(numzbins == 0) numzbins = 1;
		
	dx = xsize/numxbins;
	dy = ysize/numybins;
	dz = zsize/numzbins;

	numcells = (numxbins)*(numybins)*(numzbins);

	grid = (int**) malloc(numcells*sizeof(int*));
	gridcount = (int*) malloc(numcells*sizeof(int));

	for(i=0;i<numcells;i++)
	{
		grid[i] = (int*) malloc(1000*sizeof(int)); //TYPEDEF 1000 atoms later: 5x5A^2
		gridcount[i] = 0;
	}

	for(i=0;i<numatoms;i++)
	{
		xbin = int(floor((crd[i][0] - xmin)/dx));
		ybin = int(floor((crd[i][1] - ymin)/dy));
		zbin = int(floor((crd[i][2] - zmin)/dz));

		gridindex = zbin*numxbins*numybins + ybin*numxbins + xbin;
		grid[gridindex][gridcount[gridindex]] = i;
		gridcount[gridindex] = gridcount[gridindex] + 1;
	}

	for(xbin=0;xbin<numxbins;xbin++)
	{
		for(ybin=0;ybin<numybins;ybin++)
		{
			for(zbin=0;zbin<numzbins;zbin++)
			{
				cell = zbin*numxbins*numybins + ybin*numxbins + xbin;
				for(i=-2;i<=2;i++)
				{
					for(j=-2;j<=2;j++)
					{
						for(k=-2;k<=2;k++)
						{
							newxbin = xbin + i;
							newybin = ybin + j;
							newzbin = zbin + k;
							if(newxbin < 0 || newxbin>= numxbins || newybin < 0 || newybin >=numybins || newzbin < 0 || newzbin >= numzbins)
							{
								continue;
							}
							else
							{
								neighborcell = newzbin*numxbins*numybins + newybin*numxbins + newxbin;
								count1 = gridcount[cell];
								count2 = gridcount[neighborcell];
								tempgrid1 = grid[cell];
								tempgrid2 = grid[neighborcell];
		
								for(m=0;m<count1;m++)
								{
									for(n=0;n<count2;n++)
									{
										atomi = tempgrid1[m];
										atomj = tempgrid2[n];
		
										rij2 = distance2(mol,atomi,atomj);
										if(rij2 <= r_cutoff)
										{
											found = 0;
											currentcount = exclusioncount[atomi];
		
											for(p=0;p<currentcount;p++)
											{
												if(exclusions[atomi][p] == atomj)
												{
													found = 1;
													break;
												}
											}
											/* Band-aid to a situation where an atom is
  											   not on its own exclusion list (last atom usually).
 											*/
											if(atomi == atomj && found == 0)
											{
												/*printf("%d  %d\n", atomi, atomj);*/
												found = 1;
											}
											if(found == 0)
											{
												nblist[atomi][nbcount[atomi]] = atomj;
												nbcount[atomi] = nbcount[atomi] + 1;
											}
											nlist[atomi][ncount[atomi]] = atomj;
											ncount[atomi] = ncount[atomi] + 1;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	for(i=0;i<numcells;i++)
	{
		free(grid[i]);
	}
	free(grid);
	free(gridcount);
}
