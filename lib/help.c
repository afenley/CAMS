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

/* The followinig function prints the help menu for the stress program
 */
void help()
{
	printf("\t\t:-)  Stress analysis program: v0.1 (-:\t\t\t\n");
	printf("\t\tby Hari S. Muddana and Andrew T. Fenley\n");
	printf("DESCRIPTION:\n");
	printf("------------\n");
	printf("The stress program reads a run input file (.tpr), \n");
	printf("an index file (.ndx), a trajectory file (.trr or .xtc), and \n");
	printf("outputs the stresses in a pdb trajectory (.pdb) format. All \n");
	printf("files are mandatory. The program assumes that you have gromacs \n");
	printf("properly installed and running. The code internally calls some \n");
	printf("gromacs commands (e.g. gmxdump).\n\n");
	printf("Option       Filename       Type           Description\n");
	printf("------------------------------------------------------\n");


	printf("  -s         topol.tpr      Input, (req.)   Run input file\n");
	printf("  -n         index.ndx      Input, (req.)   Index group for output.\n");
        printf("  -f         traj.trr       Input, (req.)   Input trajectory file. Can also be .xtc file\n");
        printf("  -o         output         Output,(Opt.)   Prefix for all output files, default = output\n");
	printf("  -debug                           (Opt.)   Flag to turn on debug information\n");
	printf("  -gb                              (Opt.)   Flag to turn on GB (by default gb is turned off)\n");
	printf("  -salt      0.145                 (Opt.)   Salt concentration in M for GB stresses\n");
	printf("  -bondi                           (Opt.)   Forces the GB code to use Bondi Radii\n");
	printf("  -kinetic                         (Opt.)   Includes kinetic stresses from velocities\n");
	printf("  -split                           (Opt.)   Flag to write stress contributions (bond, angle, dihedral, etc.)\n");
	printf("					    to seperate files.\n");
	printf("\n\nNOTES:\n");
	printf("------\n");
	//printf("1) when no index file is provided the stresses will be output for \n");
	//printf("   for all the atoms.\n");
	printf("1) Periodic boundary conditions are not taken into account at this \n");
	printf("   time. So, make sure that group of atoms (protein, etc.) you are \n");
	printf("   interested in is centered in the box. You might also want to remove \n");
	printf("   translation and rotation of the atom group. \n");
	printf("2) (to be ignored for now!)The user has to provide one of the following: a) a single pdb input \n");
	printf("   file, or b) a trajectory of the system. In either case, the number of \n");
	printf("   atoms in the pdb/trr file should match the number of atoms in your tpr \n");
	printf("   file. In case, where both .trr and .pdb are provided, the stresses will \n");
	printf("   be computed for the trajectory and not the pdb file. \n");
	printf("4) This program is not efficient in reading gromacs file formats. In most \n");
	printf("   cases, we dump the gromacs file into an ascii file and parse this file for \n");
	printf("   for input. Large temporary files will be created in this process, and so \n");
	printf("   make sure there is enough space on the disk.\n");
	printf("5) In case of any already existing output file with the same filename, the program\n");
	printf("   terminates without any output.\n");
	printf("6) Do not turn on the GB flag if you have explicit waters in the simulation trajectory\n\n\n");

}
