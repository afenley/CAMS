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
#include<string.h>
#include<stdlib.h>

extern char ndxfilename[];
extern char tprfilename[];
extern char trrfilename[];
extern char outfilename[];
extern int  debugflag;
extern int  gbflag;
extern int  gbsalt;
extern int  bondiflag;
extern int  kineticflag;
extern int  splitflag;

/* The following function reads the command line arguments from the user
 */
void commandlineargs(int argc, char *argv[])
{
	int i;
	FILE* infile;
	char filename[200];

	strcpy(ndxfilename,"");
	strcpy(trrfilename,"");
	strcpy(tprfilename,"");
	strcpy(outfilename,"output");
	gbflag    = 0;
	debugflag = 0;
	gbsalt    = 0.145;

	if(argc < 2)
	{
		printf("ERROR: no input arguments were found at the input command\n");
		printf("       please double check your command line arguments\n");
		printf("       e.g.: stress.exe -n test.ndx -s test.tpr -f test.trr -o output\n");
		exit(1);
	}

	i = 1;
	while(i<argc)
	{
		if(argv[i][0] != '-')
		{
			printf("ERROR: cannot understand the command at input\n");
			printf("       please double check your command line arguments\n");
			printf("       e.g.: stress.exe -n test.ndx -s test.tpr -f test.trr -o output\n");
			exit(1);
		}
		else
		{							
			if(!strcmp(argv[i],"-n"))
			{
				strcpy(ndxfilename,argv[i+1]); i++;
			}
			else if(!strcmp(argv[i],"-s"))
			{
				strcpy(tprfilename,argv[i+1]); i++;
			}
			else if(!strcmp(argv[i],"-f"))
			{
				strcpy(trrfilename,argv[i+1]); i++;
			}	
			else if(!strcmp(argv[i],"-o"))
			{
				strcpy(outfilename,argv[i+1]); i++;
			}	
			else if(!strcmp(argv[i],"-debug"))
			{
				debugflag = 1;
			}	
			else if(!strcmp(argv[i],"-split"))
			{
				splitflag = 1;
			}	
			else if(!strcmp(argv[i],"-kinetic"))
			{
				kineticflag = 1;
			}	
			else if(!strcmp(argv[i],"-bondi"))
			{
				bondiflag = 1;
			}	
			else if(!strcmp(argv[i],"-gb"))
			{
				gbflag = 1;	
			}	
			else if(!strcmp(argv[i],"-salt"))
			{
				gbsalt = atof(argv[i+1]); i++;
			}	
			else if(!strcmp(argv[i],"-h"))
			{
				exit(1);
			}	
			i++;	
		}

	}
	if ((strlen(outfilename) == 0) || (strlen(ndxfilename) == 0) || (strlen(trrfilename) == 0) || (strlen(tprfilename) == 0))
	{
		printf("ERROR: one of the required filename was not input at the command line.\n");
		printf("       e.g.: stress.exe -n test.ndx -s test.tpr -f test.trr -o output\n");
		exit(1);
	}
	if(strcmp(trrfilename+strlen(trrfilename)-4,".xtc") && strcmp(trrfilename+strlen(trrfilename)-4,".trr"))
	{
		printf("ERROR: The input trajectory file is neither a .trr nor a .xtc file.\n");
		printf("       Usage: stress.exe -n test.ndx -s test.tpr -f test.trr -o output\n");
		exit(1);
	}

	strcpy(filename,outfilename); strcat(filename,"_total.dat");
	infile = fopen(filename,"r");
	if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	
	strcpy(filename,outfilename); strcat(filename,"_atomavg.pdb");
	infile = fopen(filename,"r");
	if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	
	strcpy(filename,outfilename); strcat(filename,"_atommsf.pdb");
	infile = fopen(filename,"r");
	if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	
	strcpy(filename,outfilename); strcat(filename,"_resavg.pdb");
	infile = fopen(filename,"r");
	if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	
	strcpy(filename,outfilename); strcat(filename,"_resmsf.pdb");
	infile = fopen(filename,"r");
	if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	
	if(splitflag)
	{
		strcpy(filename,outfilename); strcat(filename,"_bond.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);} 
	
		strcpy(filename,outfilename); strcat(filename,"_angle.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
 		
		strcpy(filename,outfilename); strcat(filename,"_dihedral.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_vdw.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_qq.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_solvent.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_kinetic.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_bonded.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}

		strcpy(filename,outfilename); strcat(filename,"_nonbonded.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
		
		strcpy(filename,outfilename); strcat(filename,"_elec.dat");
		infile = fopen(filename,"r");
		if(infile != NULL) { printf("ERROR: %s file found in the current directory.\n Delete/move this file or give different name for -output option\n",filename); fclose(infile); exit(1);}
	}
}

	
