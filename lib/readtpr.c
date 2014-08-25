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
#include<string.h>
#include<math.h>
#include "../include/xdrfile.h"
#include "../include/xdrfile_trr.h"
#include "../include/xdrfile_xtc.h"
#include"../include/gmxlib.h"
#include"../include/structdefs.h"

extern int debugflag;

/* The following function all the grunt work of reading the different system
 * and force-field terms
 */
void readtpr(char tprfilename[],struct topology *mol)
{
	FILE *infile;
	FILE *infilecopy;
	char command[1000];
	char buffer[1000];
	float mass,charge,radius,s_hct;
	char name[10],tempstring[1000];
	char *temp;
	int i,n,count,numwords,numresidues,numatomtypes;
	fpos_t fposition;

	/* generate the ascii formatted topology file of the entire system
	 * using the gmxdump program
	 */
	infile = fopen(tprfilename,"r");
	if(infile == NULL) {printf("%s file not found\n",tprfilename); exit(1);}
	fclose(infile);

	sprintf(command,"(gmxdump -s %s -sys > topol_dump.top) > /dev/null 2>&1",tprfilename);
	system(command);

	/* In the following we will read the different sections of the .tpr file
	 */
	infile = fopen("topol_dump.top","r");
	while((fgets(buffer,1000,infile))!=NULL)
	{
		/* reading the name of the system. This is purely for debug purposes. We can 
		 * safely delete the following if loop */
		if(!strncmp(buffer,"topology:",9))
		{
			fgets(buffer,1000,infile);
			for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '=') { buffer[i] = ' ';} }
			temp = strstr(buffer,"name");
		}
		/* reading the atoms section of the .tpr dump. The following information is 
		 * gathered in this section: 1) numatoms, 2) mass, 3) charge, 4) resnr, 
		 * 5) atom name.
		 */
		if(!strncmp(buffer,"   atoms:",9))
		{
			fgets(buffer,1000,infile);
			for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '(' || buffer[i] == ')') { buffer[i] = ' ';} }
			temp = strstr(buffer,"atom");
			sscanf(temp+4,"%d",&mol->numatoms);
			mol->initialize();
			for(n=0;n<mol->numatoms;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&mol->ljtypes[n]);
				temp = strstr(buffer," m ");
				sscanf(temp+3,"%f",&mol->mass[n]);
				temp = strstr(buffer," q ");
				sscanf(temp+3,"%f",&mol->charge[n]);
				temp = strstr(buffer,"resind");
				sscanf(temp+6,"%d",&mol->resnr[n]);
			}
			fgets(buffer,1000,infile);
			for(n=0;n<mol->numatoms;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"name");
				sscanf(temp+4,"%s",mol->name[n]);
			}
		}
		/* reading the force-field atom type information. This section is also optional.
		 * This is likely not to be used in any piece of the code. Maybe in the PDB
		 * output
		 */
		if((!strncmp(buffer,"      type",10)) && (buffer[strlen(buffer)-2] == ':'))
		{
			for(n=0;n<mol->numatoms;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"name");
				sscanf(temp+4,"%s",mol->type[n]);
			}
		}
		/* reading the residue names of each atoms. This information is useful
		 * in generating the output PDB file
		 */
		if((!strncmp(buffer,"      residue",13)) && (buffer[strlen(buffer)-2] == ':'))
		{
			for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '(' || buffer[i] == ')') { buffer[i] = ' ';} }
			temp = strstr(buffer,"residue");
			sscanf(temp+7,"%d",&numresidues);

			for(n=0;n<numresidues;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"name");
				sscanf(temp+4,"%s",name);
				for(i=0;i<mol->numatoms;i++)
				{
					if(mol->resnr[i] == n)
					{
						strcpy(mol->resname[i],name);
					}
				}
			}
		}
		
		/*reading atomtype radius information for the GB code
		*/
		if((!strncmp(buffer,"   atomtypes:",13)))
		{
			fgetpos(infile,&fposition);
			numatomtypes = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"atomtype")!=NULL)
			{
				numatomtypes = numatomtypes + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the bond section
			 */
			fsetpos(infile,&fposition);
			for(n=0;n<numatomtypes;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' ||
                                                buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"gb_radius");
				sscanf(temp+9,"%f",&radius);
				temp = strstr(buffer,"S_hct");
				sscanf(temp+5,"%f",&s_hct);
				for(i=0;i<mol->numatoms;i++)
				{
					if(mol->ljtypes[i] == n)
					{
						mol->radius[i] = radius*10.0;
						mol->s_hct[i]  = s_hct;
					}
				}	
			}
		}

		if((!strncmp(buffer,"   excls",8)) && (buffer[strlen(buffer)-2] == ':'))
		{
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			count = -1;
			while(count < mol->numatoms-1)
			{
				fgets(buffer,1000,infile);
				if(!strncmp(buffer,"      excls",11))
				{
					count++;
					for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '"') { buffer[i] = ' ';} }
					temp = strstr(buffer,"=") + 1;
					numwords = wordcount(temp);
					for(i=0;i<numwords;i++)
					{
						sscanf(temp,"%s",tempstring);
						temp = strstr(temp,tempstring) + strlen(tempstring);
					       	mol->exclusions[count][mol->exclusioncount[count]] = atoi(tempstring);
						mol->exclusioncount[count] = mol->exclusioncount[count] + 1;	
					}
				}
				else
				{
					for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
						buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '"') { buffer[i] = ' ';} }
					temp = buffer;
					numwords = wordcount(temp);
					for(i=0;i<numwords;i++)
					{
						sscanf(temp,"%s",tempstring);
						temp = strstr(temp,tempstring) + strlen(tempstring);
					       	mol->exclusions[count][mol->exclusioncount[count]] = atoi(tempstring);
						mol->exclusioncount[count] = mol->exclusioncount[count] + 1;	
					}
				}

			}
		}
		if(!strncmp(buffer,"   idef:",8))
		{
			fgets(buffer,1000,infile);
			temp = strstr(buffer,"atnr");
			sscanf(temp+5,"%d",&mol->numljtypes);
			fgets(buffer,1000,infile);
			for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
				buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
			temp = strstr(buffer,"ntypes") + 6;
			sscanf(temp,"%d",&mol->numffterms);
			mol->initializeff();
			count = -1;
			while(count < mol->numffterms-1)
			{
				fgets(buffer,1000,infile);
				if(!strncmp(buffer,"         functype",17))
				{
					count++;
					for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || 
						buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
					if(strstr(buffer,"LJ_SR")!=NULL)
					{
						temp = strstr(buffer,"c6") + 2;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						mol->forcefield[count][0] = mol->forcefield[count][0] * (6e6/4.184);
						temp = strstr(buffer,"c12") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
						mol->forcefield[count][1] = mol->forcefield[count][1] * (12e12/4.184);
					}
					else if(strstr(buffer,"BONDS")!=NULL)
					{
						temp = strstr(buffer,"b0A") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						temp = strstr(buffer,"cbA") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
					}
					else if(strstr(buffer,"ANGLES")!=NULL)
					{
						temp = strstr(buffer,"thA") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						temp = strstr(buffer,"ctA") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
					}
					else if(strstr(buffer,"PDIHS")!=NULL)
					{
						temp = strstr(buffer,"phiA") + 4;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						temp = strstr(buffer,"cpA") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
						temp = strstr(buffer,"mult") + 4;
						sscanf(temp,"%f",&mol->forcefield[count][2]);
					}
					else if(strstr(buffer,"PIDIHS")!=NULL)
					{
						temp = strstr(buffer,"phiA") + 4;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						temp = strstr(buffer,"cpA") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
						temp = strstr(buffer,"mult") + 4;
						sscanf(temp,"%f",&mol->forcefield[count][2]);
					}
					else if(strstr(buffer,"RBDIHS")!=NULL)
					{
						temp = strstr(buffer,"rbcA[0]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						temp = strstr(buffer,"rbcA[1]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
						temp = strstr(buffer,"rbcA[2]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][2]);
						temp = strstr(buffer,"rbcA[3]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][3]);
						temp = strstr(buffer,"rbcA[4]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][4]);
						temp = strstr(buffer,"rbcA[5]") + 7;
						sscanf(temp,"%f",&mol->forcefield[count][5]);
					}
					else if(strstr(buffer,"LJ14")!=NULL)
					{
						temp = strstr(buffer,"c6A") + 3;
						sscanf(temp,"%f",&mol->forcefield[count][0]);
						mol->forcefield[count][0] = mol->forcefield[count][0] * (6e6/4.184);
						temp = strstr(buffer,"c12A") + 4;
						sscanf(temp,"%f",&mol->forcefield[count][1]);
						mol->forcefield[count][1] = mol->forcefield[count][1] * (12e12/4.184);
					}

				}
			}
		}

		/* read the fudge factor for 1-4 electrostatic interactions
		 */
		if(!strncmp(buffer,"      fudgeQQ",13))
		{
			for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
				buffer[i] == ']' || buffer[i] == ',' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
			temp = strstr(buffer,"fudgeQQ") + 7;
			sscanf(temp,"%f",&mol->fudgeQQ);
		}
		
		/* read the bond section of the topology. This is now the bonds in the system and not the 
		 * force-field. Reading the bond section is slightly ugly, check out the comments below
		 * Note that we modified the units in the following. Also, the bond function is modelled
		 * as E = k (r-r)^2
		 */
		if(!strncmp(buffer,"      Bond:",11))
		{
			/* the nr right after the bond section does not contain the right number of bonds, 
			 * so we first parse the file to see how many bonds are there in teh bond section
			 * and then come back to read the bond data. This is required to do dynamic memory 
			 * allocation*/
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numbonds = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"BONDS")!=NULL)
			{
				mol->numbonds = mol->numbonds + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the bond section
			 */
			fsetpos(infile,&fposition);
			mol->initializebonds();	
			for(n=0;n<mol->numbonds;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->bonds[n]).type));
				temp = strstr(buffer,"BONDS") + 5;
				sscanf(temp,"%d %d",&(mol->bonds[n].atomi),&(mol->bonds[n].atomj));
				/* In the following we are converting the bond lengths and force constants to
				 * angstroms and kcal/mol.A2 units. The standard units in gromacs are nanometers
				 * and KJ/mol.nm2.*/
				mol->bonds[n].length = mol->forcefield[mol->bonds[n].type][0]*10.0;
				mol->bonds[n].fc = mol->forcefield[mol->bonds[n].type][1]/836.8;
			}

		}

		/* read the angle section of the topology. Similar to the bond section, reading the angle section 
		 * is also slightly ugly. see the comments below. Note that the units are modified. Also note
		 * that the harmonic angle function used here is: E = k(theta - theta0)^2 
		 */	
		if(!strncmp(buffer,"      Angle:",12))
		{
			/* the nr right after the angle section does not contain the right number of angles, 
			 * so we first parse the file to see how many angles are there in the angle section
			 * and then come back to read the angle data. This is required to do dynamic memory 
			 * allocation*/
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numangles = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"ANGLES")!=NULL)
			{
				mol->numangles = mol->numangles + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the angle section
			 */
			fsetpos(infile,&fposition);
			mol->initializeangles();	
			for(n=0;n<mol->numangles;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->angles[n]).type));
				temp = strstr(buffer,"ANGLES") + 6;
				sscanf(temp,"%d %d %d",&(mol->angles[n].atomi),&(mol->angles[n].atomj),&(mol->angles[n].atomk));
				/* In the following we are converting the angle and force constants to
				 * radians and kcal/mol.rad2 units.
				 */
				mol->angles[n].angle = mol->forcefield[mol->angles[n].type][0]*M_PI/180.0;
				mol->angles[n].fc = mol->forcefield[mol->angles[n].type][1]/8.368;
			}

		}
		/* read the proper dihedral section of the topology. Similar to the bond section, reading the dihedral section 
		 * is also slightly ugly. See below for comments. Note that the units are modified. The functional form used in
		 * converting the parameters is E = k (1+ cos(n*theta - phi))
		 */	
		if(!strncmp(buffer,"      Proper Dih.:",18))
		{
			/* the nr right after the dihedral section does not contain the right number of dihedrals, 
			 * so we first parse the file to see how many dihedrals are there in the dihedral section
			 * and then come back to read the dihedral data. This is required to do dynamic memory 
			 * allocation*/
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numpdihedrals = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"PDIHS")!=NULL)
			{
				mol->numpdihedrals = mol->numpdihedrals + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the proper dihedral section
			 */
			fsetpos(infile,&fposition);
			mol->initializepdihedrals();	
			for(n=0;n<mol->numpdihedrals;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->pdihedrals[n]).type));
				temp = strstr(buffer,"PDIHS") + 5;
				sscanf(temp,"%d %d %d %d",&(mol->pdihedrals[n].atomi),&(mol->pdihedrals[n].atomj),&(mol->pdihedrals[n].atomk),&(mol->pdihedrals[n].atoml));
				/* In the following we are converting the angle and force constants to
				 * radians and kcal/mol.rad2 units.
				 */
				mol->pdihedrals[n].mult = mol->forcefield[mol->pdihedrals[n].type][2];
				mol->pdihedrals[n].fc = mol->forcefield[mol->pdihedrals[n].type][1]/4.184;
				mol->pdihedrals[n].phi = mol->forcefield[mol->pdihedrals[n].type][0]*M_PI/180.0;
			}

		}

		/* read the RB dihedral section of the topology. Similar to the bond section, reading the dihedral section 
		 * is also slightly ugly. See below for comments. Note that the units are modified. The functional form used in
		 * converting the parameters is 
		 */	
		if(!strncmp(buffer,"      Ryckaert-Bell.:",21))
		{
			/* the nr right after the dihedral section does not contain the right number of dihedrals, 
			 * so we first parse the file to see how many dihedrals are there in the dihedral section
			 * and then come back to read the dihedral data. This is required to do dynamic memory 
			 * allocation*/
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numrbdihedrals = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"RBDIHS")!=NULL)
			{
				mol->numrbdihedrals = mol->numrbdihedrals + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the proper dihedral section
			 */
			fsetpos(infile,&fposition);
			mol->initializerbdihedrals();	
			for(n=0;n<mol->numrbdihedrals;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->rbdihedrals[n]).type));
				temp = strstr(buffer,"RBDIHS") + 6;
				sscanf(temp,"%d %d %d %d",&(mol->rbdihedrals[n].atomi),&(mol->rbdihedrals[n].atomj),&(mol->rbdihedrals[n].atomk),&(mol->rbdihedrals[n].atoml));
				/* In the following we are converting the angle and force constants to
				 * radians and kcal/mol.rad2 units.
				 */
				mol->rbdihedrals[n].C0 = mol->forcefield[mol->rbdihedrals[n].type][0]/4.184;
				mol->rbdihedrals[n].C1 = mol->forcefield[mol->rbdihedrals[n].type][1]/4.184;
				mol->rbdihedrals[n].C2 = mol->forcefield[mol->rbdihedrals[n].type][2]/4.184;
				mol->rbdihedrals[n].C3 = mol->forcefield[mol->rbdihedrals[n].type][3]/4.184;
				mol->rbdihedrals[n].C4 = mol->forcefield[mol->rbdihedrals[n].type][4]/4.184;
				mol->rbdihedrals[n].C5 = mol->forcefield[mol->rbdihedrals[n].type][5]/4.184;
				
				mol->rbdihedrals[n].k1 = -0.75*(mol->rbdihedrals[n].C3) - (mol->rbdihedrals[n].C1);
				mol->rbdihedrals[n].k2 = -0.5*(mol->rbdihedrals[n].C2 + mol->rbdihedrals[n].C4); 
				mol->rbdihedrals[n].k3 = -(mol->rbdihedrals[n].C3)/4.0;
				mol->rbdihedrals[n].k4 = -(mol->rbdihedrals[n].C4)/8.0;
				mol->rbdihedrals[n].k0 = mol->rbdihedrals[n].C0 - mol->rbdihedrals[n].k2 + mol->rbdihedrals[n].k4;
			}

		}
		/* read the improper dihedral section of the topology. Similar to the bond section, reading the dihedral section 
		 * is also slightly ugly. See below for comments. Note that the units are modified. The functional form used in
		 * converting the parameters is E = k (theta - theta0)^2
		 */	
		if(!strncmp(buffer,"      Improper Dih.:",20))
		{
			/* the nr right after the dihedral section does not contain the right number of dihedrals, 
			 * so we first parse the file to see how many dihedrals are there in the dihedral section
			 * and then come back to read the dihedral data. This is required to do dynamic memory 
			 * allocation*/
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numidihedrals = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"PIDIHS")!=NULL)
			{
				mol->numidihedrals = mol->numidihedrals + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the proper dihedral section
			 */
			fsetpos(infile,&fposition);
			mol->initializeidihedrals();	
			for(n=0;n<mol->numidihedrals;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->idihedrals[n]).type));
				temp = strstr(buffer,"PIDIHS") + 6;
				sscanf(temp,"%d %d %d %d",&(mol->idihedrals[n].atomi),&(mol->idihedrals[n].atomj),&(mol->idihedrals[n].atomk),&(mol->idihedrals[n].atoml));
				/* In the following we are converting the angle and force constants to
				 * radians and kcal/mol.rad2 units.
				 */
				mol->idihedrals[n].fc = mol->forcefield[mol->idihedrals[n].type][1]/4.184;
				mol->idihedrals[n].phi = mol->forcefield[mol->idihedrals[n].type][0]*M_PI/180.0;
			}

		}
		/* reading the pairs section of the topology.
		 */
		if(!strncmp(buffer,"      LJ-14:",12))
		{
			/* getting the count of the number of pairs
			 */
			fgets(buffer,1000,infile);
			fgets(buffer,1000,infile);
			fgetpos(infile,&fposition);
			mol->numpairs = 0;
			fgets(buffer,1000,infile);
			while(strstr(buffer,"LJ14")!=NULL)
			{
				mol->numpairs = mol->numpairs + 1;
				fgets(buffer,1000,infile);
			}
			/* resetting the file pointer position and reading the contents of the pairs section
			 */
			fsetpos(infile,&fposition);
			mol->initializepairs();	
			for(n=0;n<mol->numpairs;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '=' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"type") + 4;
				sscanf(temp,"%d",&((mol->pairs[n]).type));
				temp = strstr(buffer,"LJ14") + 4;
				sscanf(temp,"%d %d",&(mol->pairs[n].atomi),&(mol->pairs[n].atomj));
				/* In the following we are converting the c6 and c12 constants to
				 * angstroms and kcal/mol units.
				 */
				mol->pairs[n].c6 = mol->forcefield[mol->pairs[n].type][0];
				mol->pairs[n].c12 = mol->forcefield[mol->pairs[n].type][1];
			}

		}

		/* reading the starting coordinates in the topology file. This is just to fill the coordinates and will
		 * only the used when a trajectory file is not provided.
		 */
		if(buffer[0] == 'x' && buffer[strlen(buffer)-2] == ':')
		{
			/* reading the x, y, and z coordinates of the atoms. Note that the coordinates 
			 * are converted to angstroms units
			 */
			for(n=0;n<mol->numatoms;n++)
			{
				fgets(buffer,1000,infile);
				for(i=0;i<strlen(buffer);i++) { if(buffer[i] == '{' || buffer[i] == '}' || buffer[i] == '[' || 
					buffer[i] == ']' || buffer[i] == ',' || buffer[i] == ')' || buffer[i] == '(' || buffer[i] == '"') { buffer[i] = ' ';} }
				temp = strstr(buffer,"=") + 1;
				sscanf(temp,"%f %f %f", &mol->crd[n][0],&mol->crd[n][1],&mol->crd[n][2]);
				mol->crd[n][0]    = mol->crd[n][0]*10.0;
				mol->crd[n][1]    = mol->crd[n][1]*10.0;
				mol->crd[n][2]    = mol->crd[n][2]*10.0;
				mol->tprcrd[n][0] = mol->crd[n][0];
				mol->tprcrd[n][1] = mol->crd[n][1];
				mol->tprcrd[n][2] = mol->crd[n][2];
			}
		}
	

	}
	fclose(infile);
	
	/* deleting the temporary dump file generated by gmxdump
	 */
	system("rm topol_dump.top");
}
