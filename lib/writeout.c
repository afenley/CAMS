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
#include"../include/structdefs.h"

extern int splitflag;

/* write a pdb with the current stresses in the topology structure
 */
void writeout(char outfilename[],struct topology *mol)
{
	char  filename[200];
	float totalstress;
	int i,j,indexgroupcount;
        float *bondstress,*anglestress,*dihedralstress,*vdwstress,*qqstress,*solventstress,*kineticstress,*avgstress,*msfstress;
        float *avgbond,*avgangle,*avgdihedral,*avgvdw,*avgqq,*avgsolvent,*avgkinetic;
	float *msfbond,*msfangle,*msfdihedral,*msfvdw,*msfqq,*msfsolvent,*msfkinetic;
	FILE  *bondfile,*anglefile,*dihedralfile,*vdwfile,*qqfile,*solventfile,*kineticfile,*nonbondfile,*bondedfile,*elecfile,*totalfile;
        int *indexgroup;

        indexgroup     = mol->indexgroup;
        avgstress      = mol->avgstress;
        msfstress      = mol->msfstress;
        bondstress     = mol->bondstress;
        anglestress    = mol->anglestress;
        dihedralstress = mol->dihedralstress;
        vdwstress      = mol->vdwstress;
        qqstress       = mol->qqstress;
        solventstress  = mol->solventstress;
        kineticstress  = mol->kineticstress;

	if(splitflag)
	{
		avgbond     = mol->avgbond;
		avgangle    = mol->avgangle;
		avgdihedral = mol->avgdihedral;
		avgvdw      = mol->avgvdw;
		avgqq       = mol->avgqq;
		avgsolvent  = mol->avgsolvent;
		avgkinetic  = mol->avgkinetic;

		msfbond     = mol->msfbond;
		msfangle    = mol->msfangle;
                msfdihedral = mol->msfdihedral;
                msfvdw      = mol->msfvdw;
                msfqq       = mol->msfqq;
                msfsolvent  = mol->msfsolvent;
                msfkinetic  = mol->msfkinetic;

		strcpy(filename,outfilename); strcat(filename,"_bond.dat");
		bondfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_angle.dat");
		anglefile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_dihedral.dat");
		dihedralfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_vdw.dat");
		vdwfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_qq.dat");
		qqfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_solvent.dat");
		solventfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_kinetic.dat");
		kineticfile  = fopen(filename,"a");
	
		strcpy(filename,outfilename); strcat(filename,"_bonded.dat");
		bondedfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_nonbond.dat");
		nonbondfile  = fopen(filename,"a");
		
		strcpy(filename,outfilename); strcat(filename,"_elec.dat");
		elecfile  = fopen(filename,"a");
	}
	
	strcpy(filename,outfilename); strcat(filename,"_total.dat");
	totalfile  = fopen(filename,"a");

	indexgroupcount = mol->indexgroupcount;	
	for(j=0;j<indexgroupcount;j++)
	{
		i = mol->indexgroup[j];
		totalstress  = bondstress[i]+anglestress[i]+dihedralstress[i]+vdwstress[i]+qqstress[i]+solventstress[i]+kineticstress[i];
		avgstress[j] = avgstress[j] + totalstress;
                msfstress[j] = msfstress[j] + totalstress*totalstress;
		fprintf(totalfile,"%8.6f ",totalstress);
		
		if(splitflag)
		{
/* Adding averages to splits - ATF */
			avgbond[j]     = avgbond[j] + bondstress[i];
			avgangle[j]    = avgangle[j] + anglestress[i];
			avgdihedral[j] = avgdihedral[j] + dihedralstress[i];
			avgvdw[j]      = avgvdw[j] + vdwstress[i];
			avgqq[j]       = avgqq[j] + qqstress[i];
			avgsolvent[j]  = avgsolvent[j] + solventstress[i];
                        avgkinetic[j]  = avgkinetic[j] + kineticstress[i];

                        msfbond[j]     = msfbond[j] + bondstress[i]*bondstress[i];
                        msfangle[j]    = msfangle[j] + anglestress[i]*anglestress[i];
                        msfdihedral[j] = msfdihedral[j] + dihedralstress[i]*dihedralstress[i];
                        msfvdw[j]      = msfvdw[j] + vdwstress[i]*vdwstress[i];
                        msfqq[j]       = msfqq[j] + qqstress[i]*qqstress[i];
                        msfsolvent[j]  = msfsolvent[j] + solventstress[i]*solventstress[i];
                        msfkinetic[j]  = msfkinetic[j] + kineticstress[i]*kineticstress[i];

			fprintf(bondfile,"%8.6f ",bondstress[i]);
			fprintf(anglefile,"%8.6f ",anglestress[i]);
			fprintf(dihedralfile,"%8.6f ",dihedralstress[i]);
			fprintf(vdwfile,"%8.6f ",vdwstress[i]);
			fprintf(qqfile,"%8.6f ",qqstress[i]);
			fprintf(solventfile,"%8.6f ",solventstress[i]);
			fprintf(kineticfile,"%8.6f ",kineticstress[i]);
			fprintf(bondedfile,"%8.6f ",bondstress[i]+anglestress[i]+dihedralstress[i]);
			fprintf(nonbondfile,"%8.6f ",vdwstress[i]+qqstress[i]+solventstress[i]);
			fprintf(elecfile,"%8.6f ",qqstress[i]+solventstress[i]);
        	}        
	}
	
	fprintf(totalfile,"\n");
	fclose(totalfile);
	if(splitflag)
	{
		fprintf(bondfile,"\n");
		fprintf(anglefile,"\n");
		fprintf(dihedralfile,"\n");
		fprintf(vdwfile,"\n");
		fprintf(qqfile,"\n");
		fprintf(solventfile,"\n");
		fprintf(kineticfile,"\n");
		fprintf(bondedfile,"\n");
		fprintf(nonbondfile,"\n");
		fprintf(elecfile,"\n");
		fclose(bondfile);
		fclose(anglefile);
		fclose(dihedralfile);
		fclose(vdwfile);
		fclose(qqfile);
		fclose(solventfile);
		fclose(kineticfile);
		fclose(bondedfile);
		fclose(nonbondfile);
		fclose(elecfile);
	}
}
