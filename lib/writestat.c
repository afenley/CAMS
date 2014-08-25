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
#include<string.h>
#include"../include/mdmath.h"
#include"../include/structdefs.h"

extern int splitflag;

/* write the pdb with means and variances
 */

void writestat(char outfilename[], struct topology *mol, int numframes)
{
	int   i,j;
	/* concurrent filestreams */
	FILE  *outfile,*outfile1,*outfile2;
	/* simultaneous file output of all the split data -ATF */
	FILE  *outbond1,*outangle1,*outdihedral1,*outvdw1,*outqq1,*outsolvent1,*outkinetic1;
	FILE  *outbond2,*outangle2,*outdihedral2,*outvdw2,*outqq2,*outsolvent2,*outkinetic2;
	float *avgstress;
	float *msfstress;
        /* Adding averages for splits: allocate array pointers -ATF */
	float *avgbond,*avgangle,*avgdihedral,*avgvdw,*avgqq,*avgsolvent,*avgkinetic;
	float *msfbond,*msfangle,*msfdihedral,*msfvdw,*msfqq,*msfsolvent,*msfkinetic;
	int   indexgroupcount;
	char  filename[200];
	float meanstress,fluctuation;
	/* mean values for residue averages for all the split data -ATF */
	float meanbond,meanangle,meandihedral,meanvdw,meanqq,meansolvent,meankinetic;
        float flucbond,flucangle,flucdihedral,flucvdw,flucqq,flucsolvent,fluckinetic;
	int   count,atomindex,resindex,prevresindex;

	avgstress       = mol->avgstress;
	msfstress       = mol->msfstress;
	indexgroupcount = mol->indexgroupcount;

	/* Adding averages from splits: importing arrays from "mol" struct -ATF */
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
	}
	/* computing the average stress of the entpre trajectory by dividing the number 
	   of frames
	*/
        for(i=0;i<indexgroupcount;i++)
        {
                avgstress[i] = avgstress[i]/numframes;
                msfstress[i] = (msfstress[i]/numframes) - (avgstress[i]*avgstress[i]);
		/* computing averages from splits: divide number of frames -ATF */
		if(splitflag)
		{
			avgbond[i]     = avgbond[i]/numframes;
			avgangle[i]    = avgangle[i]/numframes;
			avgdihedral[i] = avgdihedral[i]/numframes;
			avgvdw[i]      = avgvdw[i]/numframes;
			avgqq[i]       = avgqq[i]/numframes;
			avgsolvent[i]  = avgsolvent[i]/numframes;
			avgkinetic[i]  = avgkinetic[i]/numframes;
                        /* computing mean squared fluctuations from splits -ATF */
			msfbond[i]     = msfbond[i]/numframes - (avgbond[i]*avgbond[i]);
                        msfangle[i]    = msfangle[i]/numframes - (avgangle[i]*avgangle[i]);
                        msfdihedral[i] = msfdihedral[i]/numframes - (avgdihedral[i]*avgdihedral[i]);
                        msfvdw[i]      = msfvdw[i]/numframes - (avgvdw[i]*avgvdw[i]);
                        msfqq[i]       = msfqq[i]/numframes - (avgqq[i]*avgqq[i]);
                        msfsolvent[i]  = msfsolvent[i]/numframes - (avgsolvent[i]*avgsolvent[i]);
                        msfkinetic[i]  = msfkinetic[i]/numframes - (avgkinetic[i]*avgkinetic[i]);
		}
        }

	/* writing atomic average and msf data */
	/* averages: */
	strcpy(filename,outfilename);strcat(filename,"_atomavg.pdb");
	outfile = fopen(filename,"w");
	fprintf(outfile,"REMARK    STRESS PROGRAM DEVELOPED BY HARI S. MUDDANA AND ANDREW T. FENLEY. PLEASE CITE.\n");
        for(i=0;i<indexgroupcount;i++)
        {
		fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgstress[i]);
        }
	fprintf(outfile,"END\n");
        fclose(outfile);

	/* mean square fluctuations */
        strcpy(filename,outfilename);strcat(filename,"_atommsf.pdb");
        outfile = fopen(filename,"w");
        fprintf(outfile,"REMARK    STRESS PROGRAM DEVELOPED BY HARI S. MUDDANA AND ANDREW T. FENLEY. PLEASE CITE.\n");
        for(i=0;i<indexgroupcount;i++)
        {
                fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfstress[i]);
        }
        fprintf(outfile,"END\n");
        fclose(outfile);

	/* writing atomic average data for split */
	if(splitflag)
	{
		/*bond*/
		strcpy(filename,outfilename);strcat(filename,"_bond_atomavg.pdb");
		outfile = fopen(filename,"w");
		fprintf(outfile,"REMARK    BOND STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
		for(i=0;i<indexgroupcount;i++)
		{
			fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgbond[i]);
		}
		fprintf(outfile,"END\n");
		fclose(outfile);

        	/*bond msf*/
        	strcpy(filename,outfilename);strcat(filename,"_bond_atommsf.pdb");
        	outfile = fopen(filename,"w");
        	fprintf(outfile,"REMARK    BOND STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
        	for(i=0;i<indexgroupcount;i++)
        	{
                	fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfbond[i]);
        	}
        	fprintf(outfile,"END\n");
        	fclose(outfile);

                /*angle*/
                strcpy(filename,outfilename);strcat(filename,"_angle_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    ANGLE STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgangle[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*angle msf*/
                strcpy(filename,outfilename);strcat(filename,"_angle_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    ANGLE STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfangle[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*dihedral*/
                strcpy(filename,outfilename);strcat(filename,"_dihedral_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    DIHEDRAL STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgdihedral[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*dihedral msf*/
                strcpy(filename,outfilename);strcat(filename,"_dihedral_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    DIHEDRAL STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfdihedral[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*vdw*/
                strcpy(filename,outfilename);strcat(filename,"_vdw_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    vdW STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgvdw[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*vdw msf*/
                strcpy(filename,outfilename);strcat(filename,"_vdw_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    vdW STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfvdw[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*qq*/
                strcpy(filename,outfilename);strcat(filename,"_qq_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    QQ STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgqq[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*qq msf*/
                strcpy(filename,outfilename);strcat(filename,"_qq_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    QQ STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfqq[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*solvent*/
                strcpy(filename,outfilename);strcat(filename,"_solvent_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    SOLVENT STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgsolvent[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*solvent msf*/
                strcpy(filename,outfilename);strcat(filename,"_solvent_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    SOLVENT STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfsolvent[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*kinetic*/
                strcpy(filename,outfilename);strcat(filename,"_kinetic_atomavg.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    KINETIC STRESS: ATOM AVG. PLEASE CITE JCC PAPER.\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,avgkinetic[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);

                /*kinetic msf*/
                strcpy(filename,outfilename);strcat(filename,"_kinetic_atommsf.pdb");
                outfile = fopen(filename,"w");
                fprintf(outfile,"REMARK    KINETIC STRESS: ATOM MSF. PLEASE CITE JCC PAPER\n");
                for(i=0;i<indexgroupcount;i++)
                {
                        fprintf(outfile,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",i+1,mol->name[i],mol->resname[i],1+mol->resnr[i],mol->tprcrd[i][0],mol->tprcrd[i][1],mol->tprcrd[i][2],1.0,msfkinetic[i]);
                }
                fprintf(outfile,"END\n");
                fclose(outfile);
		/* end of split average and msf atomic output */
	}

	/* writing residue average and msf stress
	*/
		
	strcpy(filename,outfilename);strcat(filename,"_resavg.pdb");
	outfile1 = fopen(filename,"w");
	fprintf(outfile1,"REMARK    STRESS PROGRAM DEVELOPED BY HARI S. MUDDANA AND ANDREW T. FENLEY. PLEASE CITE.\n");
	strcpy(filename,outfilename);strcat(filename,"_resmsf.pdb");
	outfile2 = fopen(filename,"w");
	fprintf(outfile2,"REMARK    STRESS PROGRAM DEVELOPED BY HARI S. MUDDANA AND ANDREW T. FENLEY. PLEASE CITE.\n");

        fluctuation = 0.0;
        meanstress = 0.0;
        count = 0;
        prevresindex = mol->resnr[mol->indexgroup[0]];
        for(i=0;i<indexgroupcount;i++)
        {
                atomindex = mol->indexgroup[i];
                resindex  = mol->resnr[atomindex];
                if(resindex == prevresindex)
                {
                        meanstress = meanstress + avgstress[i];
                        fluctuation = fluctuation + msfstress[i];
                        count = count + 1;
                }
                else
                {
                        for(j=i-count;j<i;j++)
                        {
		                fprintf(outfile1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanstress/count);
                		fprintf(outfile2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,fluctuation/count);
                        }
                        fluctuation = msfstress[i];
                        meanstress = avgstress[i];
                        count = 1;
                        prevresindex = resindex;
                }
        }
        for(j=i-count;j<i;j++)
        {
		fprintf(outfile1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanstress/count);
               	fprintf(outfile2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,fluctuation/count);
        }
	
	fprintf(outfile1,"END\n");
	fprintf(outfile2,"END\n");
        fclose(outfile1);
	fclose(outfile2);

	/* writing residue averages for the split data.
 	   due to the large amount of open file streams
	   we might want to consider functionalizing this.
	*/
	if(splitflag)
	{
		strcpy(filename,outfilename);strcat(filename,"_bond_resavg.pdb");
		outbond1 = fopen(filename,"w");
		fprintf(outbond1,"REMARK    BOND STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_angle_resavg.pdb");
                outangle1 = fopen(filename,"w");
                fprintf(outangle1,"REMARK    ANGLE STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_dihedral_resavg.pdb");
                outdihedral1 = fopen(filename,"w");
                fprintf(outdihedral1,"REMARK    DIHEDRAL STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_vdw_resavg.pdb");
                outvdw1 = fopen(filename,"w");
                fprintf(outvdw1,"REMARK    vdW STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_qq_resavg.pdb");
                outqq1 = fopen(filename,"w");
                fprintf(outqq1,"REMARK    QQ STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_solvent_resavg.pdb");
                outsolvent1 = fopen(filename,"w");
                fprintf(outsolvent1,"REMARK    SOLVENT STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_kinetic_resavg.pdb");
                outkinetic1 = fopen(filename,"w");
                fprintf(outkinetic1,"REMARK    KINETIC STRESS: RESIDUE AVG. PLEASE CITE JCC PAPER.\n");
		/* placeholder for MSF splits later: */
		strcpy(filename,outfilename);strcat(filename,"_bond_resmsf.pdb");
		outbond2 = fopen(filename,"w");
		fprintf(outbond2,"REMARK    BOND STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_angle_resmsf.pdb");
                outangle2 = fopen(filename,"w");
                fprintf(outangle2,"REMARK    ANGLE STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_dihedral_resmsf.pdb");
                outdihedral2 = fopen(filename,"w");
                fprintf(outdihedral2,"REMARK    DIHEDRAL STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_vdw_resmsf.pdb");
                outvdw2 = fopen(filename,"w");
                fprintf(outvdw2,"REMARK    vdW STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_qq_resmsf.pdb");
                outqq2 = fopen(filename,"w");
                fprintf(outqq2,"REMARK    QQ STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_solvent_resmsf.pdb");
                outsolvent2 = fopen(filename,"w");
                fprintf(outsolvent2,"REMARK    SOLVENT STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

                strcpy(filename,outfilename);strcat(filename,"_kinetic_resmsf.pdb");
                outkinetic2 = fopen(filename,"w");
                fprintf(outkinetic2,"REMARK    KINETIC STRESS: RESIDUE MSF. PLEASE CITE JCC PAPER.\n");

		meanbond     = 0.0;
		meanangle    = 0.0;
		meandihedral = 0.0;
		meanvdw      = 0.0;
		meanqq       = 0.0;
		meansolvent  = 0.0;
		meankinetic  = 0.0;
                flucbond     = 0.0;
                flucangle    = 0.0;
                flucdihedral = 0.0;
                flucvdw      = 0.0;
                flucqq       = 0.0;
                flucsolvent  = 0.0;
                fluckinetic  = 0.0;
		count = 0;
		prevresindex = mol->resnr[mol->indexgroup[0]];

		for(i=0;i<indexgroupcount;i++)
		{
			atomindex = mol->indexgroup[i];
			resindex  = mol->resnr[atomindex];
			if(resindex == prevresindex)
			{
				meanbond     = meanbond + avgbond[i];
				meanangle    = meanangle + avgangle[i];
				meandihedral = meandihedral + avgdihedral[i];
				meanvdw      = meanvdw + avgvdw[i];
				meanqq       = meanqq + avgqq[i];
				meansolvent  = meansolvent + avgsolvent[i];
				meankinetic  = meankinetic + avgkinetic[i];
				flucbond     = flucbond + msfbond[i];
				flucangle    = flucangle + msfangle[i];
				flucdihedral = flucdihedral + msfdihedral[i];
				flucvdw      = flucvdw + msfvdw[i];
				flucqq       = flucqq + msfqq[i];
				flucsolvent  = flucsolvent + msfsolvent[i];
				fluckinetic  = fluckinetic + msfkinetic[i];
				count = count + 1;
			}
			else
			{
				for(j=i-count;j<i;j++)
				{
					fprintf(outbond1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanbond/count);
                                        fprintf(outangle1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanangle/count);
                                        fprintf(outdihedral1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meandihedral/count);
                                        fprintf(outvdw1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanvdw/count);
                                        fprintf(outqq1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanqq/count);
                                        fprintf(outsolvent1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meansolvent/count);
                                        fprintf(outkinetic1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meankinetic/count);

					fprintf(outbond2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucbond/count);
                                        fprintf(outangle2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucangle/count);
                                        fprintf(outdihedral2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucdihedral/count);
                                        fprintf(outvdw2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucvdw/count);
                                        fprintf(outqq2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucqq/count);
                                        fprintf(outsolvent2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucsolvent/count);
                                        fprintf(outkinetic2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,fluckinetic/count);
				}

				meanbond     = avgbond[i];
				meanangle    = avgangle[i];
				meandihedral = avgdihedral[i];
				meanvdw      = avgvdw[i];
				meanqq       = avgqq[i];
				meansolvent  = avgsolvent[i];
				meankinetic  = avgkinetic[i];
				flucbond     = msfbond[i];
				flucangle    = msfangle[i];
				flucdihedral = msfdihedral[i];
				flucvdw      = msfvdw[i];
				flucqq       = msfqq[i];
				flucsolvent  = msfsolvent[i];
				fluckinetic  = msfkinetic[i];
				count = 1;
				prevresindex = resindex;
			}
		}
		for(j=i-count;j<i;j++)
		{
			fprintf(outbond1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanbond/count);
                        fprintf(outangle1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanangle/count);
                        fprintf(outdihedral1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meandihedral/count);
                        fprintf(outvdw1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanvdw/count);
                        fprintf(outqq1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meanqq/count);
                        fprintf(outsolvent1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meansolvent/count);
                        fprintf(outkinetic1,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,meankinetic/count);
			fprintf(outbond2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucbond/count);
                        fprintf(outangle2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucangle/count);
                        fprintf(outdihedral2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucdihedral/count);
                        fprintf(outvdw2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucvdw/count);
                        fprintf(outqq2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucqq/count);
                        fprintf(outsolvent2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,flucsolvent/count);
                        fprintf(outkinetic2,"ATOM  %5d %4s %-4sX%4d    %8.3f%8.3f%8.3f%6.2f %9.6f  \n",j+1,mol->name[j],mol->resname[j],1+mol->resnr[j],mol->tprcrd[j][0],mol->tprcrd[j][1],mol->tprcrd[j][2],1.0,fluckinetic/count);
		}

		fprintf(outbond1,"END\n");
                fprintf(outangle1,"END\n");
                fprintf(outdihedral1,"END\n");
                fprintf(outvdw1,"END\n");
                fprintf(outqq1,"END\n");
                fprintf(outsolvent1,"END\n");
                fprintf(outkinetic1,"END\n");
		fprintf(outbond2,"END\n");
                fprintf(outangle2,"END\n");
                fprintf(outdihedral2,"END\n");
                fprintf(outvdw2,"END\n");
                fprintf(outqq2,"END\n");
                fprintf(outsolvent2,"END\n");
                fprintf(outkinetic2,"END\n");
		fclose(outbond1);
                fclose(outangle1);
                fclose(outdihedral1);
                fclose(outvdw1);
                fclose(outqq1);
                fclose(outsolvent1);
                fclose(outkinetic1);
                fclose(outbond2);
                fclose(outangle2);
                fclose(outdihedral2);
                fclose(outvdw2);
                fclose(outqq2);
                fclose(outsolvent2);
                fclose(outkinetic2);
	}
}
