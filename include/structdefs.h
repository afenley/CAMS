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

/********************************************************************************
 * Structures for bonds, angles, dihedrals
 */
struct bond
{
	int atomi;
	int atomj;
	int type;
	float length;
	float fc;
};
struct angle
{
	int atomi;
	int atomj;
	int atomk;
	int type;
	float angle;
	float fc;
};
struct pair
{
	int atomi;
	int atomj;
	int type;
	float c6;
	float c12;
};
struct pdihedral
{
	int atomi;
	int atomj;
	int atomk;
	int atoml;
	int type;
	float fc;
	float phi;
	float mult;
};
struct rbdihedral
{
	int atomi;
	int atomj;
	int atomk;
	int atoml;
	int type;
	float C0;
	float C1;
	float C2;
	float C3;
	float C4;
	float C5;
	float k0;
	float k1;
	float k2;
	float k3;
	float k4;
};
struct idihedral
{
	int atomi;
	int atomj;
	int atomk;
	int atoml;
	int type;
	float fc;
	float phi;
};


/********************************************************************************
 */

struct topology
{

/* Definition of the different variables
 * numatoms         : Number of the atoms in the system
 * mass[]           : Array containing the masses of all atoms
 * charge[]         : Array containing the charges of all atoms
 * resnr[]          : residue number of all atoms
 * name[]           : Atom name
 * type[]           : force-field atom name
 * resname[]        : residue names of all atoms
 * exclusioncount[] : number of exclusions of each atom
 * exclusions[][]   : list of excluded atoms for each atom
 * fudeQQ	    : Fudge factor for 1-4 electrostatic interactions
 * numbonds         : number of bonds in the system
 * numangles        : number of angles in the system
 */
	int numatoms;
	float *mass;
	float *charge;
	float **chargeij;
	float *radius;
	float *s_hct;
	int  *resnr;
	char **name;
	char **type;
	char **resname;
	int *exclusioncount;
	int **exclusions;
	float fudgeQQ;
	float **forcefield;
	float (*crd)[3];
	float **tprcrd;
	float (*velocity)[3];
	int *indexgroup;
	int *ljtypes;
	int **nblist;
	int *nbcount;
	int **nlist;
	int *ncount;
	
	int numffterms;
	int numbonds;
	int numangles;
	int numpairs;
	int numpdihedrals;
	int numidihedrals;
	int numrbdihedrals;
	int indexgroupcount;
	int numljtypes;
		
	struct bond         *bonds;
	struct angle        *angles;
	struct pair         *pairs;
	struct rbdihedral   *rbdihedrals;
	struct pdihedral    *pdihedrals;
	struct idihedral    *idihedrals;
	
	float *bondstress;
	float *anglestress;
	float *dihedralstress;
	float *vdwstress;
	float *qqstress;
	float *solventstress;
	float *kineticstress;
	float *avgstress;
	float *msfstress;
        //float **force;
	
	float *avgbond;
	float *avgangle;
	float *avgdihedral;
	float *avgvdw;
	float *avgqq;
	float *avgsolvent;
	float *avgkinetic;
        float *msfbond;
        float *msfangle;
        float *msfdihedral;
        float *msfvdw;
        float *msfqq;
        float *msfsolvent;
        float *msfkinetic;
	
	float xbox;
	float ybox;
	float zbox;
	float volume;

	/* dynamic memory allocation of the index group
	 */
	void initializegroup()
	{
		int i;

		indexgroup     = (int*) malloc(indexgroupcount*sizeof(int));
		avgstress      = (float*) malloc(indexgroupcount*sizeof(float));
		msfstress      = (float*) malloc(indexgroupcount*sizeof(float));
		
		for(i=0;i<indexgroupcount;i++)
		{
			avgstress[i] = 0.0;
			msfstress[i] = 0.0;
		}
	}

	/* dynamic memory allocation of the split average stresses - ATF */
        void initializesplit()
        {
                int i;

                avgbond      = (float*) malloc(indexgroupcount*sizeof(float));
		avgangle     = (float*) malloc(indexgroupcount*sizeof(float));
		avgdihedral  = (float*) malloc(indexgroupcount*sizeof(float));
		avgvdw       = (float*) malloc(indexgroupcount*sizeof(float));
		avgqq        = (float*) malloc(indexgroupcount*sizeof(float));
		avgsolvent   = (float*) malloc(indexgroupcount*sizeof(float));
		avgkinetic   = (float*) malloc(indexgroupcount*sizeof(float));
                msfbond      = (float*) malloc(indexgroupcount*sizeof(float));
                msfangle     = (float*) malloc(indexgroupcount*sizeof(float));
                msfdihedral  = (float*) malloc(indexgroupcount*sizeof(float));
                msfvdw       = (float*) malloc(indexgroupcount*sizeof(float));
                msfqq        = (float*) malloc(indexgroupcount*sizeof(float));
                msfsolvent   = (float*) malloc(indexgroupcount*sizeof(float));
                msfkinetic   = (float*) malloc(indexgroupcount*sizeof(float));

                for(i=0;i<indexgroupcount;i++)
                {
                        avgbond[i]     = 0.0;
                        avgangle[i]    = 0.0;
			avgdihedral[i] = 0.0;
			avgvdw[i]      = 0.0;
			avgqq[i]       = 0.0;
			avgsolvent[i]  = 0.0;
			avgkinetic[i]  = 0.0;
                        msfbond[i]     = 0.0;
                        msfangle[i]    = 0.0;
                        msfdihedral[i] = 0.0;
                        msfvdw[i]      = 0.0;
                        msfqq[i]       = 0.0;
                        msfsolvent[i]  = 0.0;
                        msfkinetic[i]  = 0.0;
                }
        }
	
	/* initialize all other variables in the molecule structure
	*/
	void initialize()
	{
		int i;
		mass           = (float*) malloc(numatoms*sizeof(float));
		charge         = (float*) malloc(numatoms*sizeof(float));
		chargeij       = (float**) malloc(numatoms*sizeof(float*));
		radius         = (float*) malloc(numatoms*sizeof(float));
		s_hct          = (float*) malloc(numatoms*sizeof(float));
		resnr          = (int*) malloc(numatoms*sizeof(int));
		exclusioncount = (int*) malloc(numatoms*sizeof(int));
		exclusions     = (int**) malloc(numatoms*sizeof(int*));
		name           = (char**) malloc(numatoms*sizeof(char*));
		type           = (char**) malloc(numatoms*sizeof(char*));
		resname        = (char**) malloc(numatoms*sizeof(char*));
		crd            = (float (*)[3])calloc(numatoms,sizeof(*crd));
		tprcrd         = (float**) malloc(numatoms*sizeof(float*));
		velocity       = (float (*)[3])calloc(numatoms,sizeof(*velocity));
		ljtypes        = (int*) malloc(numatoms*sizeof(int));
		nblist         = (int**) malloc(numatoms*sizeof(int*));
		nlist          = (int**) malloc(numatoms*sizeof(int*));
		nbcount        = (int*) malloc(numatoms*sizeof(int));
		ncount         = (int*) malloc(numatoms*sizeof(int));
	
		bondstress     = (float*) malloc(numatoms*sizeof(float));
		anglestress    = (float*) malloc(numatoms*sizeof(float));
		dihedralstress = (float*) malloc(numatoms*sizeof(float));
		vdwstress      = (float*) malloc(numatoms*sizeof(float));
		qqstress       = (float*) malloc(numatoms*sizeof(float));
		solventstress  = (float*) malloc(numatoms*sizeof(float));
		kineticstress  = (float*) malloc(numatoms*sizeof(float));
		//force          = (float**) malloc(numatoms*sizeof(float*));
	
		numbonds        = 0;
		numangles       = 0;
		numpairs        = 0;
		numpdihedrals   = 0;
		numidihedrals   = 0;
		numrbdihedrals  = 0;
		numffterms      = 0;
		numljtypes      = 0;
	
		xbox           = 1.0;
		ybox           = 1.0;
		zbox           = 1.0;
		volume         = 1.0;

		for(i=0;i<numatoms;i++)
		{
			name[i]           = (char*) malloc(5*sizeof(char));
			type[i]           = (char*) malloc(5*sizeof(char));
			resname[i]        = (char*) malloc(5*sizeof(char));
			exclusions[i]     = (int*) malloc(100*sizeof(int));
			exclusioncount[i] = 0;
			//crd[i]            = (float*) malloc(3*sizeof(float));
			tprcrd[i]         = (float*) malloc(3*sizeof(float));
			//velocity[i]       = (float*) malloc(3*sizeof(float));
			nblist[i]         = (int*) malloc(1000*sizeof(int));
			nlist[i]          = (int*) malloc(1000*sizeof(int));
			radius[i]         = -1.0;
			s_hct[i]          = 1.0;
			chargeij[i]       = (float*) malloc(numatoms*sizeof(float));
			
		}
	}

	void setchargeij()
	{
		int i,j;
		for(i=0;i<numatoms;i++)
		{
			for(j=0;j<numatoms;j++)
			{
				chargeij[i][j] = charge[i] * charge[j];
			}
		}
	}	
	/*reports if the radii are missing
	*/
	int missingradii()
	{
		int i;
		for(i=0;i<numatoms;i++)
		{
			if(radius[i] <= 0.0)
			{
				return 1;
			}
		}
		return 0;
	}

	/* set atomic radii to bondi radii
	*/
	void setradii()
	{
		int i;
		for(i=0;i<numatoms;i++)
		{
			switch((int)mass[i])
			{
				case 1:
					radius[i] = 1.1;
					s_hct[i]  = 0.85;
					break;
				case 12:
					radius[i] = 1.7;
					s_hct[i]  = 0.72;
					break;
				case 14:
					radius[i] = 1.55;
					s_hct[i]  = 0.79;
					break;
				case 16:
					radius[i] = 1.52;
					s_hct[i]  = 0.85;
					break;
				case 31:
					radius[i] = 1.80;
					s_hct[i]  = 0.86;
					break;
				case 32:
					radius[i] = 1.80;
					s_hct[i]  = 0.96;
				default:
					radius[i] = 1.50;
					radius[i] = 0.85;
					break;
			}
		}
	}
	
	/* dynamic memory allocation of the bond lists
	 */
	void initializebonds()
	{
		bonds = (struct bond *) malloc(numbonds*sizeof(struct bond));
	}
	/* dynamic memory allocation of the angle lists
	 */
	void initializeangles()
	{
		angles = (struct angle *) malloc(numangles*sizeof(struct angle));
	}
	/* dynamic memory allocation of the improper dihedrals
	 */
	void initializeidihedrals()
	{
		idihedrals = (struct idihedral*) malloc(numidihedrals*sizeof(struct idihedral));
	}

	/* dynamic memory allocation of the proper dihedrals
	 */
	void initializepdihedrals()
	{
		pdihedrals = (struct pdihedral*) malloc(numpdihedrals*sizeof(struct pdihedral));
	}
	/* dynamic memory allocation of the rb dihedrals
	 */
	void initializerbdihedrals()
	{
		rbdihedrals = (struct rbdihedral*) malloc(numrbdihedrals*sizeof(struct rbdihedral));
	}
	/* dynamic memory allocation of the pairs
	 */
	void initializepairs()
	{
		pairs = (struct pair*) malloc(numpairs*sizeof(struct pair));
	}

	/* dynamic memory allocation of the forcefield table
	 */
	void initializeff()
	{
		int i,j;
		forcefield = (float**) malloc(numffterms*sizeof(float*));
		for(i=0;i<numffterms;i++)
		{
			forcefield[i] = (float*) malloc(6*sizeof(float));
			for(j=0;j<6;j++)
			{
				forcefield[i][j] = 0.0;
			}
		}

	}


	/* prints system information. Used for debug purposes
	 */
	void debug()
	{
		int i,j;
		for(i=0;i<numatoms;i++)
		{
			printf("%5s %5s %5d %6s %12.4f %12.4f %d\n",name[i],type[i],resnr[i],resname[i],mass[i],charge[i],exclusioncount[i]);
		}
		/*for(i=0;i<numatoms;i++)
		{
			for(j=0;j<exclusioncount[i];j++)
			{
				printf("%d ",exclusions[i][j]);
			}
			printf("\n");
		}
		for(i=0;i<numffterms;i++)
		{
			for(j=0;j<6;j++)
			{
				printf("%.12f ",forcefield[i][j]);
			}
			printf("\n");
		}
		printf("FudgeQQ: %f\n",fudgeQQ);
		printf("Number of bonds: %d\n",numbonds);
		for(i=0;i<numbonds;i++)
		{
			printf("%d %d %d %f %f \n",bonds[i].atomi,bonds[i].atomj,bonds[i].type,bonds[i].length,bonds[i].fc);
		}
		printf("Number of angles: %d\n",numangles);
		for(i=0;i<numangles;i++)
		{
			printf("%d %d %d %d %f %f\n",angles[i].atomi,angles[i].atomj,angles[i].atomk,angles[i].type,angles[i].angle,angles[i].fc);
		}*/
		/*for(i=0;i<numatoms;i++)
		{
			printf("%f %f %f\n",crd[i][0],crd[i][1],crd[i][2]);
		}
		printf("%d\n",indexgroupcount);
		for(i=0;i<indexgroupcount;i++)
		{
			printf("%d ",indexgroup[i]);
		}*/
		printf("%d\n",numpairs);
		for(i=0;i<numpairs;i++)
		{
			printf("%d %d %d %f %f\n",pairs[i].atomi,pairs[i].atomj,pairs[i].type,pairs[i].c6,pairs[i].c12);
		}

	}
};
