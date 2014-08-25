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

void help();
void commandlineargs(int argc, char *argv[]);
void readtpr(char tprfilename[],struct topology *mol);
int  wordcount(char *line);
void readindexgroup(char *ndxfilename,struct topology *mol);
void readpdb(char pdbfilename[],struct topology *mol);
void getnblist(struct topology *mol);
void getnblist_grid(struct topology *mol);
void getnlist_grid(struct topology *mol);
void writeout(char outfilename[],struct topology *mol);
void writestat(char outfilename[],struct topology *mol,int numframes);
void addneighbor(struct topology *mol, int i, int j);
int readframe(XDRFILE *infile, struct topology *mol);
