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

/* Functions that perform energy calculations
 * This is primarily for debug purposes
 */

float bondenergy(struct topology *mol);
float angleenergy(struct topology *mol);
float pairenergy(struct topology *mol);
float pdihedralenergy(struct topology *mol);
float idihedralenergy(struct topology *mol);
float nbenergy(struct topology *mol);
float rbdihedralenergy(struct topology *mol);
float kineticenergy(struct topology *mol);
float gbenergy(struct topology *mol);
