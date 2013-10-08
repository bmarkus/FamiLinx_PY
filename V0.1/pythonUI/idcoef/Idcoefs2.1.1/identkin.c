/* identkin.c */
/*
  Copyright (C) 2004, 2005 Mark Abney

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>
#include "ibdgraph.h"
#include "identkin.h"

void kin2ident(double *ident, Probvec_t *pvec)
{
   int i, j;
   double kincoef[NIDSTATE];
   double kin2idmat[NIDSTATE][NIDSTATE] = {
      {1.00,  0.00, -0.50,  0.00, -0.50,  0.00,  0.50,  0.25,  0.00},
      {0.00,  1.00, -0.50, -1.00, -0.50, -1.00,  0.50,  0.75,  1.00},
      {0.00,  0.00,  2.00,  0.00,  0.00,  0.00, -2.00, -1.00,  0.00},
      {0.00,  0.00,  0.00,  2.00,  0.00,  0.00,  0.00, -1.00, -2.00},
      {0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -2.00, -1.00,  0.00},
      {0.00,  0.00,  0.00,  0.00,  0.00,  2.00,  0.00, -1.00, -2.00},
      {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00,  0.00,  0.00},
      {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00,  0.00},
      {0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  4.00} };

   kincoef[0] = pvec->istate[0];
   kincoef[1] = pvec->istate[1];
   kincoef[2] = pvec->istate[2] + pvec->istate[3];
   kincoef[3] = pvec->istate[4];
   kincoef[4] = pvec->istate[5] + pvec->istate[6];
   kincoef[5] = pvec->istate[7];
   kincoef[6] = pvec->istate[8] + pvec->istate[9];
   kincoef[7] = pvec->istate[10] + pvec->istate[11] + 
      pvec->istate[12] + pvec->istate[13];
   kincoef[8] = pvec->istate[14];

   for (i=0; i<NIDSTATE; i++) {
      ident[i] = 0;
      for (j=0; j<NIDSTATE; j++)
	 ident[i] += kin2idmat[i][j] * kincoef[j];
   }
}

double kindiff(Probvec_t *x1, Probvec_t *x2)
{
   int i;
   double diff, old=0;
   static double kin1[NIDSTATE], kin2[NIDSTATE];

   kin1[0] = x1->istate[0];
   kin1[1] = x1->istate[1];
   kin1[2] = x1->istate[2] + x1->istate[3];
   kin1[3] = x1->istate[4];
   kin1[4] = x1->istate[5] + x1->istate[6];
   kin1[5] = x1->istate[7];
   kin1[6] = x1->istate[8] + x1->istate[9];
   kin1[7] = x1->istate[10] + x1->istate[11] + 
      x1->istate[12] + x1->istate[13];
   kin1[8] = x1->istate[14];

   kin2[0] = x2->istate[0];
   kin2[1] = x2->istate[1];
   kin2[2] = x2->istate[2] + x2->istate[3];
   kin2[3] = x2->istate[4];
   kin2[4] = x2->istate[5] + x2->istate[6];
   kin2[5] = x2->istate[7];
   kin2[6] = x2->istate[8] + x2->istate[9];
   kin2[7] = x2->istate[10] + x2->istate[11] + 
      x2->istate[12] + x2->istate[13];
   kin2[8] = x2->istate[14];

   for (i=0; i<NIDSTATE; i++) {
      diff = fabs(kin1[i]-kin2[i]);
      diff = (old > diff) ? old : diff;
      old = diff;
   }

   for (i=0; i<15; i++){
      diff = fabs(x1->istate[i] - x2->istate[i]);
      diff = (old > diff) ? old : diff;
      old = diff;
   }

   return diff;
}
