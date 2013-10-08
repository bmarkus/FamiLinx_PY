/* ibdgraph.c */

/* Functions for using and creating ibdgraph structures. */
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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ibdgraph.h"
#include "pedigree.h"

Ibdgraph_t ibdgr_init(int count, ...)
/* The variable argument list must be the id numbers of the individuals
   in the list. The number of individuals cannot exceed NGENES.
*/
{
   va_list ap;
   
   Ibdgraph_t gr;
   int i;

   va_start(ap, count);
   if (count > NGENES) {
      fprintf(stderr, "count in ibdgr_init too big.\n");
      fprintf(stderr, "argument values are: ");
      for (i=0; i<count; i++)
	 fprintf(stderr, "%d ", va_arg(ap, int));
      fprintf(stderr, "\n");
      exit(1);
   }

   for (i=0; i<count; i++)
      gr.genelist[i] = va_arg(ap, int);
   for (i=count; i<NGENES; i++)
      gr.genelist[i] = 0;
   for (i=0; i<NGENES; i++)
      gr.connectarr[i] = i;

   return gr;
}

void boundcond4(Ibdgraph_t *node, Probvec_t *prob)
{
   //static int state[NGENES*(NGENES-1)/2];
   //int nconn=0;

   if (areconnected(node->connectarr, 0, 1)) {
      if (areconnected(node->connectarr, 0, 2)) {
	 if (areconnected(node->connectarr, 0, 3)) {
	    prob->istate[0] = 1;
	 }
	 else { //0-1 true, 0-2 true, 0-3 false
	    prob->istate[2] = 1;
	 }
      }
      else { // 0-1 true, 0-2 false
	 if (areconnected(node->connectarr, 0, 3)) {
	    prob->istate[3] = 1;
	 }
	 else { //0-1 true, 0-2 false, 0-3 false
	    if (areconnected(node->connectarr, 2, 3)) {
	       prob->istate[1] = 1;
	    }
	    else { //0-1 true, 0-2 false, 0-3 false, 2-3 false
	       prob->istate[4] = 1;
	    }
	 }
      }
   }
   else { //0-1 false
      if (areconnected(node->connectarr, 0, 2))
	 if (areconnected(node->connectarr, 0, 3)) 
	    prob->istate[5] = 1;
	 else // 0-1 false, 0-2 true, 0-3 false
	    if (areconnected(node->connectarr, 1, 3))
	       prob->istate[8] = 1;
	    else //0-1 false, 0-2 true, 0-3 false, 1-3 false
	       prob->istate[10] = 1;
      else //0-1 false, 0-2 false
	 if (areconnected(node->connectarr, 0, 3))
	    if (areconnected(node->connectarr, 1, 2))
	       prob->istate[9] = 1;
	    else // 0-1 false, 0-2 false, 0-3 true, 1-2 false
	       prob->istate[11] = 1;
	 else //0-1 false, 0-2 false, 0-3 false
	    if (areconnected(node->connectarr, 1, 2))
	       if (areconnected(node->connectarr, 1, 3))
		  prob->istate[6] = 1;
	       else // 0-1 false, 0-2 false, 0-3 false, 1-2 true, 1-3 false
		  prob->istate[12] = 1;
	    else //0-1 false, 0-2 false, 0-3 false, 1-2 false
	       if (areconnected(node->connectarr, 1, 3))
		  prob->istate[13] = 1;
	       else //0-1 false, 0-2 false, 0-3 false, 1-2 false, 1-3 false
		  if (areconnected(node->connectarr, 2, 3))
		     prob->istate[7] = 1;
		  else
		     prob->istate[14] = 1;
   }
}

void newboundcond(Ibdgraph_t *node, Probvec_t *prob)
{
   int state = connectstate(node->connectarr);

   prob->istate[state] = 1;

   return;
}

int connectstate(int *connect)
{
   static int state[64] = { 14,  7, 13, -1, 12, -1, -1,  6, //0-7
			    11, -1, -1, -1,  9, -1, -1, -1, //8-15
			    10, -1,  8, -1, -1, -1, -1, -1, //16-23
			    -1,  5, -1, -1, -1, -1, -1, -1, //24-31
			     4,  1, -1, -1, -1, -1, -1, -1, //32-39
			    -1, -1,  3, -1, -1, -1, -1, -1, //40-47
			    -1, -1, -1, -1,  2, -1, -1, -1, //48-55
			    -1, -1, -1, -1, -1, -1, -1,  0 }; //56-63

   int ab = connect[0]==connect[1]; //areconnected(connect, 0, 1);
   int ac = connect[0]==connect[2]; //areconnected(connect, 0, 2);
   int ad = connect[0]==connect[3]; //areconnected(connect, 0, 3);
   int bc = connect[1]==connect[2]; //areconnected(connect, 1, 2);
   int bd = connect[1]==connect[3]; //areconnected(connect, 1, 3);
   int cd = connect[2]==connect[3]; //areconnected(connect, 2, 3);

   int idx = (ab<<5) + (ac<<4) + (ad<<3) + (bc<<2) + (bd<<1) + cd;

   /*if (idx < 0) {
      printf("ERROR: invalid index in newboundcond: %d\n", idx);
      //printibdgr(node);
      exit(10);
   }
   */
   return state[idx];
}

void printibdgr(Ibdgraph_t *graph)
{
   int i;

   printf("id:   ");
   for (i=0; i<NGENES; i++) {
      printf("%3d ", graph->genelist[i]);
   }
   printf("\nconn: ");
   for (i=0; i<NGENES; i++)
      printf("%3d ", graph->connectarr[i]);
   printf("\n");
   
   return;
}

/*
int areconnected(int *connectar, int idx1, int idx2)
{
   if (connectar[idx1] == connectar[idx2])
      return 1;
   else
      return 0;
   
   return (connectar[idx1] == connectar[idx2]);
}
*/

void connected2parent(Ibdgraph_t *node, int prnt, int idx)
{
   int i, id;

   id = node->genelist[idx];
   for (i=0; i<NGENES; i++) {
      if (areconnected(node->connectarr, idx, i) && node->genelist[i] == id)
	 //node->genelist[i] = (node->genelist[i]) ->parent[prnt];
	 node->genelist[i] = getparent(prnt, node->genelist[i]);
   }
   return;
}

void founder2parent(Ibdgraph_t *node, int prnt, int idx)
{
   int i, id;

   id = node->genelist[idx];
   for (i=0; i<NGENES; i++) {
      if (areconnected(node->connectarr, idx, i) && node->genelist[i] == id)
	 if (prnt)
	    node->genelist[i] = - node->genelist[i];
   }
   return;
}

void connect(Ibdgraph_t *node, int idx1, int idx2)
{
   int *conn_state = node->connectarr;
   int i, state;

   if (conn_state[idx1] == conn_state[idx2])
      return;
   for (state=conn_state[idx2], i=0; i<NGENES; i++) {
      if (conn_state[i] == state)
	 conn_state[i] = conn_state[idx1];
   }
   return;
}

void initpvec(Probvec_t *pvec)
{
   int i;
   for (i=0; i<NCOEF; i++)
      pvec -> istate[i] = 0;
   return;
}

void add2pvec(Probvec_t *orig, Probvec_t *x)
{
   int i;

   /*for (i=0; i<NCOEF; i++)
      (orig->istate[i])/(1<<(orig->coef)) += (x->istate[i])/(1<<(x->coef));
   */

   for (i=0; i<NCOEF; i++)
      orig->istate[i] += x->istate[i];
   return;
}

void mult2pvec(double coef, Probvec_t *x)
{
   int i;
   double *pr=x->istate;

   for (i=0; i<NCOEF; i++)
      //x->istate[i] *= coef;
      pr[i] *= coef;
   return;
}

void printpvec(Probvec_t *x)
{
   int i;

   puts("State\tProbability\n");
   for (i=0; i<NCOEF; i++)
      printf("%d\t%.12g\n", i, x->istate[i]);
   return;
}
      
double pvecdiff(Probvec_t *x1, Probvec_t *x2)
{
   int i;
   double diff = 0;
   double a, b, old = 0;

   for (i=0; i<NCOEF; i++) {
      a = x1->istate[i];
      b = x2->istate[i];
      diff = fabs(a-b);
      diff = (old > diff) ? old : diff;
      old = diff;
   }
   return diff;
}
