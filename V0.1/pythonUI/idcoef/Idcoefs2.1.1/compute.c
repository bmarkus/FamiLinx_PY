/* file: compute.c - Version 1*/

/* All the necessary routines to do the desired computation. In this case,
   we want to compute the identity coefficients for pairs of people.
*/
/*
  Copyright (C) 2004, 2005, 2013 Mark Abney

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

#include <stdio.h>
#include <stdlib.h>
#include "compute.h"
#include "ibdgraph.h"
#include "pedigree.h"
#include "followbranch.h"
#include "identkin.h"

static int **imatrix(long nrl, long nrh, long ncl, long nch);

static int **Pair; // All pairs. Pair[i][0 or 1]= pair i person 0 or 1
static int Npairs; // Number of pairs.

void readsample(char *sampfile)
/* Construct the sample. That is, the sample is the set of pairs for which 
   identity coefficients should be computed. We read from a file a list of
   external id's (findivs). If, within the file, there is one findiv per line,
   construct the sample by taking all possible pairs of findivs. If, on each
   line, there are two findivs, the sample consists of all the pairs in the file.

   Arguments and variables:
   sampfile - pointer to name of the file that holds the findivs from which
              to construct the sample.
   Pair - matrix which holds the sample. Pair[i][j] is the jth individual of
          pair i and returns an internal id. 0 < i <= Npairs, j = 0 or 1.
   Npairs - number of pairs for which identity coefficients will be computed.
*/
{
   FILE *sampfl;
   char line[MAXLEN];
   int i, j, k, tmp1, tmp2;
   int ninput;
   int nlines;
   int *samp; // Holds the individual ids before making pairs.
   int found = 1; // Assume we'll find everyone in the sample file.

   /* 
    * We want to read in the study sample and make sure they all exist
    * in the pedigree somewhere.
    */
    if ( (sampfl=fopen(sampfile,"r")) == NULL) {
      fprintf(stderr, "Cannot open sample file.\n");
      exit(1);
    }

    /* Count the number of pairs for whom we will compute id coefs */
   for(nlines=0; fgets(line, MAXLEN, sampfl);nlines++) {}
   rewind(sampfl);
   samp = (int *) calloc(nlines, sizeof (int));

   /*
    * If the first line has only one id, then we will construct all possible
    * pairs of ids. If  lines have two ids, then use only these pairs. After
    * getting the first line and checking the number of id's verify that these
    * id's are in the pedigree and place the corresponding internal id into
    * either the samp vector or Pair matrix.
    */
   fgets(line, MAXLEN, sampfl);
   ninput = sscanf(line, "%d %d", &tmp1, &tmp2);
   if (ninput == 1) 
      Npairs = nlines * (nlines+1)/2;
   else
      Npairs = nlines;
   Pair = imatrix(0, Npairs, 0, 1);
   /* Make sure people are in the pedigree */
   if (ninput == 1) {
      samp[0] = findid(tmp1);
      if (!samp[0])
	 found = 0;
   }
   else {
      Pair[0][0] = findid(tmp1);
      if (!Pair[0][0]) found = 0;
      Pair[0][1] = findid(tmp2);
      if (!Pair[0][1]) found = 0;
   }
   /* Now go through the rest of the lines in the file. */
   for (i=1; i<nlines; i++) {
      fgets(line, MAXLEN, sampfl);
      sscanf(line, "%d %d", &tmp1, &tmp2);
      /* Make sure people are in the pedigree */
      if (ninput == 1) {
	 samp[i] = findid(tmp1);
	 if (!samp[i])
	    found = 0;
      }
      else {
	 Pair[i][0] = findid(tmp1);
	 if (!Pair[i][0]) found = 0;
	 Pair[i][1] = findid(tmp2);
	 if (!Pair[i][1]) found = 0;
      }
   }
   
   if (!found) {
       fprintf(stderr, "Some sample individuals not found in pedigree.\n");
       exit (2);
   }

   /*
    * If the input file only had one id per line, construct all possible pairs
    * of individuals and place each pair in the Pair matrix.
    */
   if (ninput == 1) {
      for (k=0,i=0; i<nlines; i++) {
	 for (j=i; j<nlines; j++, k++) {
	    Pair[k][0] = samp[i];
	    Pair[k][1] = samp[j];
	 }
      }
   }

   /*for (i=0; i<Npairs; i++)
      printf("%d  %d\n", Pair[i][0], Pair[i][1]);
      exit(10);*/

   free(samp);
}

void computeidcoefs(const char *outfile)
{
   Ibdgraph_t node;
   Probvec_t prob;
   double idcoef[NIDSTATE];
   int i, j;
   FILE *outfl;

   outfl = fopen(outfile, "w");
   if (outfl == NULL) {
      fprintf(stderr,"Could not open %s for output.\n", outfile);
      exit(10);
   }
   for (i=0; i<Npairs; i++) {
      minimalped(2, Pair[i]);
      //printminped();
      node = ibdgr_init(4, Pair[i][0], Pair[i][0], Pair[i][1], Pair[i][1]);
      //prob = nodeprob_sym(&node);
      prob = nodeprob(&node);
      //printpvec(&prob);
      kin2ident(idcoef, &prob);
      fprintf(outfl,"%d\t%d", getfindiv(Pair[i][0]), getfindiv(Pair[i][1]));
      for (j=0; j<NIDSTATE; j++)
	 fprintf(outfl,"\t%.9g", idcoef[j]);
      fprintf(outfl, "\n");
   }

}


#define NR_END 1
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        int **m;

        /* allocate pointers to rows */
        m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
        if (!m) { 
	   fprintf(stderr, "allocation failure. Out of memory.\n");
	   exit(1);
	}
        m += NR_END;
        m -= nrl;


        /* allocate rows and set pointers to them */
        m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
        if (!m[nrl]) {
	   fprintf(stderr, "allocation failure. Out of memory.\n");
	   exit(1);
	}
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

