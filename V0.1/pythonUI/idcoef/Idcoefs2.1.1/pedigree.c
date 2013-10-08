/* pedigree.c - Version 5 */

/* Code for dealing with pedigree information.

   Version 3 differs from version 2 in that the parents within a Person_t
   structure are recorded by their internal id rather than as pointers to
   Person_t structures.

   Version 5 (no version 4) restructures some code to remove the study sample
   code. This information is now included in another file which also includes
   the code to compute the desired quantities.
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
#include "pedigree.h"

typedef struct person_tag {
   int id; /* Sequential id used internally by the programs. */
   int findiv; /* The id defined externally by the user for this person. */
   int parent[2]; // parent[0] id number of mother, parent[1] id of father.
} Person_t;

typedef struct pedigree {
   int npeople; /* Total number in the pedigree. */
   Person_t *member;
} Pedigree_t;

struct idmap {
   int id; // internal id.
   int fid; // external id (id used in input file).
};

static void addancestors(int *last_p, int *pedlist, int *inlist);
static int **imatrix(long nrl, long nrh, long ncl, long nch);
static void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
static int fidcomp(const void *e1, const void *e2);
static void printped(void);

static struct idmap *Fid2id; // Used to map external id's to internal id's.
Pedigree_t Pedx;
static int *Istruefounder;
static int *Isfounder;
static int *Inminped; // = (int *)calloc(npeop+1, sizeof(int));
static int *Isinstudy; // = (int *)calloc(npeop+1, sizeof(int));
static int *Child; // = (int *)calloc(npeop+1, sizeof(int));
static int *Founderlist; // = (int *)calloc(npeop+1, sizeof(int));
static int *Minped; // = (int *)calloc(npeop+1, sizeof(int));
static int *Noffspring; // = (int *)calloc(npeop+1, sizeof(int));

int getparent(int whichpar, int indiv)
{
   //printf("Check person.\n"); fflush(stdout);
   /*if (indiv > Pedx.npeople) {
      fprintf(stderr,"Error: person %d does not exist in pedigree.\n", indiv);
      exit(1);
   }
   */
   //return Pedx.member[indiv].parent[whichpar]->id;
   return Pedx.member[indiv].parent[whichpar];

}

int getfindiv(int id)
{
   /*if (id > Pedx.npeople) {
      fprintf(stderr,"Error: person %d does not exist in pedigree.\n", id);
      exit(1);
   }
   */
   return Pedx.member[id].findiv;
}
int getnpeop(void)
{
   return Pedx.npeople;
}

void printped(void)
{
   int i, parent0, parent1;
   for (i=1; i<=Pedx.npeople; i++) {
      parent0 = Pedx.member[i].parent[0];
      parent1 = Pedx.member[i].parent[1];
      fprintf(stdout, "%d (%d) \t%d (%d)\t%d (%d)\n", Pedx.member[i].id,
	      Pedx.member[i].findiv, parent0, Pedx.member[parent0].findiv,
	      parent1, Pedx.member[parent1].findiv);
   }
}

int fidcomp(const void *e1, const void *e2)
{
   int v1 = ((struct idmap *)e1)->fid;
   int v2 = ((struct idmap *)e2)->fid;
   return (v1<v2) ? -1 : (v1>v2) ? 1 : 0;
}

void readped(char *pedfile)
/* Read the pedigree file and construct the Pedigree structure. The Pedigree
   structure has two components, the first an integer with the number of 
   people in the pedigree, the second a pointer to an array of Person structures.
   There is one Person structure for each individual in the pedigree, ordered
   such that parents come before offspring. A Person structure contains the
   internal and external id's of the individual and a two element array containing
   the internal id's of the person's parents (element 0 is for the mother and
   element 1 is for the father).

   Arguments/variable:
   pedfile - string with the name of the file containing the pedigree.
   Pedx - stucture holding the pedigree data.
   Fid2id - array holding the mapping of external to internal id's.
*/
{
    Person_t *pedmem; // points to array holding all the people in the pedigree
    FILE *pedfl;
    int i, j, npeop;
    char line[MAXLEN];
    int *fid, **fparent; // temporary storage for external id's of people and parents.
    int found;
    int leaveprog = 0;

    if ( (pedfl=fopen(pedfile,"r")) == NULL) {
	fprintf(stderr, "Cannot open pedigree file.\n");
	exit(1);
    }

    /* Count the number of people in the pedigree. */
    for(npeop=0; fgets(line, MAXLEN, pedfl);npeop++) {}
    rewind(pedfl);
    Pedx.npeople = npeop;
    fprintf(stdout, "npeop = %d\n", Pedx.npeople); fflush(stdout);

    fid = calloc(npeop+1, sizeof(int));
    fparent = imatrix(0, npeop, 0, 1);
    pedmem = calloc(npeop+1, sizeof(Person_t));
    Pedx.member = pedmem;
    
    Fid2id = (struct idmap *)calloc(npeop+1, sizeof (struct idmap));
    
    /*
     * The following arrays are needed later when finding the minimal pedigrees
     * for a set of individuals. Doing computations on the minimal pedigree is
     * more efficient than using the whole pedigree.
     */
    Istruefounder = (int *)calloc(npeop+1, sizeof(int)); // founder as defined in pedigree file
    Isfounder = (int *) calloc(npeop+1, sizeof(int)); // founder in the derived minimal pedigree

    Inminped = (int *)calloc(npeop+1, sizeof(int));
    Isinstudy = (int *)calloc(npeop+1, sizeof(int));
    Child = (int *)calloc(npeop+1, sizeof(int));
    Founderlist = (int *)calloc(npeop+1, sizeof(int));
    Minped = (int *)calloc(npeop+1, sizeof(int));
    Noffspring = (int *)calloc(npeop+1, sizeof(int));

    /*
     * For each person in the pedigree file assign an internal id equal to the line
     * number of the file on which the person appears. Store the external and internal
     * id's of the person in proper pedigree structure element and in Fid2id to allow
     * mapping external id to internal id later. Also store the external id's of the 
     * person's parents in fparent. Later we will make sure that parents come before
     * children and fill in the parent component for each person in the pedigree.
     */

    /* Person 0 doesn't actually exist. */
    fid[0] = 0;
    fparent[0][0] = 0;
    fparent[0][1] = 0;
    for (i=1; i<=npeop; i++) {
	fgets(line, MAXLEN, pedfl);
	sscanf(line, "%d %d %d", &fid[i], &fparent[i][0], &fparent[i][1]); 
	if (fparent[i][0] == 0 && fparent[i][1] != 0) {
	    fprintf(stderr, "%d has only one parent specified. All non-founders ",fid[i]);
	    fprintf(stderr, "must have both parents specified.");
	    leaveprog = 1;
	}
	if (fparent[i][1] == 0 && fparent[i][0] != 0) {
	    fprintf(stderr, "%d has only one parent specified. All non-founders ",fid[i]);
	    fprintf(stderr, "must have both parents specified.");
	    leaveprog = 1;
	}

	pedmem[i].id = i;
	pedmem[i].findiv = fid[i];
	Fid2id[i].id = i;
	Fid2id[i].fid = fid[i];
	if (fparent[i][0] == 0 && fparent[i][1] == 0)
	    Istruefounder[i] = 1;
	else
	    Istruefounder[i] = 0;
    }
    /* Sort to make looking up a findiv quicker. */
    qsort(Fid2id, npeop+1, sizeof(struct idmap), fidcomp);

    /*
     * For the parents of individual 'i', look through the pedigree that comes
     * before 'i' and make sure the parents are there. Once the parent is found
     * at element 'j' fill in the parent array in the pedigree with 'j', since
     * this will be the internal id of the parent.
     */
    /* Make the parent structure member have the right values. Do mother first.*/
    for (i=1; i<=npeop; i++) {
	for (found=0, j=0; j<i; j++) { //Require parents to be before children
	    if (fid[j] == fparent[i][0]) {
		found = 1;
		break;
	    }
	}
	if (!found) {
	    fprintf(stderr, "Mother of %d not found in pedigree.\n", fid[i]);
	    leaveprog = 1;
	}
	pedmem[i].parent[0] = j;

	/* Now check the father. */
	for (found=0, j=0; j<i; j++) { // Require parents to be before children
	    if (fid[j] == fparent[i][1]) {
		found = 1;
		break;
	    }
	}
	if (!found) {
	    fprintf(stderr, "Father of %d not found in pedigree.\n", fid[i]);
	    leaveprog = 1;
	}
	pedmem[i].parent[1] = j;
    }

#ifdef DEBUG_PED
    puts("Debug mesg:");
    for (i=1; i<=npeop; i++)
	fprintf(stdout, "%d (%d) \t%d (%d)\t%d (%d)\n", pedmem[i].id,
		fid[i], 
		pedmem[i].parent[0], fparent[i][0], pedmem[i].parent[1], 
		fparent[i][1]);
    puts("Debug mesg:");
    for (i=1; i<=npeop; i++)
	fprintf(stdout, "%d (%d) \t%d (%d)\t%d (%d)\n", Pedx.member[i].id,
		fid[i], Pedx.member[i].parent[0], fparent[i][0], 
		Pedx.member[i].parent[1], fparent[i][1]);
#endif


    if (leaveprog) {
	puts("Fix the files and try again. Bye...\n");
	exit(1);
    }

    fclose(pedfl);
    free(fid);
    free_imatrix(fparent, 0, npeop, 0, 1);
}

int findid(int fid)
{
   struct idmap target, *result;

   target.fid = fid;
   result = bsearch(&target, Fid2id, Pedx.npeople+1, sizeof(struct idmap), fidcomp);
   if (result)
      return result->id;
   else 
      fprintf(stderr, "%d not found in pedigree.\n", fid);
   return 0;
}


int isafounder(int id)
{
   /*if (id > Pedx.npeople ) {
      fprintf(stderr, "Illegal id (%d) in isafounder.\n", id);
      exit(2);
   }
   */
   if (Istruefounder[id] || Isfounder[id])
      return 1;
   else
      return 0;
}

void minimalped(int nsample, int *samplelist)
/* Find the minimal pedigree connecting the people in the sample list.
   The minimal pedigree consists of only those people in the list and
   their ancestors. Superfluous founders are removed by assigning the 
   child of the superfluous founder pair as a founder.
*/
{
   int i, lastidx, nminped, fid, fidx, otherparent;

   Inminped[0] = 1;
   for (i=1; i<=Pedx.npeople; i++) {
      Inminped[i] = 0;
      Isinstudy[i] = 0;
      Child[i] = 0;
      Founderlist[i] = 0;
      Minped[i] = 0;
      Noffspring[i] = 0;
      Isfounder[i] = 0;
   }

   /*
    * Make sure everyone in the sample is in the minimal pedigree
    * and add all of the ancestors for each person.
    */
   for (lastidx=0, i=0; i<nsample; i++) {
      Isinstudy[samplelist[i]] = 1;
      if (!Inminped[samplelist[i]]) {
	 Inminped[samplelist[i]] = 1;
	 Minped[++lastidx] = samplelist[i];
	 addancestors(&lastidx, Minped, Inminped);
      }
   }

   /*
    * We want to eliminate founder couples that have only one child
    * and make the child the founder. If the new founder (ie the child)
    * is now part of a founder pair with only one child, repeat the process.
    * To do this, first count the number of offspring everybody has, then
    * go through each founder and check to see if he/she should be removed.
    */
   nminped = lastidx;
   for (fidx=0, i=1; i<=nminped; i++) {
      Noffspring[getparent(0, Minped[i])]++;
      Noffspring[getparent(1, Minped[i])]++;
      Child[getparent(0, Minped[i])] = Minped[i];
      Child[getparent(1, Minped[i])] = Minped[i];
      if (isafounder(Minped[i])) 
	 Founderlist[fidx++] = Minped[i];
   }
   //printminped();
   /*
    * if founder is not in the study sample and has one offspring and
    * the offspring's other parent is a founder not in the study sample
    * with only the offspring as a child, make the offspring a founder
    * and add the offspring to the end of the founder list.
    */
   for (i=0; (fid=Founderlist[i]) != 0; i++) {
       if ( !Inminped[fid]) continue;
       if ( !Isinstudy[fid] && Noffspring[fid]==1) {
           otherparent = getparent(0, Child[fid]);
           if (otherparent == fid) 
               otherparent = getparent(1, Child[fid]);
           if (isafounder(otherparent) && !Isinstudy[otherparent] && 
               Noffspring[otherparent]==1) {
               Isfounder[Child[fid]] = 1;
               Founderlist[fidx++] = Child[fid];
               Inminped[fid] = 0;
               Inminped[otherparent] = 0;
           }
       }
   }
   //printminped();
   return;
}

void printminped(void)
{
   FILE *pedfl = fopen("minimal.pedigree", "w");
   int i;
   for (i=1; Minped[i] != 0 && i<=Pedx.npeople; i++)
      if (Inminped[Minped[i]]) {
	 if (isafounder(Minped[i]))
             fprintf(pedfl, "%d\t0\t0\n", Pedx.member[Minped[i]].findiv);
         else
             fprintf(pedfl,"%d\t%d\t%d\n", Pedx.member[Minped[i]].findiv,
                     Pedx.member[getparent(0,Minped[i])].findiv,
                     Pedx.member[getparent(1, Minped[i])].findiv);
      }
   fclose(pedfl);
}

void addancestors(int *last_p, int *pedlist, int *inlist)
{
   int idx, mother, father;

   for (idx=*last_p; idx <= *last_p; idx++) {
      mother = getparent(0, pedlist[idx]);
      if (!inlist[mother]) {
	 pedlist[++(*last_p)] = mother;
	 inlist[mother] = 1;
      }
      father = getparent(1, pedlist[idx]);
      if (!inlist[father]) {
	 pedlist[++(*last_p)] = father;
	 inlist[father] = 1;
      }
   }
   return;
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

#define FREE_ARG char*
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}
