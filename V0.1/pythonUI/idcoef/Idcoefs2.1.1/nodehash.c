/* nodehash.c */

/*
 * A hash implementation to store and recall kinship probabilities
 * already computed in the nodeprob routine. This version is for
 * computations that are not conditional on genotypes.
 *
 */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "ibdgraph.h"
#include "nodehash.h"

#if ULONG_MAX == 18446744073709551615ULL
#define LONG_IS_64BITS
#endif

//void printclasses(int ncl, Equiv_class_ng *class_p);

#if defined(LONG_IS_64BITS)
typedef  unsigned long       ub8;   /* unsigned 8-byte quantities */
#else
typedef  unsigned long long  ub8;   /* unsigned 8-byte quantities */
#endif

typedef  unsigned int  ub4;   /* unsigned 4-byte quantities */
typedef  unsigned char ub1;


#define NID 4
struct hash_cell {
   Probvec_t prob;
   //int connectarr[NGENES];
   int nodeid[NID];
   int cstate;
   int isempty;
   struct hash_cell *next;
   struct hash_cell *prev; /* In the hash table, this points to the last
			      cell in the collision list. */
};

static unsigned long Ncells; // make long for 64 bits
static unsigned long Tablesize; // make long for 64 bits
static unsigned Bsize;
static struct hash_cell *Table;
static unsigned long Maxcells; // make long for 64 bits
static int Npeop;
static unsigned long Mask; // make long for 64 bits

static int Nmade;
static int Recycle=0;
/*#define ANALSZ 275000
  static int Hashanal[ANALSZ];*/
static unsigned Call[15];

#if defined(LONG_IS_64BITS)
static unsigned long hashfunc(Ibdgraph_t *node); // make long for 64 bits
#else
static unsigned hashfunc(Ibdgraph_t *node); // make long for 64 bits
#endif

static unsigned fnv_32(void *buf, size_t len, unsigned hval);
static unsigned bbhash(void *key, size_t length, unsigned initval);

static unsigned bbhash2(
			void *key,        /* the key */
			register ub4  length,   /* the length of the key, in ub4s */
			register ub4  initval  /* the previous hash, or an arbitrary value */
			);

static ub8 bbhash64( 
		    register ub8 *k,        /* the key */
		    register ub8  length,   /* the length of the key */
		    register ub8  level    /* the previous hash, or an arbitrary value */
		    );

#define COLLISION_DEPTH 1


void hashinit(int memsize, int npeop)
   /* Initialize the hash table. 
      memsize is the maximum amount of memory for hash cells that will 
      be allocated during execution of the program. memsize is given
      in megabytes.
   */
{
   long i;

   for (i=0; i<15; i++)
      Call[i] = 0;
   
   Npeop = npeop;
   Nmade = 0;
   /*for (i=0; i<ANALSZ; i++)
     Hashanal[i] = 0;*/

   Maxcells = ((unsigned long)memsize)*1000000L/sizeof(struct hash_cell);

   for (Tablesize=1, i=1; Tablesize<Maxcells/COLLISION_DEPTH; i++, Tablesize<<=1)
      ;
   Bsize = i - 2;
   Tablesize >>= 1;
   /* 
    * We want to make sure that the Table size is small enough that there is
    * enough memory available to make a collision list when collisions occur.
    */
   if (Maxcells - Tablesize < Tablesize>>1) {
      Tablesize >>= 1;
      Bsize--;
   }
   Ncells = Tablesize;
   Mask = (1UL << Bsize) - 1; // make long for 64 bits

#if defined(DEBUG) || defined(DEBUG_HASH)
   printf("sizeof hash_cell = %lu\n", sizeof(struct hash_cell));
   printf("Maxcells = %lu\n", Maxcells);
   printf("Tablesize = %lu\n", Tablesize);
   printf("Bsize = %u\n", Bsize);
   printf("Mask = %x\n", Mask);
#endif

   Table = calloc(Tablesize, sizeof(struct hash_cell));
   if (Table == NULL) {
      fprintf(stderr, "Not enough memory for hash table.\n");
      exit(2);
   }

   for (i=0; i<Tablesize; i++) {
      Table[i].isempty = 1;
      Table[i].next = NULL;
      Table[i].prev = NULL;
   }
}

#define FNV_INIT ((unsigned)0x811c9dc5)
#if defined(LONG_IS_64BITS)
unsigned long hashfunc(Ibdgraph_t *node) // make long for 64 bits
#else
unsigned hashfunc(Ibdgraph_t *node) // make into long for 64 bits
#endif
{
   int i, cstate;
   unsigned hash32; // make into long for 64 bits
   unsigned long hash64;
#if defined(LONG_IS_64BITS)
   static unsigned long geneid[NID]; // static to avoid reallocation with each call.
   static unsigned long cstate_p[1];
   cstate_p[0] = (unsigned long) (node->connstate);
#else
   static unsigned int geneid[NID];
   static unsigned int cstate_p[1];
   cstate_p[0] = (unsigned int) (node->connstate);
#endif

   /*int id0 = node->genelist[0];
     int id1 = node->genelist[1];
     int id2 = node->genelist[2];
     int id3 = node->genelist[3];
   */
   /*
    * The two component geneid is for the symmetric version
    * where we don't care about the order of the first two ids
    * of the last two ids. If we do care, then just assign each 
    * geneid to the values in node->genelist.
    */
   /*geneid[0] = (Npeop+1) * (id0 + id1) + abs(id0 - id1);
   geneid[1] = (Npeop+1) * (id2 + id3) + abs(id2 - id3);
   */
   for (i=0; i<NID; i++)
      geneid[i] = (unsigned) (node->genelist[i]);

   //cstate = connectstate(node->connectarr);
   //hash32 = fnv_32(node->connectarr, sizeof(node->connectarr), hash32);

   //hash32 = fnv_32(geneid, sizeof(geneid), FNV_INIT);
   //hash32 = fnv_32(&node->connstate, sizeof(node->connstate), hash32);

#if defined(LONG_IS_64BITS)
   hash64 = bbhash64(geneid, NID, (unsigned long)FNV_INIT);
   hash64 = bbhash64(cstate_p, 1, hash64);
   hash64 &= Mask;

   return hash64;
#else

   hash32 = bbhash2(geneid, NID, FNV_INIT);
   hash32 = bbhash2(&node->connstate, 1, hash32);

   //printf("hash32 = %u\n", hash32);
   /*
    * hash32 now has a 32 bit hash value for the node. Now we need to 
    * fold that value down to the size of the hash table.
    */
   //hash32 = (hash32 >> Bsize) ^ (hash32 & Mask); // For FNV hash
   hash32 &= Mask;

   //printf("folded hash32 = %u\n", hash32);
   return hash32;
#endif
}

void hashstore(Ibdgraph_t *node, Probvec_t *prob)
{
   int i;
   struct hash_cell *newcell_p;
#if defined(LONG_IS_64BITS)
   unsigned long hv;
#else
   unsigned hv; // make long for 64 bits
#endif
   static int geneid[NID];
   /*int id0 = node->genelist[0];
     int id1 = node->genelist[1];
     int id2 = node->genelist[2];
     int id3 = node->genelist[3];
   */
   int flag=0;

   /*
    * The two component geneid is for the symmetric version
    * where we don't care about the order of the first two ids
    * of the last two ids. If we do care, then just assign each 
    * geneid to the values in node->genelist.
    */
   /*geneid[0] = (Npeop+1) * (id0 + id1) + abs(id0 - id1);
   geneid[1] = (Npeop+1) * (id2 + id3) + abs(id2 - id3);
   */
   for (i=0; i<NID; i++)
      geneid[i] = node->genelist[i];
   
   hv = hashfunc(node);
   if (hv > Tablesize){
       fprintf(stderr,"hv = %lu, too big in hashstore. Tablesize only %lu.\n", hv, Tablesize);
       exit(1);
   }
   /* Goto entry in hash table and place there, if empty. */
   if (Table[hv].isempty) {
      for (i=0; i<NCOEF; i++)
	 Table[hv].prob.istate[i] = prob->istate[i];

      Table[hv].cstate = node->connstate;
      for (i=0; i<NID; i++)
	 Table[hv].nodeid[i] = geneid[i];
      Table[hv].isempty = 0;

      return;
   }
   /*
    * If we haven't hit the max number of cells, create a new one and add
    * it to the beginning of the collision list. 
    */
   if (Ncells < Maxcells) {
      newcell_p = malloc(sizeof(struct hash_cell));
      if (!newcell_p) {
	 fprintf(stderr, "Out of memory.\n");
	 exit(1);
      }
      Ncells++;
      Nmade++;

      newcell_p->next = Table[hv].next;
      for (i=0; i<NCOEF; i++)
	 newcell_p->prob.istate[i] = prob->istate[i];
      /*for (i=0; i<NGENES; i++)
	newcell_p->connectarr[i] = node->connectarr[i];*/
      newcell_p->cstate = node->connstate;
      for (i=0; i<NID; i++)
	 newcell_p->nodeid[i] = geneid[i];
      newcell_p->prev = &Table[hv];
      Table[hv].next = newcell_p;
      if (newcell_p->next == NULL) /* newcell is the end of collision list */
	 Table[hv].prev = newcell_p;
      else
	 newcell_p->next->prev = newcell_p;
      //if (flag) printf("Added to collision list.\n");
   }
   else { /* Are at the max number of cells */
      Recycle = 1;
      /* Place the new cell info at end of list and move it to the front. */
      if (Table[hv].prev != NULL && Table[hv].prev != Table[hv].next) {
	 /* The case where there are at least two cells in the list. */
	 for (i=0; i<NCOEF; i++)
	    Table[hv].prev->prob.istate[i] = prob->istate[i];
	 /*for (i=0; i<NGENES; i++)
	   Table[hv].prev->connectarr[i] = node->connectarr[i];*/
	 Table[hv].prev->cstate = node->connstate;
	 for (i=0; i<NID; i++)
	    Table[hv].prev->nodeid[i] = geneid[i];
	 Table[hv].prev->next = Table[hv].next;
	 Table[hv].prev->prev->next = NULL;
	 Table[hv].next = Table[hv].prev;
	 Table[hv].prev = Table[hv].prev->prev;
	 Table[hv].next->prev = &Table[hv];
	 Table[hv].next->next->prev = Table[hv].next;
	 //if (flag) printf("Replaced cell in long list.\n");
      }
      else if (Table[hv].prev != NULL) {
	 /* Now case where there is only one cell in the collision list. */
	 for (i=0; i<NCOEF; i++)
	    Table[hv].next->prob.istate[i] = prob->istate[i];
	 /*for (i=0; i<NGENES; i++)
	   Table[hv].next->connectarr[i] = node->connectarr[i];*/
	 Table[hv].next->cstate = node->connstate;
	 for (i=0; i<NID; i++)
	    Table[hv].next->nodeid[i] = geneid[i];
	 //if (flag) printf("Replaced cell in short list.\n");
      }
      else {
	 /* No collision list, put everything in the table. */
	 for (i=0; i<NCOEF; i++)
	    Table[hv].prob.istate[i] = prob->istate[i];
	 /*for (i=0; i<NGENES; i++)
	   Table[hv].connectarr[i] = node->connectarr[i];*/
	 Table[hv].cstate = node->connstate;
	 for (i=0; i<NID; i++)
	    Table[hv].nodeid[i] = geneid[i];
	 //if (flag) printf("Replaced cell in table.\n");
      }
   }
   return;
}

int hashfind(Ibdgraph_t *node, Probvec_t *prob)
   /* Search through hash to find the cell with equivalence classes
      specified in class_p. Returns the probability for the equivalence
      classes. If the search was unsuccessful, returns 0.
   */
{
   int i, connstate=node->connstate;
   int cst, cst2=0;
   int *nid, *nid2=NULL;
   struct hash_cell *cell_p, *cell_p2;
#if defined(LONG_IS_64BITS)
   unsigned long hv;
#else
   unsigned hv; // make long for 64 bits
#endif
   static int geneid[NID];
   struct hash_cell *table;
   table = Table;

   /*int id0 = node->genelist[0];
     int id1 = node->genelist[1];
     int id2 = node->genelist[2];
     int id3 = node->genelist[3];

     int flag=0;
   */
   /*
    * The two component geneid is for the symmetric version
    * where we don't care about the order of the first two ids
    * of the last two ids. If we do care, then just assign each 
    * geneid to the values in node->genelist.
    */
   /*geneid[0] = (Npeop+1) * (id0 + id1) + abs(id0 - id1);
   geneid[1] = (Npeop+1) * (id2 + id3) + abs(id2 - id3);
   */
   for (i=0; i<NID; i++)
      geneid[i] = node->genelist[i];
   
   //connstate = connectstate(node->connectarr);
   Call[connstate]++;

   hv = hashfunc(node);
   if (hv > Tablesize) {
       printf("hv = %lu too big in hashfind. Tablesize only %lu.\n", hv, Tablesize);
      exit(1);
   }
   /* Check to see if node is in the table. */
   if (connstate == Table[hv].cstate && 
       (memcmp(geneid,Table[hv].nodeid,sizeof(int)*NID)==0) ) {
      /*Hashanal[dumidx++] = dumcnt;*/
      /*printf("found in table.\n");*/
      /*if (flag) {
	printf("Found in table:\n");
	printclasses(Table[hv].nclasses, Table[hv].class);
	printf("prob = %lf\n", Table[hv].prob);
	}
      */
      for (i=0; i<NCOEF; i++)
	 prob->istate[i] = Table[hv].prob.istate[i];
      /*if (id0 != id1) {
	 printf("found in hashfind:\n");
	 printibdgr(node);
	 printpvec(prob);
	 }*/
      return 1;
   }
   else { /* Not in table so check collision list. */
      cell_p = Table[hv].next;
      if (cell_p) {
	 cst = cell_p->cstate;
	 nid = cell_p->nodeid;
      }
      for ( ; cell_p != NULL; cell_p = cell_p2, cst = cst2, nid = nid2) {
	 cell_p2 = cell_p->next; //Start getting info for next iteration.
	 if (cell_p2) {
	    cst2 = cell_p2->cstate;
	    nid2 = cell_p2->nodeid;
	 }

	 if (connstate == cst) {
	    if (memcmp(geneid, cell_p->nodeid, sizeof(int)*NID) == 0) {
	    //if (key == nid) {
	       for (i=0; i<NCOEF; i++)
		  prob->istate[i] = cell_p->prob.istate[i];
	       return 1;
	    }
	 }
      }

   }
   /*Hashanal[dumidx++] = -dumcnt;*/
   return 0;
}


void hashdone(void)
{
   long i, len;
   int nempty = 0;
   struct hash_cell *cell_p;

   FILE *out;
   out = fopen("hashanal.output", "w");

   for (i=0; i<Tablesize; i++) {
      if (Table[i].isempty) nempty++;
      else {
	 for (len=0, cell_p = Table[i].next; cell_p != NULL; 
	      cell_p = cell_p->next, len++) {;}
	 fprintf(out, "%ld\n", len);
      }
   }
   fclose(out);

   printf("nempty = %d\n", nempty);
   printf("nfull = Tablesize - nempty = %ld\n", Tablesize - nempty);
   printf("Ncells = %lu\n", Ncells);
   printf("Nmade = %d\n", Nmade);
   printf("Recycle = %d\n", Recycle);

   for (i=0; i<15; i++)
      printf("Call[%ld] = %u\n", i, Call[i]);

   /*out = fopen("hashanal.output", "w");
     for (i=0; i<ANALSZ; i++)
     if (Hashanal[i])
     fprintf(out,"%d\n", Hashanal[i]);
     fclose(out);*/
}

/* 
 * The FNV hash is originally due to Glenn Fowler, Landon Curt Noll
 * and Phong Vo.  The code for it has been placed in the public
 * domain. (http://isthe.com/chongo/tech/comp/fnv/). It is included
 * here as an alternative to the bbhash. In my testing the bbhash was,
 * on average, slightly faster, although this may not always hold
 * true.
 */

#define FNV_32_PRIME ((unsigned)0x01000193)
unsigned fnv_32(void *buf, size_t len, unsigned hval)
{
    unsigned char *bp = (unsigned char *)buf;	/* start of buffer */
    unsigned char *be = bp + len;		/* beyond end of buffer */

    /*
     * FNV-1a hash each octet in the buffer
     */
    while (bp < be) {

	/* xor the bottom with the current octet */
	hval ^= (unsigned)*bp++;

	/* multiply by the 32 bit FNV magic prime mod 2^32 */
#if defined(NO_FNVOPT)
	hval *= FNV_32_PRIME;
#else
	hval += (hval<<1) + (hval<<4) + (hval<<7) + (hval<<8) + (hval<<24);
#endif
    }

    /* return our new hash value */
    return hval;
}

/*
 * The following hash functions are written by Bob Jenkins (http://burtleburtle.net/bob/hash/index.html)
 * They have been placed in the public domain and I gratefully make use of them here.
 */

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bit set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k     : the key (the unaligned variable-length array of bytes)
  len   : the length of the key, counting by bytes
  level : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 36+6len instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burlteburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

ub4 bbhash(void *key, size_t length, unsigned initval)
//register ub1 *k;        /* the key */
//register ub4  length;   /* the length of the key */
//register ub4  initval;    /* the previous hash, or an arbitrary value */
{
   register ub4 a,b,c,len;
   register ub1 *k = (ub1 *)key;

   /* Set up the internal state */
   len = (ub4)length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;           /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((ub4)k[1]<<8) +((ub4)k[2]<<16) +((ub4)k[3]<<24));
      b += (k[4] +((ub4)k[5]<<8) +((ub4)k[6]<<16) +((ub4)k[7]<<24));
      c += (k[8] +((ub4)k[9]<<8) +((ub4)k[10]<<16)+((ub4)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((ub4)k[10]<<24);
   case 10: c+=((ub4)k[9]<<16);
   case 9 : c+=((ub4)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((ub4)k[7]<<24);
   case 7 : b+=((ub4)k[6]<<16);
   case 6 : b+=((ub4)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((ub4)k[3]<<24);
   case 3 : a+=((ub4)k[2]<<16);
   case 2 : a+=((ub4)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

/*
--------------------------------------------------------------------
 This works on all machines.  hash2() is identical to hash() on 
 little-endian machines, except that the length has to be measured
 in ub4s instead of bytes.  It is much faster than hash().  It 
 requires
 -- that the key be an array of ub4's, and
 -- that all your machines have the same endianness, and
 -- that the length be the number of ub4's in the key
--------------------------------------------------------------------
*/
ub4 bbhash2( //k, length, initval)
	  void *key,        /* the key */
	  register ub4  length,   /* the length of the key, in ub4s */
	  register ub4  initval  /* the previous hash, or an arbitrary value */
	  )
{
   register ub4 *k = (ub4 *)key;
   register ub4 a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;           /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 3)
   {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a,b,c);
      k += 3; len -= 3;
   }

   /*-------------------------------------- handle the last 2 ub4's */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
     /* c is reserved for the length */
   case 2 : b+=k[1];
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

/*
--------------------------------------------------------------------
mix -- mix 3 64-bit values reversibly.
mix() takes 48 machine instructions, but only 24 cycles on a superscalar
  machine (like Intel's new MMX architecture).  It requires 4 64-bit
  registers for 4::2 parallelism.
All 1-bit deltas, all 2-bit deltas, all deltas composed of top bits of
  (a,b,c), and all deltas of bottom bits were tested.  All deltas were
  tested both on random keys and on keys that were nearly all zero.
  These deltas all cause every bit of c to change between 1/3 and 2/3
  of the time (well, only 113/400 to 287/400 of the time for some
  2-bit delta).  These deltas all cause at least 80 bits to change
  among (a,b,c) when the mix is run either forward or backward (yes it
  is reversible).
This implies that a hash using mix64 has no funnels.  There may be
  characteristics with 3-bit deltas or bigger, I didn't test for
  those.
--------------------------------------------------------------------
*/
#define mix64(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>43); \
  b -= c; b -= a; b ^= (a<<9); \
  c -= a; c -= b; c ^= (b>>8); \
  a -= b; a -= c; a ^= (c>>38); \
  b -= c; b -= a; b ^= (a<<23); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>35); \
  b -= c; b -= a; b ^= (a<<49); \
  c -= a; c -= b; c ^= (b>>11); \
  a -= b; a -= c; a ^= (c>>12); \
  b -= c; b -= a; b ^= (a<<18); \
  c -= a; c -= b; c ^= (b>>22); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 64-bit value
  k     : the key (the unaligned variable-length array of bytes)
  len   : the length of the key, counting by bytes
  level : can be any 8-byte value
Returns a 64-bit value.  Every bit of the key affects every bit of
the return value.  No funnels.  Every 1-bit and 2-bit delta achieves
avalanche.  About 41+5len instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 64 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, Jan 4 1997.  bob_jenkins@burtleburtle.net.  You may
use this code any way you wish, private, educational, or commercial,
but I would appreciate if you give me credit.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^64
is acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

/*
--------------------------------------------------------------------
 This works on all machines, is identical to hash() on little-endian 
 machines, and it is much faster than hash(), but it requires
 -- that the key be an array of ub8's, and
 -- that all your machines have the same endianness, and
 -- that the length be the number of ub8's in the key
--------------------------------------------------------------------
*/
ub8 bbhash64( // k, length, level)
	     register ub8 *k,        /* the key */
	     register ub8  length,   /* the length of the key */
	     register ub8  level    /* the previous hash, or an arbitrary value */
	     )
{
  register ub8 a,b,c,len;

  /* Set up the internal state */
  len = length;
  a = b = level;                         /* the previous hash value */
  c = 0x9e3779b97f4a7c13LL; /* the golden ratio; an arbitrary value */

  /*---------------------------------------- handle most of the key */
  while (len >= 3)
  {
    a += k[0];
    b += k[1];
    c += k[2];
    mix64(a,b,c);
    k += 3; len -= 3;
  }

  /*-------------------------------------- handle the last 2 ub8's */
  c += (length<<3);
  switch(len)              /* all the case statements fall through */
  {
    /* c is reserved for the length */
  case  2: b+=k[1];
  case  1: a+=k[0];
    /* case 0: nothing left to add */
  }
  mix64(a,b,c);
  /*-------------------------------------------- report the result */
  return c;
}
