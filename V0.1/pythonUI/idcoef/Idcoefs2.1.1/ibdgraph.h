/* ibdgraph.h */
#define IBDGRH


//#ifndef PEDH
//#include "pedigree.h"
//#endif


/* Functions for using ibdgraphs. */

#define NGENES 4 /* For generalized kinship coefs up to four people. */
#define NCOEF 15  /* 15 is number of detailed identity states for 4 genes.*/

typedef struct  {
   /*int coef;*//* all elements are multiplied by 1/2^coef to get prob of that state.*/
   double istate[NCOEF];
} Probvec_t;

typedef struct ibdgraph {
   int genelist[NGENES]; /* ID's of people in the list. */
   int connectarr[NGENES]; /* Identical entries are connected. */
   int connstate; // Which of the 15 possible connection states connectarr is in.
}Ibdgraph_t;

#define areconnected(c,i,j) ( (c)[i] == (c)[j] )

Ibdgraph_t ibdgr_init(int count, ...);
void boundcond4(Ibdgraph_t *node, Probvec_t *prob);
void newboundcond(Ibdgraph_t *node, Probvec_t *prob);
void printibdgr(Ibdgraph_t *graph);

//int areconnected(int *connectar, int idx1, int indx2);
void connected2parent(Ibdgraph_t *node, int prnt, int idx);
void founder2parent(Ibdgraph_t *node, int prnt, int idx);
void connect(Ibdgraph_t *node, int idx1, int idx2);
int connectstate(int *connect);

void initpvec(Probvec_t *pvec);
void add2pvec(Probvec_t *orig, Probvec_t *x);
void mult2pvec(double coef, Probvec_t *x);
void printpvec(Probvec_t *x);
double pvecdiff(Probvec_t *x1, Probvec_t *x2);
