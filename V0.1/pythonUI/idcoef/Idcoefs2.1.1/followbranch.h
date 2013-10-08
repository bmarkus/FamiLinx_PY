/* followbranch5.h */

Probvec_t nodeprob(Ibdgraph_t *genenode);
Probvec_t nodeprob_sym(Ibdgraph_t *genenode);
void sortindex(int nelem, int *array, int *idxar);
int arerepeatedfounders(Ibdgraph_t *node, int *slot);
int onlyfounders(int *id);
