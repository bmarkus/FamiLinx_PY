/* nodehash.h */

/* Functions needed to use a hash to store and recall computed probabilities
   for sets of equivalence classes.
*/


void hashinit(int size, int npeop);
void hashstore(Ibdgraph_t *node, Probvec_t *prob);
int hashfind(Ibdgraph_t *node, Probvec_t *prob);
void hashdone(void);
