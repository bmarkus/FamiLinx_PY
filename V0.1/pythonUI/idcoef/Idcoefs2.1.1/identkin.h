/* identkin.h */

#define NIDSTATE 9

void kin2ident(double *ident, Probvec_t *pvec);
double kindiff(Probvec_t *x1, Probvec_t *x2);
