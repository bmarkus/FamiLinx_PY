/* pedigree.h - Version 5 */
#define PEDH

#define MAXLEN 1024

/* General header file for pedigree structures. */

/**#define ISAFOUNDER(id) ((id)<0 ? 1 : (getparent(0,id) || getparent(1,id)) ? \
 **			0 : 1 )
 **/
#define ISAFOUNDER(id) ((id)<0 ? 1 : isafounder(id))

int getparent(int whichpar, int indiv);
int getnpeop(void);
void readped(char *pedfile/*, char *sampfile*/ );
int findid(int fid);
int isafounder(int id);
void minimalped(int nsample, int *samplelist);
void printminped(void);
int getfindiv(int id);
