/* main.c */

/* Parse the command line arguments for input and output file names,
   then do the computations. The following command line arguments are
   required:
   -o outputfile
   -p pedigree file
   -s study sample file
   -r RAM available for the computations (in MB)
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
#include "pedigree.h"
#include "ibdgraph.h"
#include "followbranch.h"
#include "nodehash.h"
#include "compute.h"

int main (int argc, char **argv)
{
   char pedfile[MAXLEN] = "pedigree";
   char samplefile[MAXLEN] = "study.sample";
   char output[MAXLEN] = "output";
   int ramsize=100; /* Available RAM in MB */
   int arg;


   if (argc > 1) {
      for (arg=1; arg < argc && argv[arg][0] == '-'; arg++) {
         switch (argv[arg][1]) 
            {
            case 'r':
               ramsize = atoi(argv[++arg]);  
               printf("ramsize = %d\n", ramsize);
               break;
            case 'p':
               strncpy(pedfile, argv[++arg], MAXLEN);
               printf("pedigree file: %s\n", pedfile);
               fflush(stdout);
               break;
            case 's':
               strncpy(samplefile, argv[++arg], MAXLEN);
               printf("sample file: %s\n", samplefile);
               break;
            case 'o':
               strncpy(output, argv[++arg], MAXLEN);
               printf("output file: %s\n", output);
               break;
	    case 'h':
	      fprintf(stdout, "Options are:\n");
	      fprintf(stdout, "-p filename (pedigree file)\n");
	      fprintf(stdout, "-s filename (study sample file)\n");
	      fprintf(stdout, "-o filename (output file)\n");
	      fprintf(stdout, "-r number (available RAM in megabytes)\n");
	      exit(0);
            default:
               fprintf (stderr, "Unknown option \"%s\"\n", argv[arg]);
               exit(1);
            }
      }
   }

   readped(pedfile);
   readsample(samplefile);
   hashinit(ramsize, getnpeop());

   computeidcoefs(output);

#  if defined(DEBUG) || defined(DEBUG_HASH)
   hashdone();
#  endif

   printf("Successful completion.\n");
   return 0;
}
