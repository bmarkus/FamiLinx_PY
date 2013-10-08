CC = gcc

CCFLAGS = -g -O3
GCC33FLAGS = -g -O3 -fno-strength-reduce

CODE = main.c pedigree.c compute.c
HEADERS = pedigree.h compute.h

bugcheck: $(CODE)
	$(CC) $(CCFLAGS) -o bugcheck $(CODE) 

idcoefs33: $(CODE) $(HEADERS)
	$(CC) $(GCC33FLAGS) -o idcoefs $(CODE) -lm

debug: $(CODE) $(HEADERS)
	$(CC) $(CCFLAGS) -o idcoefs.debug -D DEBUG $(CODE) -lm

64: $(CODE) $(HEADERS)
	$(CC) -m64  $(CCFLAGS) -o idcoefs $(CODE) -lm
