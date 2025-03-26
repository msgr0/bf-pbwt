CC?=cc
CFLAGS?=-O3

VPATH=lib/quadsort
ALL=bf-pbwt 2bfpbwt

debug: CFLAGS=-O0 -g

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

bf-pbwt: bf-pbwt.o
	${CC} $^ -o $@

2bfpbwt: 2bfpbwt.o
	${CC} $< -o $@
	
pgen: pgen.c
	${CC} -O3 $^ -o $@

clean:
	-/bin/rm -f *.o

remove: clean
	-/bin/rm -i $(ALL)
