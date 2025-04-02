CC?=cc
CFLAGS?=-O3

VPATH=lib/quadsort
ALL=bf-pbwt 2bfpbwt
LIBOMP?=/opt/homebrew/opt/libomp

debug: CFLAGS=-O0 -g

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

bf-pbwt: bf-pbwt.o
	${CC} $^ -o $@

2bfpbwt: 2bfpbwt.c
	${CC} -o $@ ${CFLAGS} -I ${LIBOMP}/include -L ${LIBOMP}/lib -lomp  $^ 

gen: gen.c
	${CC} -O3 $^ -o $@

pgen: pgen.c
	${CC} -O3 $^ -o $@

rrsort_test: rrsort_test.c
	${CC} -O3 $^ -o $@ -march=native

clean:
	-/bin/rm -f *.o

remove: clean
	-/bin/rm -i $(ALL)
