CC?=cc
CFLAGS?=-O3

VPATH=lib/quadsort
ALL=2bfpbwt pgen
LIBOMP?=/opt/homebrew/opt/libomp

APPLE_CLANG:=$(shell $(CC) --version | grep "Apple clang" > /dev/null && echo 1 || echo 0)

ifeq ($(APPLE_CLANG),1)
    CFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
else
    CFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

debug: CFLAGS=-O0 -g

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

2bfpbwt: 2bfpbwt.c
	${CC} -o $@ ${CFLAGS} -I ${LIBOMP}/include -L ${LIBOMP}/lib $(LDFLAGS) $^ 

pgen: pgen.c
	${CC} -O3 $^ -o $@

clean:
	-rm *.o $(ALL)
