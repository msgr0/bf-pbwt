CC?=cc
CFLAGS?=-O3

VPATH=lib/quadsort
ALL=gen 2bfpbwt-bm 2bfpbwt-bcf
LIBOMP?=/opt/homebrew/opt/libomp
HTSLIB?=/opt/htslib

LIBOMP_INCL=-I ${LIBOMP}/include -L ${LIBOMP}/lib
HTSLIB_INCL=-I ${HTSLIB}/include -L ${HTSLIB}/lib

APPLE_CLANG:=$(shell $(CC) --version | grep "Apple clang" > /dev/null && echo 1 || echo 0)
ifeq ($(APPLE_CLANG),1)
    CFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
else
    CFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

# CCINCL+= -I ${HTSLIB}/include -L ${HTSLIB}/lib

.PHONY: all
all: ${ALL}

debug: CFLAGS=-fstrict-aliasing -Wstrict-aliasing -fsanitize=address -O0 -g
debug: ${ALL}

leaks: CFLAGS=-fstrict-aliasing -Wstrict-aliasing -O0 -g
leaks: ${ALL}

%.o: %.c %.h
	${CC} -c ${CFLAGS} ${CCINCL} $< -o $@

2bfpbwt-enc: CCINCL=${LIBOMP_INCL}
2bfpbwt-enc: 2bfpbwt.c ioenc.o
	${CC} -o $@ ${CFLAGS} -DBF2IOMODE_ENC ${CCINCL} $(LDFLAGS) $^ 

2bfpbwt-bm: CCINCL=${LIBOMP_INCL}
2bfpbwt-bm: 2bfpbwt.c iobm.o
	${CC} -o $@ ${CFLAGS} -DBF2IOMODE_BM ${CCINCL} $(LDFLAGS) $^ 

iobcf.o: iobcf.c
	${CC} -c ${CFLAGS} ${HTSLIB_INCL} $< -o $@ -lhts

2bfpbwt-bcf: CCINCL=${LIBOMP_INCL} ${HTSLIB_INCL}
2bfpbwt-bcf: 2bfpbwt.c iobcf.o
	${CC} -o $@ ${CFLAGS} -DBF2IOMODE_BCF ${CCINCL} $(LDFLAGS) -lhts $^ 

gen: gen.c
	${CC} -O3 $^ -o $@

bcf2enc: CCINCL=${HTSLIB_INCL}
bcf2enc: bcftoenc.c 
	${CC} -o $@ ${CCINCL} -lhts $^

clean:
	-rm *.o $(ALL)


CHECKPANEL?=panel.r10.c800.txt
CHECKMODES?=lin bli blis blim ars bar bars barm prs bpr spr
check-bm-%: 2bfpbwt-bm
	@./$^ $* ${CHECKPANEL} 1>/dev/null 2>&1
check-bm: $(foreach m,$(CHECKMODES),check-bm-$(m))
.PHONY: check-bm

SHELL:=/bin/bash
PANEL?=panel.r5k.c50k.txt
LOOPS?=10
MODES?=bar bli blim
loop-%: 2bfpbwt-bm
	@tmpfile=$$(mktemp); \
	for ((i=1; i <= ${LOOPS}; ++i)) do \
		echo -e "[$*] running $$i" >&2; \
		\time -al ./$^ $* ${PANEL} 2>&1 \
			| grep "real" | tr -s ' ' | cut -f2 -d' ' >> $$tmpfile ;\
	done; \
	sum=$$(paste -sd+ $$tmpfile | bc); \
	avg=$$(echo "scale=4; $$sum / ${LOOPS}" | bc); \
	echo "[$*] avg_time: $$avg s"; \

loop: $(foreach m,$(MODES),loop-$(m))

.PHONY: loop
