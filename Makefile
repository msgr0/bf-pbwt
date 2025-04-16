CC?=cc
CFLAGS?=-O3

VPATH=lib/quadsort
ALL=2bfpbwt-bm gen
ALL+=2bfpbwt # NOTE: this will be removed
LIBOMP?=/opt/homebrew/opt/libomp

APPLE_CLANG:=$(shell $(CC) --version | grep "Apple clang" > /dev/null && echo 1 || echo 0)

ifeq ($(APPLE_CLANG),1)
    CFLAGS += -Xpreprocessor -fopenmp
    LDFLAGS += -lomp
else
    CFLAGS += -fopenmp
    LDFLAGS += -fopenmp
endif

.PHONY: all
all: ${ALL}

debug: CFLAGS=-O0 -g

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@


# NOTE: this will be removed
2bfpbwt: 2bfpbwt.c iobm.o
	${CC} -o $@ ${CFLAGS} -I ${LIBOMP}/include -L ${LIBOMP}/lib $(LDFLAGS) $^ 

2bfpbwt-bm: 2bfpbwt.c iobm.o
	${CC} -o $@ ${CFLAGS} -I ${LIBOMP}/include -L ${LIBOMP}/lib $(LDFLAGS) $^ 
gen: gen.c
	${CC} -O3 $^ -o $@
	
clean:
	-rm *.o $(ALL)


SHELL:=/bin/bash
PANEL?=panel.r5k.c50k.txt
LOOPS?=10
MODES?=bar bli blim
loop-%: 2bfpbwt
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
