CC?=gcc
CFLAGS?=-O3

ALL=plot
LIBPLPLOT=/opt/homebrew/include/plplot
HBLIB=/opt/homebrew/lib


APPLE_CLANG:=$(shell $(CC) --version | grep "Apple clang" > /dev/null && echo 1 || echo 0)

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

plot: plots.c
	${CC} -o $@ ${CFLAGS} -I${LIBPLPLOT} -L${HBLIB} -lplplot $^ 

