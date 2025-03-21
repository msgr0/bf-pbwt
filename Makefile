CC = gcc

CFLAGS = -O0 -g

ALL = bf-pbwt

%.0: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

bf-pbwt: bf-pbwt.o
	${CC} $^ -o $@

clean:
	-/bin/rm -f *.o

remove: clean
	-/bin/rm -i $(ALL)
