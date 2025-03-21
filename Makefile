CC=cc
CFLAGS=-O0 -g

ALL = bf-pbwt

%.o: %.c %.h
	${CC} -c ${CFLAGS} $< -o $@

all: ${ALL}

bf-pbwt: bf-pbwt.o
	${CC} $^ -o $@

2bfpbwt: 2bfpbwt.o
	${CC} $^ -o $@
	

clean:
	-/bin/rm -f *.o

remove: clean
	-/bin/rm -i $(ALL)
