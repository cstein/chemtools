CC=g++
OB_BASE=/home/casper/Programs
OB_INCLUDE=-I${OB_BASE}/include/openbabel-2.0/
OB_LIBRARY=-L${OB_BASE}/lib -lopenbabel

INCLUDE=${OB_INCLUDE}
LIB=${OB_LIBRARY}

CPPFLAGS=-Wextra -Wall

opthydrogens: main.o
	${CC} $^ ${LIB} -o $@

main.o: main.cpp
	${CC} ${CPPFLAGS} ${INCLUDE} -c $<

clean:
	rm -f opthydrogens *.o
