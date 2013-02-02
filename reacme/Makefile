CC=g++
OB_BASE=/home/casper/Programs
OB_INCLUDE=-I${OB_BASE}/include/openbabel-2.0
OB_LIBRARY=-L${OB_BASE}/lib -lopenbabel

INCLUDE=${OB_INCLUDE}
LIB=${OB_LIBRARY}

CPPFLAGS=-Wextra -Wall

reacme: main.o reactor.o
	${CC} $^ ${LIB} -o $@

main.o: main.cpp
	${CC} ${CPPFLAGS} ${INCLUDE} -c $<

reactor.o: reactor.cpp reactor.hpp
	${CC} ${CPPFLAGS} ${INCLUDE} -c $<

clean:
	rm -f reacme *.o
	rm -rf html latex
