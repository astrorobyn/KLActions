CC = gcc

#set the location of the GSL libraries and header files if necessary
#GSLLIBDIR =
#GSLINCDIR =

#compile and link flags if GSL location must be specified
#CFLAGS = -I${GSLINCDIR}
#LDFLAGS = -L${GSLLIBDIR} -lm -lgsl -lgslcblas

#if GSL is in the path use these lines instead
CFLAGS =
LDFLAGS =  -lm -lgsl -lgslcblas

# Location of mathlink stuff for galileo
MATHEMATICA_INSTALL_DIR = /Software/users/math/Mathematica_8.0.1
MLINKDIR = ${MATHEMATICA_INSTALL_DIR}/SystemFiles/Links/MathLink/DeveloperKit/
SYS = Linux-x86-64
MLIBS = -lML64i3 -lpthread -lrt -lstdc++

# Location of mathlink stuff for my laptop
#MATHEMATICA_INSTALL_DIR = /Applications/Mathematica.app/
#MLINKDIR = ${MATHEMATICA_INSTALL_DIR}/SystemFiles/Links/MathLink/DeveloperKit/
#SYS = 
#MLIBS = -lMLi3 -lstdc++ -framework Foundation

CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = "${CADDSDIR}/mprep"
MCC = "${CADDSDIR}/mcc"

#if GSL is in the path use this one
CFLAGS = -g
LDFLAGS =  -g -lm -lgsl -lgslcblas

#names of executable and objects
EXEC = cartToAA
ALLC = main.o isochroneFunctions.o bodyAux.o io.o 
ALLMATH = main_math.o ${EXEC}tm.o bodyAux.o isochroneFunctions.o

all: ${EXEC}_math ${EXEC}_c

${EXEC}_c: ${ALLC}
	${CC} ${ALLC} ${LDFLAGS} -o $@

${EXEC}_math: ${ALLMATH}
	${CC} ${ALLMATH} ${LDFLAGS} -L${LIBDIR} ${MLIBS} -o $@ 

main_math.o: main_math.c
	${CC} -c -I${INCDIR} ${CFLAGS} $? -o $@

${EXEC}tm.o: ${EXEC}tm.c
	${CC} -c -I${INCDIR} -g ${EXEC}tm.c 

${EXEC}tm.c: ${EXEC}.tm
	${MPREP} $? -o ${EXEC}tm.c

isochroneFunctions.o: isochroneFunctions.c isochroneFunctions.h
	$(CC) -c ${CFLAGS} isochroneFunctions.c -o $@

bodyAux.o: bodyAux.c bodyAux.h
	$(CC) -c ${CFLAGS} bodyAux.c -o $@

io.o: io.c io.h
	$(CC) -c ${CFLAGS} io.c -o $@

clean: 
	rm -f ${ALLC} ${ALLMATH} ${EXEC} cartToAA.tm.* *~