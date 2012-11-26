#Makes the version of the converter from action-angle vars to cartesian
#that works with MathLink

CC = gcc -g

# for galileo
MATHEMATICA_INSTALL_DIR = /Software/users/math/Mathematica_8.0.1
MLINKDIR = ${MATHEMATICA_INSTALL_DIR}/SystemFiles/Links/MathLink/DeveloperKit/
SYS = Linux-x86-64
MLIBS = -lML64i3 -lpthread -lrt -lstdc++

#for my laptop
#MATHEMATICA_INSTALL_DIR = /Applications/Mathematica.app/
#MLINKDIR = ${MATHEMATICA_INSTALL_DIR}/SystemFiles/Links/MathLink/DeveloperKit/
#SYS = 
#MLIBS = -lMLi3 -lstdc++ -framework Foundation

CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions

INCDIR = ${CADDSDIR}
LIBDIR = ${CADDSDIR}

MPREP = "${CADDSDIR}/mprep"
MCC = "${CADDSDIR}/mcc"


EXTRALIBS = -lm -lgsl -lgslcblas

AATOCART = aa_to_cart
ALLMATH = main_math.o ${AATOCART}tm.o ${AATOCART}.o eta_from_aa.o 
ALLC = main_c.o ${AATOCART}.o eta_from_aa.o

all: ${AATOCART}_math ${AATOCART}_c

${AATOCART}_c: ${ALLC}
	${CC} ${ALLC} ${EXTRALIBS} -o $@

${AATOCART}_math: ${ALLMATH}
	${CC} -I${INCDIR} ${ALLMATH} -L${LIBDIR} ${MLIBS} ${EXTRALIBS} -o $@ 

eta_from_aa.o: eta_from_aa.c
	${CC} -c $?

${AATOCART}.o: ${AATOCART}.c
	${CC} -c $? 

${AATOCART}tm.o:  ${AATOCART}.tm
	${MPREP} $? -o ${AATOCART}tm.c
	${CC} -c -I${INCDIR} ${AATOCART}tm.c 

main_c.o: main_c.c
	${CC} -c $?

main_math.o: main_math.c
	${CC} -c $? -I${INCDIR}

clean: 
	rm -f ${ALLMATH} ${AATOCART}_math main_c.o ${AATOCART}_c *~