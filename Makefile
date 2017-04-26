CC = gfortran
CFLAG = -fbacktrace -c -g -ffree-line-length-500 
CCFLAG= -g -fbacktrace
CLIB = -llapack -lblas 

all : main


main : types.o constants.o input.o output.o utils.o interpol1D.o main2.o 
	$(CC) types.o constants.o input.o output.o utils.o interpol1D.o main2.o -o main $(CCFLAG) $(CLIB)

main2.o : main2.f90
	$(CC) main2.f90 $(CFLAG) $(CLIB)
	
types.o : types.f90
	$(CC) types.f90 $(CFLAG) $(CLIB)
	
constants.o : constants.f90
	$(CC) constants.f90 $(CFLAG) $(CLIB)

input.o : input.f90
	$(CC) input.f90 $(CFLAG) $(CLIB)
	
output.o : output.f90
	$(CC) output.f90 $(CFLAG) $(CLIB)
	
utils.o : utils.f90
	$(CC) utils.f90 $(CFLAG) $(CLIB)
	
interpol1D.o : interpol1D.f90
	$(CC) interpol1D.f90 $(CFLAG) $(CLIB)
	

clean :
	rm *.o *.mod main

	
