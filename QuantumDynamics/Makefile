CC = g++ -std=c++11

LIBS = -lblas -llapack -lsuperlu

DEBUG = -g

CFLAGS = -c $(DEBUG) 

LFLAGS = $(DEBUG)

OPTIMAZATION=  -O3 -fopenmp 

OBJS = main.o solver.o solver2d.o

PROGRAM_NAME = program

program: $(OBJS)
	$(CC)  $(LFLAGS) $(OBJS)  -o $(PROGRAM_NAME) $(LIBS)

main.o: main.cpp solver.h 
	$(CC) $(CFLAGS) main.cpp 

solver.o: solver.cpp solver.h math_module.h
	$(CC) $(CFLAGS) solver.cpp


solver2d.o: solver2d.cpp solver.h math_module.h
	$(CC) $(CFLAGS) solver2d.cpp




clean:
	\rm *.o $(PROGRAM_NAME) *.dat *.txt


run :  $(PROGRAM_NAME)
	./$(PROGRAM_NAME)
	
