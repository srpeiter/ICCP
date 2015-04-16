CC = g++ 

DEBUG = -g

CFLAGS = -c $(DEBUG) 

LFLAGS = $(DEBUG)

OPTIMAZATION=  -O3 -fopenmp 

OBJS = particle_func.o minimize.o main.o

PROGRAM_NAME = program

program: $(OBJS)
	$(CC)  $(LFLAGS) $(OBJS) -o $(PROGRAM_NAME)

test.o: main.cpp allheaders.h
	$(CC) $(CFLAGS) main.cpp 

particle_func.o:particle_func.cpp allheaders.h
	$(CC) $(CFLAGS) particle_func.cpp

minimize.o:minimize.cpp montecarlo.h
	$(CC) $(CFLAGS) minimize.cpp



clean:
	\rm *.o $(PROGRAM_NAME)


run :  $(PROGRAM_NAME)
	./$(PROGRAM_NAME)
	python data.py



	

