CC = g++ -std=gnu++11 

DEBUG = -g

CFLAGS = -c $(DEBUG) 

LFLAGS = $(DEBUG)

OPTIMAZATION=  -O3 -fopenmp 

OBJS = energy_comp.o particle_func.o plot_func.o minimize.o generate_metropolis.o test.o

PROGRAM_NAME = program

program: $(OBJS)
	$(CC)  $(LFLAGS) $(OBJS) -o $(PROGRAM_NAME)

test.o: test.cpp allheaders.h
	$(CC) $(CFLAGS) test.cpp 



energy_comp.o: energy_comp.cpp allheaders.h
	$(CC)  $(CFLAGS) energy_comp.cpp 

plot_func.o:plot_func.cpp allheaders.h
	$(CC) $(CFLAGS) plot_func.cpp

particle_func.o:particle_func.cpp allheaders.h
	$(CC) $(CFLAGS) particle_func.cpp

minimize.o:minimize.cpp montecarlo.h
	$(CC) $(CFLAGS) minimize.cpp

generate_metropolis.o:generate_metropolis.cpp montecarlo.h
	$(CC) $(CFLAGS) generate_metropolis.cpp




clean:
	\rm *.o  $(PROGRAM_NAME)


run :  $(PROGRAM_NAME)
	./$(PROGRAM_NAME)

