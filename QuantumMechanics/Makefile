CC = g++ -std=c++11

DEBUG = -g

CFLAGS = -c $(DEBUG) 

LFLAGS = $(DEBUG)

OPTIMAZATION=  -O3 -fopenmp 

OBJS = energy_comp.o particle_func.o plot_func.o test.o

PROGRAM_NAME = program

program: $(OBJS)
	$(CC)  $(LFLAGS) $(OBJS) -o $(PROGRAM_NAME)

test.o: test.cpp allheaders.h
	$(CC) $(OPTIMAZATION) $(CFLAGS) test.cpp 



energy_comp.o: energy_comp.cpp allheaders.h
	$(CC) $(OPTIMAZATION) $(CFLAGS) energy_comp.cpp 

plot_func.o:plot_func.cpp allheaders.h
	$(CC) $(OPTIMAZATION) $(CFLAGS) plot_func.cpp

particle_func.o:particle_func.cpp allheaders.h
	$(CC) $(OPTIMAZATION) $(CFLAGS) particle_func.cpp


clean:
	\rm *.o  $(PROGRAM_NAME)


run :  $(PROGRAM_NAME)
	./$(PROGRAM_NAME)

