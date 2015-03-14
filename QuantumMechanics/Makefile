CC = g++ -std=c++11

DEBUG = -g

CFLAGS = -c $(DEBUG)

LFLAGS = $(DEBUG)

OBJS = file.o plot_func.o

PROGRAM_NAME = program

program: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(PROGRAM_NAME)

file.o: file.cpp plotting.h
	$(CC) $(CFLAGS) file.cpp 

plot_func.o:plot_func.cpp plotting.h
	$(CC) $(CFLAGS) plot_func.cpp

clean:
	\rm *.o  program

run :  $(PROGRAM_NAME)
	./$(PROGRAM_NAME)

