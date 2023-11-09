CC=gcc
CFLAGS=-Wextra -Wall -lm
EXEC=RO

all: $(EXEC)

RO: main.o functions.o
	$(CC) -o $@ $^ $(CFLAGS)
	@echo "\nRun the program with ./RO"

main.o: main.c
	$(CC) -o $@ -c $< $(CFLAGS)

functions.o: Functions.c Functions.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean: 
	rm -f *.o

rmExec: clean
	rm -f $(EXEC)