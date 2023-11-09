CC=gcc
CFLAGS=-Wextra -Wall -lm
EXEC=RO

all: $(EXEC)

RO: TP1.o functions.o
	$(CC) -o $@ $^ $(CFLAGS)
	@echo "\nRun the program with ./RO"

TP1.o: TP1.c
	$(CC) -o $@ -c $< $(CFLAGS)

functions.o: TP1Functions.c TP1Functions.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean: 
	rm -f *.o

rmExec: clean
	rm -f $(EXEC)