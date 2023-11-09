CC=gcc
CFLAGS=-Wextra -Wall -lm
EXEC=RO

all: $(EXEC)

RO: RO.o functions.o
	$(CC) -o $@ $^ $(CFLAGS)
	@echo "\nRun the program with ./RO"

RO.o: RO.c
	$(CC) -o $@ -c $< $(CFLAGS)

functions.o: functions.c functions.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean: 
	rm -f *.o

rmExec: clean
	rm -f $(EXEC)