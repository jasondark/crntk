CC := gcc
CFLAGS := -fopenmp -std=c99 -pedantic -Wall -mtune=native -march=native -O2

obj/crntk.o: src/crntk.c include/crntk.h
	$(CC) $(CFLAGS) -c src/crntk.c -o obj/crntk.o

example: src/example.c obj/crntk.o
	$(CC) $(CFLAGS) src/example.c obj/crntk.o -lm -o bin/example

clean:
	rm -f obj/* && rm -f bin/*

