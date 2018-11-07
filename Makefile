CFLAGS = -O3 -std=c99 
LDFLAGS = -lm
CC = gcc

integrator : integrator.c
	gcc -o integrator integrator.c -lm $(CFLAGS)
