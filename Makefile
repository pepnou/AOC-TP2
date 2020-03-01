CC=gcc

CFLAGS=-std=c99 -g3 -g

OFLAGS=-O3 -march=skylake -funroll-loops

LFLAGS=-lm -pthread

all: diffusion_X11 diffusion

diffusion_X11:
	$(CC) $(CFLAGS) $(OFLAGS) flame.c -DWITH_X11=1 diffusion.c -o diffusion_X11 $(LFLAGS) -lX11

diffusion:
	$(CC) $(CFLAGS) $(OFLAGS) diffusion.c -o diffusion $(LFLAGS)

clean:
	rm -Rf *~ diffusion_X11 diffusion
