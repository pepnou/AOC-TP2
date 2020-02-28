CC=gcc

CFLAGS=-std=c99 -g3

OFLAGS=-O0

LFLAGS=-lm

all: diffusion_X11 diffusion

diffusion_X11:
	$(CC) $(CFLAGS) $(OFLAGS) flame.c -DWITH_X11=1 diffusion.c -o diffusion_X11 $(LFLAGS) -lX11

diffusion:
	$(CC) $(CFLAGS) $(OFLAGS) diffusion.c -o diffusion $(LFLAGS)

clean:
	rm -Rf *~ diffusion_X11 diffusion