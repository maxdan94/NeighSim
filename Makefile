CC=gcc
CFLAGS=-O9

all: sim cosine jaccard rmhub

sim : sim.c
	$(CC) $(CFLAGS) sim.c -o sim -lm -fopenmp

cosine : cosine_opt.c
	$(CC) $(CFLAGS) cosine_opt.c -o cosine_opt -lm -fopenmp

jaccard : jaccard_opt.c
	$(CC) $(CFLAGS) jaccard_opt.c -o jaccard_opt -lm -fopenmp

rmhub : rmhub.c
	$(CC) $(CFLAGS) rmhub.c -o rmhub

clean:
	rm sim cosine_opt jaccard_opt rmhub
