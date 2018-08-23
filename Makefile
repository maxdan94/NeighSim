CC=gcc
CFLAGS=-O9

all: sim sim2 cosine jaccard jaccard2 rmhub

sim : sim.c
	$(CC) $(CFLAGS) sim.c -o sim -lm -fopenmp

sim2 : sim_nohub.c
	$(CC) $(CFLAGS) sim_nohub.c -o sim_nohub -lm -fopenmp

cosine : cosine_opt.c
	$(CC) $(CFLAGS) cosine_opt.c -o cosine_opt -lm -fopenmp

jaccard : jaccard_opt.c
	$(CC) $(CFLAGS) jaccard_opt.c -o jaccard_opt -lm -fopenmp

jaccard2 : jaccard_opt_nohub.c
	$(CC) $(CFLAGS) jaccard_opt_nohub.c -o jaccard_opt_nohub -lm -fopenmp

rmhub : rmhub.c
	$(CC) $(CFLAGS) rmhub.c -o rmhub

clean:
	rm sim sim_nohub cosine_opt jaccard_opt jaccard_opt_nohub rmhub
