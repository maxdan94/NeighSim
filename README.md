## Info:

These C programs computes some neighborhood-based similarities between each pair of nodes in a very large sparse graph.

The cosine similarity between nodes $u$ and $v$ is defined as 
$$\frac{|\Delta(u) \cap \Delta(v)|}{\sqrt{|\Delta(u)|\cdot |\Delta(v)|}},$$ 
where $|\Delta(u)|$ is the set of neighbors of node u.

The jaccard similarity between nodes $u$ and $v$ is defined as 
$$\frac{|\Delta(u) \cap \Delta(v)|}{|\Delta(u) \cup \Delta(v)|},$$

The F1 similarity between nodes $u$ and $v$ is defined as 
$$\frac{2 |\Delta(u) \cap \Delta(v)|}{|\Delta(u)| + |\Delta(v)|},$$

In practice, the sim.c is quite scalable as it avoids to compute the cosine similarity between pairs of nodes having no neighbors in common.

"cosine_opt.c" and "jaccard_opt.c" are more scalable. They aim at computing each pair of nodes with similarity higher than a threshold $\alpha$ (input parameter). Using this threshold, computing the similarities between pairs of nodes that have too different degrees is avoided (as these pairs of nodes will have a similarity lower than the threshold).

## To compile:
type "Make", or type
- gcc sim.c -O3 -o sim -lm -fopenmp
- gcc cosine_opt.c -O3 -o cosine_opt -lm -fopenmp
- gcc jaccard_opt.c -O3 -o jaccard_opt -lm -fopenmp
- gcc rmhub.c -O3 -o rmhub

## To execute:

./sim p net.txt
- p is the number of threads to use (nearly optimal degree of parallelism)
- net.txt is the input directed graph "source target" on each line. Node's IDs should be integers, preferably from 0 to n-1. 
It will print values in the terminal to plot a histogram with 0.1 bucket-size.

./cosine_opt p a net.txt 
./jaccard_opt p a net.txt
- p is the number of threads to use (nearly optimal degree of parallelism)
- a is the input threshold : only similarities higher than this threshold
- net.txt is the input directed graph "source target" on each line. Node's IDs should be integers, preferably from 0 to n-1. 
It will print values in the terminal to plot a histogram with 0.1 bucket-size.


The programs will be faster if the input graph has small degrees. Indeed, the running time is in $O(\sum_{u\in V} d(u)^2)$.

If the program does not scale, because there are too many nodes with a very high degree, then just remove these hubs:

./rmhub max-degree neti.txt neto.txt

- max-degree is the maximum allowed degree. For instance 10,000.
- neti.txt is the input directed graph: "source target" on each line. node's IDs should be integer, preferably from 0 to n-1.
- neto.txt is the output directed graph: (with hubs removed).

## Modification:

The code can be modified to obtain any other wished output, such as each pair of nodes with the similarity larger than a threshold (the output consisting of all pairs of nodes with nonzero similarity might be too large and not so useful).

The code can be modified to compute any similarity between nodes $u$ and $v$ of the form 
$$f(|\Delta(u)|,|\Delta(v)|, |\Delta(u)\cup \Delta(v)|, |\Delta(u)\cap \Delta(v)|).$$ 
Computing Adammic-Adar is also possible: https://it.wikipedia.org/wiki/Coefficiente_Adamic/Adar

## Performance:
Note that the programs are embarrassingly parallel and using k threads will divide the time by k (up to some large k).

On a commodity machine using a single thread and without removing hubs' edges:

### sim.c
- On https://snap.stanford.edu/data/com-Amazon.html (1M edges): 1 second
- On https://snap.stanford.edu/data/com-Youtube.html (3M edges): 1 minute
- On https://snap.stanford.edu/data/com-Orkut.html (117M edges): 45 minutes

### jaccard.c
- On https://snap.stanford.edu/data/com-Orkut.html (117M edges) with a threshold of 0.5: 10 minutes
- On https://snap.stanford.edu/data/com-Orkut.html (117M edges) with a threshold of 0.8: 3 minutes

## Reference:

The program shows that a "smart" brute-force approach is relatively scalable for this problem. The only problem being the RAM: it does not scale if the input graph does not fit in RAM (i.e., if 2 integers for each edge in the graph cannot be stored in RAM).
 
Graph compression à la Boldi-Vigna could be a solution: http://law.di.unimi.it/datasets.php

## Initial contributors:

Maximilien Danisch 
August 2018 
http://bit.ly/danisch 
maximilien.danisch@gmail.com