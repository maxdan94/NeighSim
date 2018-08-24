/*
gcc cosine.c -O9 -o cosine -lm -fopenmp
./cosine n_threads a_threshold net
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <omp.h>


#define NLINKS 500000000 //Maximum number of links, will automatically increase if needed.

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	//edge list structure:
	unsigned n;//number of nodes
	unsigned e;//number of edges
	edge *edges;//list of edges

	//relabel in degree ordering order
	unsigned *rank;
	unsigned *map;

	//neighborhoods:
	unsigned *d0; //original degrees
	unsigned *d; //degrees
	unsigned *cd; //cumulative degrees: (start with 0) length=dim+1
	unsigned *adj; //list of neighbors
} graph;


//compute the maximum of three unsigned
unsigned max3(unsigned a,unsigned b,unsigned c){
	a=(a>b) ? a : b;
	return (a>c) ? a : c;
}

//used in qsort
int cmpfunc(void const *a, void const *b){
	unsigned const *pa = a;
	unsigned const *pb = b;
	if (*pa<*pb)
		return -1;
	return 1;
}

//reading the edgelist from file
graph* readedgelist(char* edgelist){
	unsigned e1=NLINKS;
	graph *g=malloc(sizeof(graph));
	FILE *file;

	g->n=0;
	g->e=0;
	file=fopen(edgelist,"r");
	g->edges=malloc(e1*sizeof(edge));
	while (fscanf(file,"%u %u\n", &(g->edges[g->e].s), &(g->edges[g->e].t))==2) {
		g->n=max3(g->n,g->edges[g->e].s,g->edges[g->e].t);
		if (g->e++==e1) {
			e1+=NLINKS;
			g->edges=realloc(g->edges,e1*sizeof(edge));
		}
	}
	fclose(file);
	g->n++;
	g->edges=realloc(g->edges,g->e*sizeof(edge));

	return g;
}

typedef struct {
	unsigned node;
	unsigned deg;
} nodedeg;

int compare_nodedeg(void const *a, void const *b){
	nodedeg const *pa = a;
	nodedeg const *pb = b;
	return (pa->deg < pb->deg) ? -1 : 1;
}

void degord(graph *g, unsigned dmax) {
	unsigned i,j;
	unsigned *d=calloc(g->n,sizeof(unsigned));
	nodedeg *nodedeglist=malloc(g->n*sizeof(nodedeg));
	for (i=0;i<g->e;i++) {
		d[g->edges[i].s]++;
		d[g->edges[i].t]++;
	}
	for (i=0;i<g->n;i++) {
		nodedeglist[i].node=i;
		nodedeglist[i].deg=0;
	}
	for (i=0;i<g->e;i++) {
		if (d[g->edges[i].t]<=dmax){
			nodedeglist[g->edges[i].s].deg++;
		}
		if (d[g->edges[i].s]<=dmax){
			nodedeglist[g->edges[i].t].deg++;
		}
	}
	qsort(nodedeglist,g->n,sizeof(nodedeg),compare_nodedeg);
	//printf("%u %u %u\n",nodedeglist[0].deg,nodedeglist[100].deg,nodedeglist[g->n-1].deg);
	g->rank=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++) {
			g->rank[nodedeglist[i].node]=i;
	}
	free(nodedeglist);
}

void relabel(graph *g) {
	unsigned i,j,source,target;
	g->map=malloc(g->n*sizeof(unsigned));
	for (i=0;i<g->n;i++) {
		g->map[g->rank[i]]=i;
	}
	for (i=0;i<g->e;i++) {
		g->edges[i].s=g->rank[g->edges[i].s];
		g->edges[i].t=g->rank[g->edges[i].t];
	}
}



//Building the special graph structure
void mkgraph(graph *g,unsigned dmax){
	unsigned i,s,t,max,tmp;

	g->d0=calloc(g->n,sizeof(unsigned));
	g->d=calloc(g->n,sizeof(unsigned));
	g->cd=malloc((g->n+1)*sizeof(unsigned));
	g->adj=malloc(2*g->e*sizeof(unsigned));
	for (i=0;i<g->e;i++) {
		g->d0[g->edges[i].s]++;
		g->d0[g->edges[i].t]++;
	}
	max=0;
	for (i=0;i<g->n;i++) {
		max=(g->d0[i]>max)?g->d0[i]:max;
	}
	printf("Maximum degree: %u\n",max);
	for (i=0;i<g->e;i++) {
		if (g->d0[g->edges[i].t]<=dmax){
			g->d[g->edges[i].s]++;
		}
		if (g->d0[g->edges[i].s]<=dmax){
			g->d[g->edges[i].t]++;
		}
	}
	g->cd[0]=0;
	for (i=1;i<g->n+1;i++) {
		g->cd[i]=g->cd[i-1]+g->d0[i-1];
	}
	bzero(g->d0,(g->n)*sizeof(unsigned));
	for (i=0;i<g->e;i++) {
		s=g->edges[i].s;
		t=g->edges[i].t;
		g->adj[g->cd[s] + g->d0[s]++ ]=t;
		g->adj[g->cd[t] + g->d0[t]++ ]=s;
	}

	for (i=0;i<g->n;i++) {
		qsort(g->adj+g->cd[i],g->d0[i],sizeof(unsigned),cmpfunc);
	}

}


void freegraph(graph *g){
	free(g->edges);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g);
}

unsigned binsearch(unsigned *tab, unsigned l, unsigned r, unsigned x){
	//printf("%u %u\n",l,r);
	unsigned i;

	unsigned mid = l + (r - l + 1)/2;
	if (tab[mid] == x)//x always in tab
		return mid+1;
	if (tab[mid] > x)
		return binsearch(tab, l, mid-1, x); 
	return binsearch(tab, mid+1, r, x);

}


//histogram of cosine values
unsigned long long* cosine(graph *g,double a,unsigned dmax){
	unsigned i,j,k,u,v,w,n;
	double val;
	unsigned long long *hist_p,*hist=calloc(10,sizeof(unsigned long long));
	bool *tab;
	unsigned *list,*inter;
	#pragma omp parallel private(i,j,k,u,v,w,val,tab,hist_p,inter,list,n)
	{
	hist_p=calloc(10,sizeof(unsigned long long));
	tab=calloc(g->n,sizeof(bool));
	list=malloc(g->n*sizeof(unsigned));
	inter=calloc(g->n,sizeof(unsigned));

	#pragma omp for schedule(dynamic, 1) nowait
	for (u=0;u<g->n;u++){//embarrassingly parallel...
		n=0;
		for (i=g->cd[u];i<g->cd[u+1];i++){
			v=g->adj[i];
			if (g->d0[v]==0)///////
				continue;
			if (g->d0[v]>dmax)
				continue;
			for (j=binsearch(g->adj,g->cd[v],g->cd[v+1]-1,u);j<g->cd[v+1];j++){
				w=g->adj[j];
				if (((double)g->d[u])/((double)(g->d[w]))<a){
					break;
				}
				if(tab[w]==0){
					list[n++]=w;
					tab[w]=1;
				}
				inter[w]++;
			}
			//printf("outb\n");
		}
		for (i=0;i<n;i++){
			w=list[i];
			val=((double)inter[w])/((double)(g->d[u]+g->d[w]-inter[w]));
			//printf("%u %u %le\n",u,w,val);//to print the pairs and similarity
			if (val>0.9){
				hist_p[9]++;
			}
			else {
				hist_p[(int)(floor(val*10))]++;
			}
			tab[w]=0;
			inter[w]=0;
		}
	}
	free(tab);
	free(list);
	free(inter);
	#pragma omp critical
	{
		for (i=0;i<10;i++){
			hist[i]+=hist_p[i];
		}
		free(hist_p);
	}
	}
	return hist;
}


int main(int argc,char** argv){
	graph* g;
	unsigned i,dmax;
	unsigned long long tot=0;
	unsigned long long *hist;
	double a;
	time_t t0,t1,t2;
	t1=time(NULL);
	t0=t1;

	omp_set_num_threads(atoi(argv[1]));
	printf("Similarities greater than: %s\n",argv[2]);
	a=atof(argv[2]);
	printf("Only taking into account common neighbors with degree > %s\n",argv[3]);
	dmax=atoi(argv[3]);
	printf("Reading edgelist from file %s\n",argv[4]);
	g=readedgelist(argv[4]);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Number of nodes: %u\n",g->n);
	printf("Number of edges: %u\n",g->e);

	printf("Degree Ordering\n");
	degord(g,dmax);
	relabel(g);

	printf("Building Graph\n");

	mkgraph(g,dmax);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	printf("Computing jaccard similarities\n");

	hist=cosine(g,a,dmax);

	t2=time(NULL);
	printf("- Time = %ldh%ldm%lds\n",(t2-t1)/3600,((t2-t1)%3600)/60,((t2-t1)%60));
	t1=t2;

	freegraph(g);

	printf("- Overall time = %ldh%ldm%lds\n",(t2-t0)/3600,((t2-t0)%3600)/60,((t2-t0)%60));

	printf("ONLY JACCARD SIMILARITIES GREATER THAN %lf ARE CORRECT\n",ceil(a*10)/10);
	printf("Number of jaccard similarities in\n");
	for (i=0;i<9;i++){
		printf("]0.%u, 0.%u] = %llu\n",i,i+1,hist[i]);
		tot+=hist[i];
	}
	printf("]0.9, 1.0] = %llu\n",hist[9]);
	tot+=hist[9];
	printf("Number of non-zero jaccard similarities = %llu\n",tot);

	return 0;
}


