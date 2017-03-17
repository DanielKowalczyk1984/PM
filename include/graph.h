#ifndef GRAPH_H
#define GRAPH_H

#ifdef __cplusplus
extern "C" {
#endif  // GRAPH_H

typedef struct adjNode {
    int *adj;
    int  node;
    int  weight;
    int  residual_weight;
    int  active;
    int  degree;
    int  color;
    int  profit;
} adjNode;

typedef struct adjGraph {
    adjNode *nodelist;
    int *    adjspace;
    int **   adjMatrix;
    int      vcount;
    int      ecount;

    int *perm;
    int *weightorder;
    int *makespan;
    int  makespancolor;
    int  ncolor;
    int  totweight;
    int  flag;
} adjGraph;

/* Build and Free adjGraph and adjNode */
int adjGraph_build(adjGraph *G,
                   int       vcount,
                   int       ecount,
                   const int elist[],
                   const int weights[]);
int adjGraph_buildquick(adjGraph *G, int vcount, int ecount, int *elist);
int adjGraph_copy(adjGraph *Gdst, const adjGraph *Gsrc);
void adjGraph_init(adjGraph *G);
void adjGraph_free(adjGraph *G);
void adjGraph_freequick(adjGraph *G);
int adjGraph_simplify(adjGraph *G);
int adjGraph_simplifyquick(adjGraph *G);
int adjGraph_complement(adjGraph *Gc, const adjGraph *G);
int adjGraph_get_elist(int *ecount, int *elist[], const adjGraph *G);
int read_adjlist(
    char *f, int *pvcount, int *pecount, int **pelist, int **weightlist);

/* Some usefull functions*/
int adjGraph_delete_unweighted(adjGraph *G,
                               int **    new_nweights,
                               const int nweights[]);
int adjGraph_connectedness(const adjGraph *G);
void graph_print(int ecount, const int elist[]);
int adjGraph_reset_schedule(adjGraph *G);
int adjGraph_adjust_schedule(adjGraph *Gdst, const adjGraph *Gsrc);
void reset_color(adjGraph *G);

/* Some Sorting functions on adjGraph*/
void adjGraph_quicksort_perm(adjNode *nodelist,
                             int *    perm,
                             int      vcount,
                             int (*compareFunction)(adjNode *, adjNode *));
void adjGraph_quicksort(adjNode *nodelist,
                        int      vcount,
                        int (*compareFunction)(adjNode *, adjNode *));
void adjGraph_sort_adjlists_by_id(adjGraph *G);
int adjGraph_degree(adjNode *n1, adjNode *n2);
int adjGraph_weight(adjNode *n1, adjNode *n2);
int adjGraph_invdegree(adjNode *n1, adjNode *n2);
int adjGraph_print(int ecount, const int elist[]);

#ifdef __cplusplus
}
#endif  // GRAPH_H

#endif  // GRAPH_H
