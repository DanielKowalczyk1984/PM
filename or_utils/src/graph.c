////////////////////////////////////////////////////////////////
//                                                            //
//  graph.c                                                   //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 21/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#include <defs.h>
#include <graph.h>
#include <util.h>

static int permute_nodes(int*  invorder,
                         int   nb_nodes,
                         int   nb_edges,
                         int*  elist,
                         int*  weights,
                         int** pielist,
                         int** pi_weights);

static void adjNode_SWAP(adjNode* n1, adjNode* n2, adjNode* temp);
#define fill(x)                                 \
    {                                           \
        list[x].degree = 0;                     \
        list[x].node = i;                       \
        list[x].color = -1;                     \
        list[x].weight = weights[x];            \
        G->totweight += weights[x];             \
        for (int j = i; j < G->nb_nodes; j++) { \
            G->adjMatrix[i][j] = 0;             \
            G->adjMatrix[j][i] = 0;             \
        }                                       \
        i++;                                    \
    }

#define fill1()                                 \
    {                                           \
        list->degree = 0;                       \
        list->node = i;                         \
        list->color = -1;                       \
        list->weight = *weights;                \
        G->totweight += *weights;               \
        for (int j = i; j < G->nb_nodes; j++) { \
            G->adjMatrix[i][j] = 0;             \
            G->adjMatrix[j][i] = 0;             \
        }                                       \
        i++;                                    \
        list++;                                 \
        weights++;                              \
    }

#define fill2()                                 \
    {                                           \
        list[*elist].degree++;                  \
        list[*(elist + 1)].degree++;            \
        G->adjMatrix[*elist][*(elist + 1)] = 1; \
        G->adjMatrix[*(elist + 1)][*elist] = 1; \
        elist += 2;                             \
    }

#define fill3(x)                                          \
    {                                                     \
        list[elist[2 * x]].degree++;                      \
        list[elist[2 * x + 1]].degree++;                  \
        G->adjMatrix[elist[2 * x]][elist[2 * x + 1]] = 1; \
        G->adjMatrix[elist[2 * x + 1]][elist[2 * x]] = 1; \
    }

static void fill_graph(adjGraph* G, const int* weights) {
    int      i = 0;
    int      nb_nodes = G->nb_nodes;
    adjNode* list = G->nodelist;

    if (nb_nodes & 1) {
        fill1();
    }

    nb_nodes >>= 1;

    if (nb_nodes & 1) {
        fill1();
        fill1();
    }

    nb_nodes >>= 1;

    while (nb_nodes--) {
        fill(0);
        fill(1);
        fill(2);
        fill(3);
        weights += 4;
        list += 4;
    }
}

static void fill_list(adjGraph* G, const int* elist) {
    adjNode* list = G->nodelist;
    int      nb_edges = G->nb_edges;

    if (nb_edges & 1) {
        fill2();
    }

    nb_edges >>= 1;

    if (nb_edges & 1) {
        fill2();
        fill2();
    }

    nb_edges >>= 1;

    while (nb_edges--) {
        fill3(0);
        fill3(1);
        fill3(2);
        fill3(3);
        elist += 8;
    }
}

void reset_color(adjGraph* G) {
    int      nb_nodes = G->nb_nodes;
    adjNode* list = G->nodelist;
    G->nb_colors = 0;

    if (nb_nodes & 1) {
        list->color = -1;
        list++;
    }

    nb_nodes >>= 1;

    if (nb_nodes & 1) {
        list->color = -1;
        list++;
        list->color = -1;
        list++;
    }

    nb_nodes >>= 1;

    while (nb_nodes--) {
        list[0].color = -1;
        list[1].color = -1;
        list[2].color = -1;
        list[3].color = -1;
        list += 4;
    }
}

int adjGraph_build(adjGraph* G,
                   int       nb_nodes,
                   int       nb_edges,
                   const int elist[],
                   const int weights[]) {
    int      val = 0;
    int      i;
    int*     p;
    adjNode* nodelist;
    adjGraph_init(G);
    G->nb_nodes = nb_nodes;
    G->nb_edges = nb_edges;
    G->nodelist = CC_SAFE_MALLOC(G->nb_nodes, adjNode);
    CCcheck_NULL(G->nodelist, "Out of memory for G->nodelist");
    G->perm = CC_SAFE_MALLOC(G->nb_nodes, int);
    CCcheck_NULL(G->nodelist, "out of memory for G->nodelist");
    G->weightorder = CC_SAFE_MALLOC(nb_nodes, int);
    CCcheck_NULL(G->weightorder, "Out of memory for G->weightorder");
    G->adjMatrix = CC_SAFE_MALLOC(G->nb_nodes, int*);

    for (i = 0; i < G->nb_nodes; i++) {
        G->adjMatrix[i] = (int*)NULL;
        G->adjMatrix[i] = CC_SAFE_MALLOC(G->nb_nodes, int);
        G->perm[i] = i;
        G->weightorder[i] = i;
    }

    nodelist = G->nodelist;

    if (G->nb_edges != 0) {
        G->adjspace = CC_SAFE_MALLOC(2 * G->nb_edges, int);
        CCcheck_NULL(G->adjspace, "out of memory for G->adjspace");
    }

    fill_graph(G, weights);
    fill_list(G, elist);
    p = G->adjspace;

    for (i = 0; i < nb_nodes; ++i) {
        nodelist[i].adj = p;
        p += nodelist[i].degree;
        nodelist[i].degree = 0;
    }

    for (i = 0; i < nb_edges; i++) {
        nodelist[elist[2 * i]].adj[nodelist[elist[2 * i]].degree++] =
            elist[2 * i + 1];
        nodelist[elist[2 * i + 1]].adj[nodelist[elist[2 * i + 1]].degree++] =
            elist[2 * i];
    }

    adjGraph_quicksort_perm(nodelist, G->weightorder, G->nb_nodes,
                            (adjGraph_weight));
    return val;
}

int adjGraph_buildquick(adjGraph* G, int nb_nodes, int nb_edges, int* elist) {
    int      val = 0;
    int      i;
    int*     p;
    adjNode* nodelist;
    adjGraph_init(G);
    G->nb_nodes = nb_nodes;
    G->nb_edges = nb_edges;
    G->nodelist = CC_SAFE_MALLOC(G->nb_nodes, adjNode);
    CCcheck_NULL_2(G->nodelist, "Out of memory for G->nodelist");
    nodelist = G->nodelist;

    if (G->nb_edges) {
        G->adjspace = CC_SAFE_MALLOC(2 * G->nb_edges, int);
        CCcheck_NULL(G->adjspace, "out of memory for G->adjspace");
    }

    for (i = 0; i < nb_nodes; ++i) {
        nodelist[i].degree = 0;
    }

    for (i = 0; i < nb_edges; i++) {
        nodelist[elist[2 * i]].degree++;
        nodelist[elist[2 * i + 1]].degree++;
    }

    p = G->adjspace;

    for (i = 0; i < nb_nodes; ++i) {
        nodelist[i].adj = p;
        p += nodelist[i].degree;
        nodelist[i].degree = 0;
    }

    for (i = 0; i < nb_edges; i++) {
        nodelist[elist[2 * i]].adj[nodelist[elist[2 * i]].degree++] =
            elist[2 * i + 1];
        nodelist[elist[2 * i + 1]].adj[nodelist[elist[2 * i + 1]].degree++] =
            elist[2 * i];
    }

CLEAN:
    return val;
}

int adjGraph_copy(adjGraph* graph_dst, const adjGraph* graph_src) {
    int  val = 0;
    int* elist = (int*)NULL;
    int  nb_edges;
    int* tmp_weightlist = CC_SAFE_MALLOC(graph_src->nb_nodes, int);
    CCcheck_NULL_2(tmp_weightlist, "Failed allocate memory to tmp_weightlist");

    for (int i = 0; i < graph_src->nb_nodes; ++i) {
        tmp_weightlist[i] = graph_src->nodelist[i].weight;
    }

    val = adjGraph_get_elist(&nb_edges, &elist, graph_src);
    CCcheck_val_2(val, "adjGraph_get_elist failed");
    val = adjGraph_build(graph_dst, graph_src->nb_nodes, nb_edges, elist,
                         tmp_weightlist);
    CCcheck_val_2(val, "adjGraph_build failed")

        CLEAN : if (val) {
        adjGraph_free(graph_dst);
    }

    if (elist) {
        free(elist);
    }

    if (tmp_weightlist) {
        free(tmp_weightlist);
    }

    return val;
}

int adjGraph_adjust_schedule(adjGraph* graph_dst, const adjGraph* graph_src) {
    int val = 0;
    int i;

    if (graph_dst->nb_nodes != graph_src->nb_nodes ||
        graph_dst->nb_edges != graph_src->nb_edges) {
        fprintf(stderr, "We are considering here two different graphs\n");
        val = 1;
        goto CLEAN;
    }

    CC_IFFREE(graph_dst->makespan, int);
    graph_dst->nb_colors = graph_src->nb_colors;
    graph_dst->makespan = CC_SAFE_MALLOC(graph_dst->nb_colors, int);
    CCcheck_NULL_2(graph_dst->makespan,
                   "Out of memory for graph_dst->makespan");

    for (i = 0; i < graph_dst->nb_colors; i++) {
        graph_dst->makespan[i] = graph_src->makespan[i];
    }

    for (i = 0; i < graph_dst->nb_nodes; i++) {
        graph_dst->nodelist[i].color = graph_src->nodelist[i].color;
    }

    graph_dst->makespancolor = graph_src->makespancolor;
CLEAN:
    return val;
}

int adjGraph_reset_schedule(adjGraph* G) {
    int val = 0;
    CC_IFFREE(G->makespan, int);
    G->nb_colors = 0;
    reset_color(G);
    G->makespancolor = 0;
    return val;
}

static int comp_node_ids(const void* v1, const void* v2) {
    int id1 = *(const int*)v1;
    int id2 = *(const int*)v2;
    return id1 - id2;
}

static void swap_nodes(int* v1, int* v2) {
    int tmp = *v1;
    *v1 = *v2;
    *v2 = tmp;
}

static int unify_adjlist(int* adjlist, int degree, int* tmp_adjlist) {
    int new_degree = 0;

    if (degree) {
        tmp_adjlist[0] = adjlist[0];
        new_degree++;

        for (int j = 1; j < degree; ++j) {
            if (adjlist[j] != adjlist[j - 1]) {
                tmp_adjlist[new_degree++] = adjlist[j];
            }
        }

        for (int j = 0; j < new_degree; ++j) {
            adjlist[j] = tmp_adjlist[j];
        }
    }

    return new_degree;
}

int adjGraph_simplify(adjGraph* G) {
    int  val = 0;
    int  i, j;
    int  nb_nodes, nb_edges;
    int* tmp_adjlist = CC_SAFE_MALLOC(G->nb_edges, int);
    int* tmp_weightlist = CC_SAFE_MALLOC(G->nb_nodes, int);
    CCcheck_NULL_2(tmp_adjlist, "Failed to allocate memory to tmp_adjlist");
    CCcheck_NULL_2(tmp_weightlist,
                   "Failed to allocate memory to tmp_weightlist");

    for (i = 0; i < G->nb_nodes; i++) {
        tmp_weightlist[i] = G->nodelist[i].weight;
    }

    for (i = 0; i < G->nb_nodes; ++i) {
        int new_degree;
        int nloops = 0;
        qsort(G->nodelist[i].adj, G->nodelist[i].degree, sizeof(int),
              comp_node_ids);
        new_degree = unify_adjlist(G->nodelist[i].adj, G->nodelist[i].degree,
                                   tmp_adjlist);

        if (dbg_lvl() > 1 && new_degree != G->nodelist[i].degree) {
            printf("Removed %d edges from node %d\n",
                   G->nodelist[i].degree - new_degree, i);
        }

        G->nodelist[i].degree = new_degree;

        for (j = 0; j < G->nodelist[i].degree; ++j) {
            if (G->nodelist[i].adj[j] == i) {
                nloops++;
                swap_nodes(&(G->nodelist[i].adj[j]),
                           &(G->nodelist[i].adj[G->nodelist[i].degree - 1]));
                --G->nodelist[i].degree;
                --j;
            }
        }

        if (dbg_lvl() > 1 && nloops) {
            printf("Removed %d loop(s) from  node %d\n", nloops, i);
        }
    }

    if (tmp_adjlist) {
        free(tmp_adjlist);
        tmp_adjlist = (int*)NULL;
    }

    nb_nodes = G->nb_nodes;
    val = adjGraph_get_elist(&nb_edges, &tmp_adjlist, G);
    CCcheck_val(val, "Failed in adjGraph_get_elist");
    adjGraph_freequick(G);
    val = adjGraph_build(G, nb_nodes, nb_edges, tmp_adjlist, tmp_weightlist);
    CCcheck_val(val, "Failed adjGraph_build");

CLEAN:
    if (tmp_adjlist) {
        free(tmp_adjlist);
    }

    if (tmp_weightlist) {
        free(tmp_weightlist);
    }

    return val;
}

int adjGraph_simplifyquick(adjGraph* G) {
    int  val = 0;
    int  i, j;
    int* tmp_adjlist = (int*)NULL;
    int  nb_nodes, nb_edges;

    if (G->nb_edges != 0) {
        tmp_adjlist = CC_SAFE_MALLOC(G->nb_edges, int);
        CCcheck_NULL(tmp_adjlist, "Failed to allocate memory to tmp_adjlist");
    }

    for (i = 0; i < G->nb_nodes; ++i) {
        int new_degree;
        int nloops = 0;
        qsort(G->nodelist[i].adj, G->nodelist[i].degree, sizeof(int),
              comp_node_ids);
        new_degree = unify_adjlist(G->nodelist[i].adj, G->nodelist[i].degree,
                                   tmp_adjlist);

        if (dbg_lvl() > 1 && new_degree != G->nodelist[i].degree) {
            printf("Removed %d edges from node %d\n",
                   G->nodelist[i].degree - new_degree, i);
        }

        G->nodelist[i].degree = new_degree;

        for (j = 0; j < G->nodelist[i].degree; ++j) {
            if (G->nodelist[i].adj[j] == i) {
                nloops++;
                swap_nodes(&(G->nodelist[i].adj[j]),
                           &(G->nodelist[i].adj[G->nodelist[i].degree - 1]));
                --G->nodelist[i].degree;
                --j;
            }
        }

        if (dbg_lvl() > 1 && nloops) {
            printf("Removed %d loop(s) from  node %d\n", nloops, i);
        }
    }

    if (tmp_adjlist) {
        free(tmp_adjlist);
        tmp_adjlist = (int*)NULL;
    }

    nb_nodes = G->nb_nodes;
    val = adjGraph_get_elist(&nb_edges, &tmp_adjlist, G);
    CCcheck_val(val, "Failed in adjGraph_get_elist");
    adjGraph_freequick(G);
    val = adjGraph_buildquick(G, nb_nodes, nb_edges, tmp_adjlist);
    CCcheck_val(val, "Failed adjGraph_build");

    if (tmp_adjlist) {
        free(tmp_adjlist);
    }

    return val;
}

int adjGraph_complement(adjGraph* Gc, const adjGraph* G) {
    int  val = 0;
    int  nb_edges = 0;
    int  nb_edges_chk = 0;
    int  v_i, a_i, na;
    int* elist = (int*)NULL;
    int* tmp_weightlist = CC_SAFE_MALLOC(G->nb_nodes, int);
    CCcheck_NULL_2(tmp_weightlist, "No memory for tmp_weightlist");

    for (int i = 0; i < G->nb_nodes; ++i) {
        tmp_weightlist[i] = G->nodelist[i].weight;
    }

    val = adjGraph_copy(Gc, G);
    CCcheck_val_2(val, "Failed adjGraph_copy");
    val = adjGraph_simplify(Gc);
    CCcheck_val_2(val, "Failed in adjGraph_simplify");
    nb_edges_chk = (Gc->nb_nodes * (Gc->nb_nodes - 1)) / 2 - Gc->nb_edges;

    if (nb_edges_chk) {
        elist = CC_SAFE_MALLOC(2 * nb_edges_chk, int);
        CCcheck_NULL_2(elist, "Failed to allocate memory to elist");

        for (v_i = 0; v_i < Gc->nb_nodes; ++v_i) {
            adjNode* v = &(Gc->nodelist[v_i]);
            int      a = -1;
            a_i = 0;

            for (na = v_i + 1; na < Gc->nb_nodes; ++na) {
                while (a_i < v->degree && a < na) {
                    a = v->adj[a_i];
                    ++a_i;
                }

                if (na != a) {
                    elist[2 * nb_edges] = v_i;
                    elist[2 * nb_edges + 1] = na;
                    ++nb_edges;
                }
            }
        }

        assert(nb_edges == nb_edges_chk);
    }

    adjGraph_free(Gc);
    val = adjGraph_build(Gc, G->nb_nodes, nb_edges_chk, elist, tmp_weightlist);
    CCcheck_val(val, "Failed adjGraph_build");
CLEAN:

    if (val) {
        adjGraph_free(Gc);
    }

    if (elist) {
        free(elist);
    }

    if (tmp_weightlist) {
        free(tmp_weightlist);
    }

    return val;
}

void adjGraph_init(adjGraph* G) {
    if (G) {
        G->nodelist = (adjNode*)NULL;
        G->adjspace = (int*)NULL;
        G->makespan = (int*)NULL;
        G->perm = (int*)NULL;
        G->weightorder = (int*)NULL;
        G->adjMatrix = (int**)NULL;
        G->nb_colors = 0;
        G->nb_nodes = 0;
        G->nb_edges = 0;
        G->totweight = 0;
        G->makespancolor = 0;
        G->flag = 0;
    }
}

void adjGraph_free(adjGraph* G) {
    if (G) {
        for (int i = 0; i < G->nb_nodes; i++) {
            CC_IFFREE(G->adjMatrix[i], int);
        }

        CC_IFFREE(G->adjMatrix, int*)
        CC_IFFREE(G->nodelist, adjNode);
        CC_IFFREE(G->adjspace, int);
        CC_IFFREE(G->makespan, int);
        CC_IFFREE(G->perm, int);
        CC_IFFREE(G->weightorder, int);
        adjGraph_init(G);
    }
}

void adjGraph_freequick(adjGraph* G) {
    if (G) {
        CC_IFFREE(G->nodelist, adjNode);
        CC_IFFREE(G->adjspace, int);
        adjGraph_init(G);
    }
}

int read_adjlist(char* f,
                 int*  pvcount,
                 int*  pecount,
                 int** pelist,
                 int** weightlist) {
    int         val = 0;
    int         nb_nodes = 0, nb_edges = 0, count = 0, prob = 0;
    int         curnode, curweight, nodeadj;
    int*        elist = (int*)NULL;
    int*        weight = (int*)NULL;
    char        buf[256], *p;
    int         bufsize;
    const char* delim = " \n";
    char*       data = (char*)NULL;
    char*       buf2 = (char*)NULL;
    int*        perm = (int*)NULL;
    int*        inv_perm = (int*)NULL;
    FILE*       in = fopen(f, "r");

    if (in != (FILE*)NULL) {
        if (fgets(buf, 254, in) != NULL) {
            p = buf;

            if (p[0] == 'p') {
                if (prob) {
                    fprintf(stderr, "ERROR: in this file we have to p lines\n");
                    val = 1;
                    goto CLEAN;
                }

                prob = 1;
                strtok(p, delim);
                data = strtok(NULL, delim);
                sscanf(data, "%d", &nb_nodes);
                data = strtok(NULL, delim);
                sscanf(data, "%d", &nb_edges);
                bufsize =
                    2 * nb_nodes * (2 + (int)ceil(log((double)nb_nodes + 10)));
                buf2 = (char*)CC_SAFE_MALLOC(bufsize, char);
                CCcheck_NULL_2(buf2, "Failed to allocate buf2");

                if (nb_edges != 0) {
                    elist = CC_SAFE_MALLOC(2 * nb_edges, int);
                    CCcheck_NULL_2(elist, "out of memory for elist");
                }

                weight = CC_SAFE_MALLOC(nb_nodes, int);
                CCcheck_NULL_2(weight, "out of memory for weight");
            } else {
                fprintf(
                    stderr,
                    "File has to give first the number vertices and edges.\n");
                val = 1;
                goto CLEAN;
            }
        } else {
            val = 1;
            goto CLEAN;
        }

        while (fgets(buf2, bufsize, in) != (char*)NULL) {
            p = buf2;

            if (p[0] == 'p') {
                if (prob) {
                    fprintf(stderr, "ERROR: in this file we have to p lines\n");
                    val = 1;
                    goto CLEAN;
                }
            } else if (p[0] == 'n') {
                if (!prob) {
                    fprintf(stderr, "ERROR n before p in file\n");
                    val = 1;
                    goto CLEAN;
                }

                if (count > nb_edges) {
                    fprintf(stderr, "ERROR: too many edges in file\n");
                    val = 1;
                    goto CLEAN;
                }

                strtok(p, delim);
                data = strtok(NULL, delim);
                sscanf(data, "%d", &curweight);
                data = strtok(NULL, delim);
                sscanf(data, "%d", &curnode);
                weight[curnode] = curweight;
                data = strtok(NULL, delim);

                while (data != NULL) {
                    elist[2 * count] = curnode;
                    sscanf(data, "%d", &nodeadj);
                    elist[2 * count + 1] = nodeadj;
                    data = strtok(NULL, delim);
                    count++;
                }
            }
        }

        perm = CC_SAFE_MALLOC(nb_nodes, int);
        CCcheck_NULL_2(perm, "Failed to allocate memory");

        for (int i = 0; i < nb_nodes; i++) {
            perm[i] = i;
        }

        CCutil_int_perm_quicksort_0(perm, weight, nb_nodes);
        inv_perm = CC_SAFE_MALLOC(nb_nodes, int);
        CCcheck_NULL_2(inv_perm, "Failed to allocate memory");

        for (int i = 0; i < nb_nodes; ++i) {
            inv_perm[perm[i]] = i;
        }

        permute_nodes(inv_perm, nb_nodes, nb_edges, elist, weight, pelist,
                      weightlist);
        *pvcount = nb_nodes;
        *pecount = count;
    } else {
        fprintf(stderr, "Unable to open file %s\n", f);
        val = 1;
        goto CLEAN;
    }

CLEAN:

    if (val) {
        CC_IFFREE(*pelist, int);
        CC_IFFREE(*weightlist, int);
    }

    CC_IFFREE(elist, int);
    CC_IFFREE(weight, int);
    CC_IFFREE(buf2, char);
    CC_IFFREE(perm, int);
    CC_IFFREE(inv_perm, int);

    if (in) {
        fclose(in);
    }

    return val;
}

void graph_print(int nb_edges, const int elist[]) {
    for (int i = 0; i < nb_edges; ++i) {
        printf("e %d %d\n", elist[2 * i], elist[2 * i + 1]);
    }
}

int adjGraph_get_elist(int* nb_edges, int* elist[], const adjGraph* G) {
    int val = 0;
    int i;
    CC_IFFREE(*elist, int);
    *nb_edges = 0;

    for (i = 0; i < G->nb_nodes; i++) {
        *nb_edges += G->nodelist[i].degree;
    }

    assert(*nb_edges % 2 == 0);

    if (*nb_edges != 0) {
        /* code */
        (*elist) = CC_SAFE_MALLOC((*nb_edges), int);
        CCcheck_NULL_2(*elist, "Out of memory for elist");
        *nb_edges = 0;

        for (i = 0; i < G->nb_nodes; i++) {
            int j;

            for (j = 0; j < G->nodelist[i].degree; ++j) {
                if (G->nodelist[i].adj[j] > i) {
                    (*elist)[(*nb_edges) * 2] = i;
                    (*elist)[(*nb_edges) * 2 + 1] = G->nodelist[i].adj[j];
                    (*nb_edges)++;
                }
            }
        }
    }

CLEAN:
    return val;
}

int adjGraph_delete_unweighted(adjGraph* G,
                               int**     new_node_weights,
                               const int node_weights[]) {
    int  val = 0;
    int  i, a_i;
    int  nb_nodes = 0;
    int  nb_edges = 0;
    int* nodes_map = CC_SAFE_MALLOC(G->nb_nodes, int);
    int* newelist = CC_SAFE_MALLOC(2 * G->nb_edges, int);
    int* new_weightlist = CC_SAFE_MALLOC(G->nb_nodes, int);

    CCcheck_NULL_2(nodes_map, "Failed to allocate nodes_map");
    CCcheck_NULL_2(nodes_map, "Failed to allocate newelist");
    CCcheck_NULL_2(new_weightlist,
                   "Failed to allocate memory to new_weightlist");

    for (i = 0; i < G->nb_nodes; ++i) {
        if (node_weights[i] == 0) {
            nodes_map[i] = -1;
            G->nodelist[i].degree = 0;
        } else {
            int a = 0;
            nodes_map[i] = nb_nodes++;

            for (a_i = 0; a_i < G->nodelist[i].degree && a < i; ++a_i) {
                a = G->nodelist[i].adj[a_i];

                if (a < i && nodes_map[a] != -1) {
                    newelist[2 * nb_edges] = nodes_map[a];
                    newelist[2 * nb_edges + 1] = nodes_map[i];
                    ++nb_edges;
                }
            }
        }
    }

    *new_node_weights = CC_SAFE_MALLOC(nb_nodes, int);
    CCcheck_NULL_2(*new_node_weights, "Failed to allocate nne");

    for (i = 0; i < G->nb_nodes; ++i) {
        int ni = nodes_map[i];

        if (ni != -1) {
            (*new_node_weights)[ni] = node_weights[i];
            new_weightlist[ni] = G->nodelist[i].weight;
        }
    }

    adjGraph_free(G);
    val = adjGraph_build(G, nb_nodes, nb_edges, newelist, new_weightlist);
    CCcheck_val_2(val, "Failed in COLORadjgraph_build");

    if (dbg_lvl() > 1) {
        printf("Reduced graph has %d nodes and %d edges.\n", nb_nodes,
               nb_edges);
    }

CLEAN:

    if (nodes_map) {
        free(nodes_map);
    }

    if (newelist) {
        free(newelist);
    }

    if (new_weightlist) {
        free(new_weightlist);
    }

    if (val) {
        if (*new_node_weights) {
            free(*new_node_weights);
        }

        *new_node_weights = (int*)NULL;
    }

    return val;
}

static void adjNode_SWAP(adjNode* n1, adjNode* n2, adjNode* temp) {
    if (n1 != n2) {
        memcpy(temp, n2, sizeof(adjNode));
        memcpy(n2, n1, sizeof(adjNode));
        memcpy(n1, temp, sizeof(adjNode));
    }
}

void adjGraph_quicksort(adjNode* nodelist,
                        int      nb_nodes,
                        int (*compareFunction)(adjNode*, adjNode*)) {
    int     i, j;
    adjNode t, temp;

    if (nb_nodes <= 1) {
        return;
    }

    adjNode_SWAP(&nodelist[0], &nodelist[(nb_nodes - 1) / 2], &temp);
    i = 0;
    j = nb_nodes;
    t = nodelist[0];

    while (1) {
        do {
            i++;
        } while (i < nb_nodes && (*compareFunction)(&nodelist[i], &t));

        do {
            j--;
        } while ((*compareFunction)(&t, &nodelist[j]));

        if (j < i) {
            break;
        }

        adjNode_SWAP(&nodelist[i], &nodelist[j], &temp);
    }

    adjNode_SWAP(&nodelist[0], &nodelist[j], &temp);
    adjGraph_quicksort(nodelist, j, (*compareFunction));
    adjGraph_quicksort(nodelist + i, nb_nodes - i, (*compareFunction));
    return;
}

void adjGraph_sort_adjlists_by_id(adjGraph* G) {
    int i;

    for (i = 0; i < G->nb_nodes; ++i) {
        qsort(G->nodelist[i].adj, G->nodelist[i].degree, sizeof(int),
              comp_node_ids);
    }
}

void adjGraph_quicksort_perm(adjNode* nodelist,
                             int*     perm,
                             int      nb_nodes,
                             int (*compareFunction)(adjNode*, adjNode*)) {
    int     i, j;
    int     temp;
    adjNode t;

    if (nb_nodes <= 1) {
        return;
    }

    CC_SWAP(perm[0], perm[(nb_nodes - 1) / 2], temp);
    i = 0;
    j = nb_nodes;
    t = nodelist[perm[0]];

    while (1) {
        do {
            i++;
        } while (i < nb_nodes && (*compareFunction)(&nodelist[perm[i]], &t));

        do {
            j--;
        } while ((*compareFunction)(&t, &nodelist[perm[j]]));

        if (j < i) {
            break;
        }

        CC_SWAP(perm[i], perm[j], temp);
    }

    CC_SWAP(perm[0], perm[j], temp);
    adjGraph_quicksort_perm(nodelist, perm, j, (*compareFunction));
    adjGraph_quicksort_perm(nodelist, perm + i, nb_nodes - i,
                            (*compareFunction));
    return;
}

int adjGraph_degree(adjNode* n1, adjNode* n2) {
    if (n1->degree != n2->degree) {
        return n1->degree > n2->degree;
    }

    return 0;
}

int adjGraph_weight(adjNode* n1, adjNode* n2) {
    if (n1->weight != n2->weight) {
        return n1->weight > n2->weight;
    }

    return 0;
}

int adjGraph_invdegree(adjNode* n1, adjNode* n2) {
    if (n1->degree != n2->degree) {
        return n1->degree < n2->degree;
    }

    return 0;
}

int adjGraph_print(int nb_edges, const int elist[]) {
    int i;

    for (i = 0; i < nb_edges; ++i) {
        printf("e %d %d\n", elist[2 * i], elist[2 * i + 1]);
    }

    return 0;
}

static int permute_nodes(int*  invorder,
                         int   nb_nodes,
                         int   nb_edges,
                         int*  elist,
                         int*  weights,
                         int** pielist,
                         int** pi_weights) {
    int  i, val = 0;
    int *ielist = (int*)NULL, *inv_weights = (int*)NULL;
    *pielist = (int*)NULL;
    *pi_weights = (int*)NULL;

    if (nb_edges != 0) {
        ielist = CC_SAFE_MALLOC(2 * nb_edges, int);
        CCcheck_NULL_2(pielist, "out of memory for pielist");
    }

    inv_weights = CC_SAFE_MALLOC(nb_nodes, int);
    CCcheck_NULL_2(inv_weights, "out of memory for inv_weights");

    for (i = 0; i < nb_edges; i++) {
        if (invorder[elist[2 * i]] < invorder[elist[2 * i + 1]]) {
            ielist[2 * i] = invorder[elist[2 * i]];
            ielist[2 * i + 1] = invorder[elist[2 * i + 1]];
        } else {
            ielist[2 * i] = invorder[elist[2 * i + 1]];
            ielist[2 * i + 1] = invorder[elist[2 * i]];
        }
    }

    for (i = 0; i < nb_nodes; i++) {
        inv_weights[invorder[i]] = weights[i];
    }

    *pielist = ielist;
    *pi_weights = inv_weights;
CLEAN:

    if (val) {
        CC_IFFREE(ielist, int);
        CC_IFFREE(inv_weights, int);
    }

    return val;
}
