#ifndef INTERVAL_H
#define INTERVAL_H

#include <glib.h>
#include "job.h"
// #include "util.h"

typedef struct _interval {
    int        a;
    int        b;
    int        begin;
    int        key;
    GPtrArray* sigma;
} interval;

typedef struct _interval_pair {
    Job* a;
    Job* b;
    int  left;
    int  right;
} interval_pair;

typedef struct _job_interval_pair {
    Job*      j;
    interval* I;
} job_interval_pair;

void      interval_init(interval*  p,
                        int        a,
                        int        b,
                        int        key,
                        GPtrArray* jobarray,
                        int        nb_jobs);
interval* interval_alloc(int        a,
                         int        b,
                         int        key,
                         GPtrArray* jobarray,
                         int        nb_jobs);
gpointer  g_copy_interval(gconstpointer, gpointer);
gpointer  g_copy_interval_pair(gconstpointer src, gpointer data);
void      interval_free(interval* p);
void      g_interval_free(void* p);

gint g_compare_interval_data(gconstpointer a, gconstpointer b, gpointer data);
// gint g_compare_duration(gconstpointer a, gconstpointer b);
void g_print_interval(gpointer data, gpointer user_data);
void print_interval_pair(GPtrArray* ordered_jobs);

void interval_pair_free(void* p);

#endif  // INTERVAL_H
