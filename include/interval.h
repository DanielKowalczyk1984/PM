#ifndef INTERVAL_H
#define INTERVAL_H

#include <glib.h>
#include <util.h>
#include <solution.h>

typedef struct _interval{
	int a;
	int b;
	int begin;
    int key;
	GPtrArray *sigma;
} interval;

typedef struct _interval_pair {
	Job *a;
	Job *b;
	int left;
	int right;
} interval_pair;

typedef struct _job_interval_pair {
    Job *j;
    interval *I;
    int take;
} job_interval_pair;

void interval_init(interval *p, int a, int b,int key, GPtrArray* jobarray, int njobs);
interval *interval_alloc(int a, int b, int key, GPtrArray *jobarray, int njobs);
interval *interval_copy(interval *src);
void interval_free(interval *p);
void g_interval_free(void *p);

gint compare_interval(gconstpointer a, gconstpointer b, gpointer data);
void g_print_interval(gpointer data, gpointer user_data);
void print_interval_pair(GPtrArray *ordered_jobs);
void count_jobs_interval_pair(GPtrArray *ordered_jobs);

void interval_pair_free(void *p);


#endif // INTERVAL_H


