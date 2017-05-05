#ifndef INTERVAL_H
#define INTERVAL_H

#include <glib.h>
#include <util.h>
#include <solution.h>

typedef struct _interval{
	int a;
	int b;
	GPtrArray *sigma;
} interval;

void interval_init(interval *p, int a, int b, GPtrArray* jobarray, int njobs);
void interval_free(interval *p);
interval *interval_alloc(int a, int b, GPtrArray *jobarray, int njobs);
void intervals_free(void *p);

gint compare_interval(gconstpointer a, gconstpointer b, gpointer data);
void g_print_interval(gpointer data, gpointer user_data);


#endif // INTERVAL_H


