#include <interval.h>

void g_print_interval(gpointer data, gpointer user_data)
{
    interval *a = (interval *)data;
    printf("%d %d: ", a->a, a->b);
    g_ptr_array_foreach(a->sigma, g_print_job, NULL);
    printf("\n");
}

gint compare_interval(gconstpointer a, gconstpointer b, gpointer data)
{
    const Job *x = *(Job * const *)a;
    const Job *y = *(Job * const *)b;
    interval *user_data = (interval *)data;
    int        diff = user_data->b - user_data->a;
    double w_x = (x->duetime <= user_data->a)
                 ? (double) x->weight / x->processingime
                 : 0.0;
    double w_y = (y->duetime <= user_data->a)
                 ? (double) y->weight / y->processingime
                 : 0.0;

    if (x->processingime > diff) {
        if (y->processingime <= diff) {
            return -1;
        } else {
            if (w_x > w_y) {
                return -1;
            } else if (w_y > w_x) {
                return 1;
            } else if (x->processingime > y->processingime) {
                return -1;
            } else if (y->processingime > x->processingime) {
                return 1;
            }

            return x->job - y->job;
        }
    } else {
        if (y->processingime > diff) {
            return 1;
        } else {
            if (w_x > w_y) {
                return -1;
            } else if (w_y > w_x) {
                return 1;
            } else if (x->processingime > y->processingime) {
                return -1;
            } else if (y->processingime > x->processingime) {
                return 1;
            }

            return x->job - y->job;
        }
    }
}

void interval_init(interval *p, int a, int b, GPtrArray *jobarray, int njobs)
{
    p->a = a;
    p->b = b;
    p->sigma = g_ptr_array_new();
    p->begin = 0;
    Job *j;

    for (int i = 0; i < njobs; ++i) {
        g_ptr_array_add(p->sigma, g_ptr_array_index(jobarray, i));
    }

    g_ptr_array_sort_with_data(p->sigma, compare_interval, p);

    j = (Job *) g_ptr_array_index(p->sigma, p->begin);
    while(p->b - p->a <= j->processingime && p->begin < (int)jobarray->len) {
        j = (Job *) g_ptr_array_index(p->sigma,p->begin);
        p->begin++;
    }
}

interval *interval_alloc(int a, int b, GPtrArray *jobarray, int njobs)
{
    interval *p = CC_SAFE_MALLOC(1, interval);
    CCcheck_NULL_3(p, "Failed to allocate memory")
    interval_init(p, a, b, jobarray, njobs);
CLEAN:
    return p;
}

void interval_free(interval *p)
{
    if (p != (interval *)NULL) {
        g_ptr_array_free(p->sigma, TRUE);
    }
}

void g_interval_free(void *p)
{
    if (p != NULL) {
        interval *part = (interval *)p;
        g_ptr_array_free(part->sigma, TRUE);
        CC_FREE(part, interval);
    }
}

void interval_pair_free(void *p){
    if(p != NULL) {
        interval_pair *tmp = (interval_pair *) p;
        CC_FREE(tmp, interval_pair);
    }
}
