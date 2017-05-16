////////////////////////////////////////////////////////////////
//                                                            //
//  util.h                                                    //
//  PMC                                                       //
//                                                            //
//  Created by Daniel on 20/02/14.                            //
//  Copyright (c) 2014 Daniel Kowalczyk. All rights reserved. //
//                                                            //
////////////////////////////////////////////////////////////////

#ifndef __UTIL_H
#define __UTIL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <defs.h>
#include <sys/utsname.h>

#define CCutil_MAXDBL (1e30)
#define CCutil_MAXINT (2147483647)

#define MAX_PNAME_LEN 128

#define CC_SWAP(a, b, t) (((t) = (a)), ((a) = (b)), ((b) = (t)))

#define CC_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define CC_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define CC_MIN_PTR(x, y) (((x) < (y)) ? (&x) : (&y))
#define CC_MAX_PTR(x, y) (((x) > (y)) ? (&x) : (&y))
#define CC_OURABS(a) (((a) >= 0) ? (a) : -(a))

#define CCcheck_val_2(val, msg)                                                \
    {                                                                          \
        if ((val)) {                                                           \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
            goto CLEAN;                                                        \
        }                                                                      \
    }

#define CCcheck_val(val, msg)                                                  \
    {                                                                          \
        if ((val)) {                                                           \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
        }                                                                      \
    }

#define CCcheck_NULL(item, msg)                                                \
    {                                                                          \
        if (!(item)) {                                                         \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
        }                                                                      \
    }

#define CCcheck_NULL_3(item, msg)                                              \
    {                                                                          \
        if (!(item)) {                                                         \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
            goto CLEAN;                                                        \
        }                                                                      \
    }

#define CCcheck_NULL_2(item, msg)                                              \
    {                                                                          \
        if (!(item)) {                                                         \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
            val = 1;                                                           \
            goto CLEAN;                                                        \
        }                                                                      \
    }

#define CCcheck_fileio(rval, msg)                                              \
    {                                                                          \
        if (pval < 0) {                                                        \
            fflush(stdout);                                                    \
            fprintf(stderr, "%s at %s, line %d\n", (msg), __FILE__, __LINE__); \
            val = 1;                                                           \
            goto CLEAN;                                                        \
        }                                                                      \
    }

#define CC_SBUFFER_SIZE (4000)
#define CC_SFNAME_SIZE (32)

typedef struct CC_SFILE {
    int status;
    int desc;
    int type;
    int chars_in_buffer;
    int current_buffer_char; /* only used for reading */
    int bits_in_last_char;   /* writing: number of empty bits in
                              * buffer[chars_in_buffer];
                              * reading: number of full bits in
                              * buffer[?] */
    int           pos;
    char          fname[CC_SFNAME_SIZE];
    char          hname[CC_SFNAME_SIZE];
    unsigned char buffer[CC_SBUFFER_SIZE];
} CC_SFILE;

typedef struct CCrandstate {
    int a;
    int b;
    int arr[55];
} CCrandstate;

int  dbg_lvl(void);
void set_dbg_lvl(int dbglvl);

/****************************************************************************/
/*                                                                          */
/*                             allocrus.c                                   */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                   MEMORY ALLOCATION MACROS                               */
/*                                                                          */
/*                           TSP CODE                                       */
/*                                                                          */
/*                                                                          */
/*  Written by:  Applegate, Bixby, Chvatal, and Cook                        */
/*  Date: February 24, 1995 (cofeb24)                                       */
/*                                                                          */
/*                                                                          */
/****************************************************************************/

#define CC_SAFE_MALLOC(nnum, type) \
    (type *)CCutil_allocrus(((size_t)(nnum)) * sizeof(type))

#define CC_SAFE_REALLOC(ptr, nnum, type) \
    (type *)CCutil_reallocrus((void *)ptr, ((size_t)(nnum)) * sizeof(type))

#define CC_FREE(object, type)             \
    {                                     \
        CCutil_freerus((void *)(object)); \
        object = (type *)NULL;            \
    }

#define CC_IFFREE(object, type)                \
    {                                          \
        if ((object)) CC_FREE((object), type); \
    }

void *CCutil_allocrus(size_t size);
void *CCutil_reallocrus(void *ptr, size_t size);
void CCutil_freerus(void *p);

int CCutil_reallocrus_scale(
    void **pptr, int *pnnum, int count, double scale, size_t size);
int CCutil_reallocrus_count(void **pptr, int count, size_t size);

int pmcfile_exists(const char *filename);

int pmcdir_exists(const char *dirname);

int pmcdir_create(const char *dirname);

/****************************************************************************/
/*                                                                          */
/*                             safe_io.c                                    */
/*                                                                          */
/****************************************************************************/

CC_SFILE
*CCutil_sopen(const char *f, const char *s),
    *CCutil_sdopen(int d, const char *s);

int CCutil_swrite(CC_SFILE *f, char *buf, int size),
    CCutil_swrite_bits(CC_SFILE *f, int x, int xbits),
    CCutil_swrite_ubits(CC_SFILE *f, unsigned int x, int xbits),
    CCutil_swrite_char(CC_SFILE *f, char x),
    CCutil_swrite_string(CC_SFILE *f, const char *x),
    CCutil_swrite_short(CC_SFILE *f, short x),
    CCutil_swrite_ushort(CC_SFILE *f, unsigned short x),
    CCutil_swrite_int(CC_SFILE *f, int x),
    CCutil_swrite_uint(CC_SFILE *f, unsigned int x),
    CCutil_swrite_double(CC_SFILE *f, double x),
    CCutil_sread(CC_SFILE *f, char *buf, int size),
    CCutil_sread_bits(CC_SFILE *f, int *x, int xbits),
    CCutil_sread_ubits(CC_SFILE *f, unsigned int *x, int xbits),
    CCutil_sread_char(CC_SFILE *f, char *x),
    CCutil_sread_string(CC_SFILE *f, char *x, int maxlen),
    CCutil_sread_short(CC_SFILE *f, short *x),
    CCutil_sread_ushort(CC_SFILE *f, unsigned short *x),
    CCutil_sread_short_r(CC_SFILE *f, short *x),
    CCutil_sread_int(CC_SFILE *f, int *x),
    CCutil_sread_uint(CC_SFILE *f, unsigned int *x),
    CCutil_sread_int_r(CC_SFILE *f, int *x),
    CCutil_sread_double(CC_SFILE *f, double *x),
    CCutil_sread_double_r(CC_SFILE *f, double *x), CCutil_sflush(CC_SFILE *f),
    CCutil_stell(CC_SFILE *f), CCutil_sseek(CC_SFILE *f, int offset),
    CCutil_srewind(CC_SFILE *f), CCutil_sclose(CC_SFILE *f),
    CCutil_sbits(unsigned int x), CCutil_sdelete_file(const char *fname),
    CCutil_sdelete_file_backup(const char *fname);

/****************************************************************************/
/*                                                                          */
/*                             sortrus.c                                    */
/*                                                                          */
/****************************************************************************/

void CCutil_int_array_quicksort(int *len, int n);
void CCutil_int_array_quicksort_0(int *len, int n);
void CCutil_int_perm_quicksort(int *perm, int *len, int n);
void CCutil_double_perm_quicksort(int *perm, double *len, int n);
void CCutil_int_perm_quicksort_0(int *perm, int *len, int n);

int CCutil_quickselect(int *len, int p, int n, int m);
int quickselect(int *V, int N, int k);

int bisearch(int *       sorted,
             const void *target,
             int         size,
             int (*compare)(const void *key1, const void *key2));

/****************************************************************************/
/*                                                                          */
/*                             util.c                                       */
/*                                                                          */
/****************************************************************************/

char *CCutil_strchr(char *s, int c), *CCutil_strrchr(char *s, int c),
    *CCutil_strdup(const char *s), *CCutil_strdup2(const char *s);

const char *CCutil_strchr_c(const char *s, int c),
    *CCutil_strrchr_c(const char *s, int c);

unsigned int CCutil_nextprime(unsigned int x);

int CCutil_our_gcd(int a, int b), CCutil_our_lcm(int a, int b),
    CCutil_print_command(int ac, char **av);

void CCutil_readstr(FILE *f, char *s, int len), CCutil_printlabel(void);

/* Construction of the header of the program*/
int program_header(int ac, char **av);
void dump_uname(void);

void print_line(void);

/****************************************************************************/
/*                                                                          */
/*                             zeit.c                                       */
/*                                                                          */
/****************************************************************************/

typedef struct CCutil_timer {
    double szeit;
    double cum_zeit;
    char   name[40];
    int    count;
} CCutil_timer;

double CCutil_zeit(void), CCutil_real_zeit(void),
    CCutil_stop_timer(CCutil_timer *t, int printit),
    CCutil_total_timer(CCutil_timer *t, int printit);

void CCutil_init_timer(CCutil_timer *t, const char *name),
    CCutil_start_timer(CCutil_timer *t), CCutil_suspend_timer(CCutil_timer *t),
    CCutil_resume_timer(CCutil_timer *t),
    CCutil_start_resume_time(CCutil_timer *t);
double getRealTime(void);
double getCPUTime(void);

/****************************************************************************/
/*                                                                          */
/*                             random.c                                     */
/*                                                                          */
/****************************************************************************/

int CCrandom_int_perm(int *perm, int n);
int CCrandom_dbl_perm(double *perm, int n);

/****************************************************************************/
/*                                                                          */
/*                             copy.c                                       */
/*                                                                          */
/****************************************************************************/

void fill_int(int *dst, int n, int v);
void fill_dbl(double *dst, int n, double v);
void fill_float(float *dst, int n, float v);
void fill_char(char *dst, int n, char v);
void acopy_int(const int *src, int *dst, int n);
void acopy_dbl(const double *src, double *dst, int n);

int bin_coef(int n, int r);

#ifdef __cplusplus
}
#endif

#endif
