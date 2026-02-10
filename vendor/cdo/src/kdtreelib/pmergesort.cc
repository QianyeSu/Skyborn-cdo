#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kdtree.h"


typedef struct param_t {
    struct kd_point *a;
    struct kd_point *b;
    size_t first;
    size_t nmemb;
    int axis;
    int max_threads;
} param_t;

extern void qsortR(void *base0, size_t n, int axis);

void pm_buildparams(struct param_t *p, struct kd_point *a, struct kd_point *b, size_t first,
                    size_t nmemb, int axis, int max_threads);
int pmergesort(struct kd_point *base, size_t nmemb, int axis, int max_threads);
void *mergesort_t(void *args);


int pmergesort(struct kd_point *base, size_t nmemb, int axis, int max_threads)
{
    struct kd_point *tmp = NULL;

    if ( max_threads > 1 )
      if ((tmp = (struct kd_point*)calloc(nmemb, sizeof(struct kd_point))) == NULL) {
        perror("malloc");
        return 0;
      }

    param_t args;
    args.a = base;
    args.b = tmp;
    args.first = 0;
    args.nmemb = nmemb;
    args.axis = axis;
    args.max_threads = max_threads;

    mergesort_t(&args);

    if (tmp) free(tmp);

    return 1;
}

void
pm_buildparams(struct param_t *p, struct kd_point *a, struct kd_point *b, size_t first,
               size_t nmemb, int axis, int max_threads)
{
    p->a = a;
    p->b = b;
    p->first = first;
    p->nmemb = nmemb;
    p->axis = axis;
    p->max_threads = max_threads;
}


void *
mergesort_t(void *args)
{
    size_t i, li, ri;
    struct param_t larg, rarg;
    pthread_t thr[2];
    param_t *mya = (param_t *) args;

    if (mya->max_threads < 2) {
        /*
         * Reached maximum number of threads allocated to this
         * branch. Proceed with sequential sort of this chunk. 
         */
        qsortR(mya->a + mya->first, mya->nmemb, mya->axis);
    } else {
        /*
         * Start two new threads, each sorting half of array a 
         */
        pm_buildparams(&larg, mya->a, mya->b, mya->first, mya->nmemb / 2, mya->axis, mya->max_threads / 2);
        /*
         * Recursively sort the left half 
         */
        if (pthread_create(&thr[0], NULL, mergesort_t, (void *) &larg)) {
            perror("pthread_create");
            return NULL;
        }

        pm_buildparams(&rarg, mya->a, mya->b, mya->first + mya->nmemb / 2,
                       mya->nmemb - mya->nmemb / 2, mya->axis, mya->max_threads / 2);
        /*
         * Recursively sort the right half 
         */
        if (pthread_create(&thr[1], NULL, mergesort_t, (void *) &rarg)) {
            perror("pthread_create");
            return NULL;
        }

        pthread_join(thr[0], NULL);
        pthread_join(thr[1], NULL);

        /*
         * Merge the two sorted chunks of array a into array b 
         */
        li = larg.first;
        ri = rarg.first;
        for(i = mya->first; i < mya->first + mya->nmemb; i++) {
            if (li >= larg.first + larg.nmemb) {
                /*
                 * We already copied everything from the left chunk,
                 * now copy from the right 
                 */
                memcpy(mya->b + i, mya->a + ri, sizeof(struct kd_point));
                ri++;
            } else if (ri >= rarg.first + rarg.nmemb) {
                /*
                 * We already copied everything from the right chunk,
                 * now copy from the left 
                 */
                memcpy(mya->b + i, mya->a + li, sizeof(struct kd_point));
                li++;
            }
            /*
             * We can still copy from both chunks, copy the smaller element 
             */
            else if ( qcmp(mya->a + li, mya->a + ri, mya->axis) < 1) {
                memcpy(mya->b + i, mya->a + li, sizeof(struct kd_point));
                li++;
            } else {
                memcpy(mya->b + i, mya->a + ri, sizeof(struct kd_point));
                ri++;
            }
        }
        /*
         * Now b is sorted, copy it back to a 
         */
        memcpy(mya->a + mya->first, mya->b + mya->first, mya->nmemb*sizeof(struct kd_point));
    }
    return NULL;
}
