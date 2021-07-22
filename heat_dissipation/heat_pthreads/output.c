#include "output.h"
//#include "config.h"
#include <stdio.h>

#define FPOPS_PER_POINT_PER_ITERATION (                 \
        1     /* current point 1 mul */ +               \
        4 + 1 /* direct neighbors 4 adds + 1 mul */ +   \
        4 + 1 /* diagonal neighbors 4 adds + 1 mul */ + \
        3     /* add current, direct and diagonal */ +  \
        1     /* add cumulative for average */ +        \
        1     /* difference old/new */ +                \
        1     /* comparison for maxdiff */            \
        )

void report_results(const struct parameters *p, const struct results *r)
{
    static volatile int init = 0;
    if (!init) {
        printf("Output:\n\n"
               "%13s %13s %13s %13s %13s %13s %13s\n",
               //PACKAGE_NAME, PACKAGE_VERSION, PACKAGE_BUGREPORT,
               "Iterations",
               "T(min)", "T(max)", "T(diff)", "T(avg)", "Time", "FLOP/s");
        init = 1;
    }

    printf("%-13zu % .6e % .6e % .6e % .6e % .6e % .6e\n",
           r->niter,
           r->tmin,
           r->tmax,
           r->maxdiff,
           r->tavg,
           r->time,
           (double)p->N * (double)p->M * 
           (double)(r->niter * FPOPS_PER_POINT_PER_ITERATION +
                    (double)r->niter / p->period) / r->time);
}

