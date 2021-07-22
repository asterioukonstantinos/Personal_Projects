#ifndef COMPUTE_H
#define COMPUTE_H

#include "input.h"
#include "output.h"
extern struct parameters p;
extern struct results r;

void do_compute(const struct parameters *p, 
                struct results *r);

#endif
