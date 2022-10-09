
#ifndef GMX_MATH_VEC_H
#define GMX_MATH_VEC_H

#include "functions.h"

static inline void dvec_add(const double* a, const double* b, double* c)
{
    double x, y, z;

    x = a[XX]+b[XX];
    y = a[YY]+b[YY];
    z = a[ZZ]+b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void dvec_sub(const double* a, const double* b, double* c)
{
    double x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline void dsvmul(double a, const double* v1, double* v2)
{
    v2[XX] = a*v1[XX];
    v2[YY] = a*v1[YY];
    v2[ZZ] = a*v1[ZZ];
}

static inline double diprod(const double* a, const double* b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}

static inline double dnorm2(const double* a)
{
    return a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ];
}

static inline double dnorm(const double* a)
{
    return sqrt(a[XX]*a[XX]+a[YY]*a[YY]+a[ZZ]*a[ZZ]);
}


static inline void dcprod(const double* a, const double* b, double* c)
{
    c[XX] = a[YY]*b[ZZ]-a[ZZ]*b[YY];
    c[YY] = a[ZZ]*b[XX]-a[XX]*b[ZZ];
    c[ZZ] = a[XX]*b[YY]-a[YY]*b[XX];
}

static inline double dvang(double* a, double* b)
{
    double w[3];
    double wlen, s;

    dcprod(a, b, w);
    wlen  = dnorm(w);
    s     = diprod(a, b);

    return atan2(wlen, s);
}


#endif

