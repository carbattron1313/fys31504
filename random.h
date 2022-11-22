//
//  random.h
//  
//
//  Created by Morten Hjorth-Jensen
//  Edited by Elena Muñoz, David Martínez, Antonio Gómez and Alejandro Carballido.
//

#ifndef random_h
#define random_h

using namespace std;

/*
** The function
**           ran1()
** is an "Minimal" random number generator of Park and Miller
** (see Numerical recipe page 280) with Bays-Durham shuffle and
** added safeguards. Call with idum a negative integer to initialize;
** thereafter, do not alter idum between sucessive deviates in a
** sequence. RNMX should approximate the largest floating point value
** that is less than 1.
** The function returns a uniform deviate between 0.0 and 1.0
** (exclusive of end-point values).
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran(long *seed)
{
    int             j;
    long            k;
    static long     iy=0;
    static long     iv[NTAB];
    double          temp;

    if (*seed <= 0 || !iy) {
        if (-(*seed) < 1) *seed=1;
        else *seed = -(*seed);
        for(j = NTAB + 7; j >= 0; j--) {
            k     = (*seed)/IQ;
            *seed = IA*(*seed - k*IQ) - IR*k;
            if(*seed < 0) *seed += IM;
            if(j < NTAB) iv[j] = *seed;
        }
        iy = iv[0];
    }
    k     = (*seed)/IQ;
    *seed = IA*(*seed - k*IQ) - IR*k;
    if(*seed < 0) *seed += IM;
    j     = iy/NDIV;
    iy    = iv[j];
    iv[j] = *seed;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// End: function ran()

