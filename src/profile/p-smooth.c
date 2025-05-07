#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fft_small.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "dirichlet.h"
#include "profiler.h"

/* precomputed smoothed table */
typedef struct {
    slong m;
    slong pe;
} smooth_struct;
typedef smooth_struct * smooth_ptr;
typedef const smooth_ptr smooth_srcptr;

/* identify pmax-smooth numbers in complete table of numbers less than len */
void
smooth_table_init_sieve(smooth_ptr tab, slong pmax, slong len)
{
    slong p, m;
    n_primes_t iter;
    n_primes_init(iter);
    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;
    /* TODO: could optimize */
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                    if (tab[m].m == 0) tab[pem].m = m, tab[pem].pe = pe;
    }
    n_primes_clear(iter);
}

/* set n=pe*m with smallest p */ 
/* all even numbers, ie 2^e * (1 mod 2)
 * then all 3^e * prime to 2*3 = 6, ie 1, 5 mod 6.
 * then all 5^e * prime to 2*3*5 = 30 ie 1, 7, 11, 13, ...
 * etc.
 * ie shifts for prime p are 1 and prime numbers between p and primorial p#
 * necessary only if p < sqrt(len)
 */
void
smooth_table_init_mod(smooth_ptr tab, slong pmax, slong len)
{
    slong i, m, num;
    ulong mod = 1;
    const ulong * prime;
    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;

    num = n_prime_pi(pmax);
    prime = n_primes_arr_readonly(num);

    for (i = 0; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        mod *= p;
        for (pe = p; pe < len; pe *= p)
        {
            slong j;
            ulong m, pem, pemod = pe * mod;
            /* congruence jumps */
            m = 1; pem = pe;
            j = i + 1;
            do {
                for (; pem < len; m +=mod, pem += pemod)
                {
                    tab[pem].m = m;
                    tab[pem].pe = pe;
                }
                m = prime[j++];
                pem = pe * m;
            } while (m < mod && pem < len);
        }
    }
    n_cleanup_primes();
}

smooth_ptr
smooth_table_create(slong * size, slong pmax, slong len)
{
    slong i, t, num, alloc;
    ulong mod = 1;
    const ulong * prime;
    smooth_ptr tab;

    num = n_prime_pi(pmax);
    alloc = len;

    flint_printf("create table of length %wd\n", alloc);

    prime = n_primes_arr_readonly(num);

    t = 0;
    tab = flint_malloc(alloc * sizeof(smooth_struct));
    for (i = 0; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        mod *= p;
        for (pe = p; pe < len; pe *= p)
        {
            slong j;
            ulong m, pem, pemod = pe * mod;
            /* congruence jumps */
            m = 1; pem = pe;
            j = i + 1;
            do {
                for (; pem < len; m += mod, pem += pemod)
                    //tab[t++] = (smooth_struct){ m, pe};
                    tab[t].m = m, tab[t].pe = pe, t++;
                m = prime[j++];
                pem = pe * m;
            } while (m < mod && pem < len);
        }
    }
    n_cleanup_primes();

    flint_printf("computed table of length %wd\n", t);

    *size = t;
    //tab = flint_realloc(tab, t * sizeof(smooth_struct));
    return tab;
}
#if 0
typedef struct {
    slong m;
    slong pe;
    slong n;
} pem_struct;
typedef pem_struct * pem_ptr;
typedef const pem_ptr pem_srcptr;

void
pem_merge(pem_ptr tab, pem_ptr a, slong na, pem_ptr b, slong nb)
{
    slong ia, ib;
    for (ia = 0, ib = 0; ia < na && ib < nb; tab++)
        *tab = (a[ia].n < b[ib].n) ? a[ia++] : b[ib++];
    if (ia < na)
        for (; ia < na; tab++) *tab = a[ia++];
    else
        for (; ib < nb; tab++) *tab = b[ib++];
}

pem_ptr
smooth_table_init_pmax(slong pmax, slong len)
{
    slong p, m;
    slong size = pmax;
    pem_ptr a, b, c;
    n_primes_t iter;
    n_primes_init(iter);

    tab = (pem_ptr) flint_malloc(size * sizeof(pem_struct));

    tab[0] = (pem_struct){ 1, 1, 1};
    n0 = 1;

    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe;
        for (pe = p; pe < len; pe *= p)
        {
            for 
        }
    }


    /* numbers < pmax go to tab1, others to tab2 */
    tab1 = flint_malloc(pmax * sizeof(pem_struct));
    tab2 = flint_malloc(size * sizeof(smooth_struct));


    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
        {
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
            {
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                {
                    if (pem < pmax)
                        tab1[pem].m = m, tab1[pem].pe = pe;
                    else if 
                    if (tab[m].m) tab[pem].m = m, tab[pem].pe = pe;




    
    tab2 = flint_malloc(size * sizeof(smooth_struct);
    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;

        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
        {
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
            {
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                {
                    if 
                    if (tab[m].m) tab[pem].m = m, tab[pem].pe = pe;

                    }

    do
    {


    tab = flint_realloc(tab, 2 * size * sizeof(smooth_struct));
    for (m = size; m < 2 * size; m++)
        tab[m].m = tab[m].pe = 0;
    size *= 2;
    } while(1);

    /* TODO: could optimize */
    {
    }
    n_primes_clear(iter);

}
#endif

int main(int argc, char * argv[])
{
    slong i, len;
    ulong len_min, len_max;
    int opt_min = 15, opt_max = 28;
    slong opt_print = 0;

    /* options */
    for (i = 1; i < argc;)
    {
        if (strcmp(argv[i], "--min") == 0 && i + 1 < argc)
            opt_min = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--max") == 0 && i + 1 < argc)
            opt_max = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--print") == 0 && i + 1 < argc)
            opt_print = atol(argv[i+1]), i +=2 ;
        else break;
    }

    if (opt_print)
    {
        slong i, size;
        smooth_ptr tab1, tab2, tab3;

        len = opt_print;

        tab1 = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_sieve(tab1, n_sqrt(len), len);
        tab2 = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_mod(tab2, n_sqrt(len), len);
        tab3 = smooth_table_create(&size, n_sqrt(len), len);

        flint_printf("size: %ld\n", size);
        for (i = 0; i < 100 && i < size; i++)
        {
            slong pe = tab3[i].pe, m = tab3[i].m;
            flint_printf("%2ld: %2ld * %2ld \n", pe*m, pe, m);
        }
        for (i = 1; i < len; i++)
            flint_printf("%2ld: %2ld * %2ld, %2ld * %2ld \n",
                    i, tab1[i].pe, tab1[i].m, tab2[i].pe, tab2[i].m);

        flint_free(tab1);
        flint_free(tab2);
        flint_free(tab3);

        return 0;
    }

    len_min = 1UL << opt_min;
    len_max = 1UL << opt_max;
    
    for (len = len_min; len < len_max; len <<= 1)
    {
        timeit_t t1, t2, t3, t4;
        slong size;
        smooth_ptr tab;

        /* complete table */
        timeit_start(t1);
        tab = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_sieve(tab, len, len);
        flint_free(tab);
        timeit_stop(t1);

        /* only primes p^2 < len */
        timeit_start(t2);
        tab = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_sieve(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t2);

        /* mod table */
        timeit_start(t3);
        tab = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_mod(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t3);

        /* mod table */
        timeit_start(t3);
        tab = flint_malloc(len * sizeof(smooth_struct));
        smooth_table_init_mod(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t3);

        /* create mod table */
        timeit_start(t4);
        tab = smooth_table_create(&size, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t4);


        flint_printf("len = %9wd, %s = %3wd, %s = %3wd, %s = %3wd, %s = %3wd ms\n",
                len, "sieve all", t1->wall,
                     "sieve sqrt", t2->wall,
                     "sieve mod", t3->wall,
                     "sieve create", t4->wall);

    }
}
