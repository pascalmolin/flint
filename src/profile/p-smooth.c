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
    slong pe;
    slong m;
} pem_struct;
typedef pem_struct * pem_ptr;
typedef const pem_ptr pem_srcptr;

/* smooth number */
typedef struct {
    ulong n;
    ulong pe;
    ulong m;
} smooth_struct;
typedef smooth_struct * smooth_ptr;
typedef const smooth_ptr smooth_srcptr;

/* identify pmax-pem numbers in complete table of numbers less than len */
void
pem_table_init_sieve(pem_ptr tab, slong pmax, slong len)
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
                    if (tab[pem].m == 0) tab[pem].m = m, tab[pem].pe = pe;
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
 *
 * Store numbers in complete table tab of length len. All numbers
 * having a prime factor > pmax can be ignored.
 *
 */
void
pem_table_init_mod(pem_ptr tab, slong pmax, slong len)
{
    slong i, m, num;
    ulong mod;
    const ulong * prime;
    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;

    num = n_prime_pi(pmax);
    prime = n_primes_arr_readonly(num);

    /* first loop: use modulus */
    mod = 1;
    for (i = 0; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        if ((mod *= p) * p > len) break;
        for (pe = p; pe < len; pe *= p)
        {
            slong m, j = i + 1, pemod = pe * mod;
            slong lim1 = len / pe; /* max value for cofactor m */
            slong qmax = FLINT_MIN(lim1, FLINT_MIN(mod, pmax)); /* max prime */

            //ulong pem;
            //for (m = len + 1, pem = pe * m; m <= lim1; m += mod, pem += pemod)
            //    tab[pem].pe = pe, tab[pem].m = m; 

            //for (m = prime[j++]; m <= qmax; m = prime[j++])
            for (m = 1; m <= qmax; m = prime[j++])
            {
                ulong pem = pe * m;
                for (; m <= lim1; m += mod, pem += pemod)
                    tab[pem].pe = pe, tab[pem].m = m; 
            }
        }
    }
    /* no longer need to shift by mod */
    for (; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        for (pe = p; pe < len; pe *= p)
        {
            slong m, j = i + 1;
            slong lim = FLINT_MIN(len / pe, pmax);
            //for (m = prime[j++]; m <= lim; m = prime[j++])
            for (m = 1; m <= lim; m = prime[j++])
            {
                slong pem = pe * m;
                tab[pem].pe = pe, tab[pem].m = m; 
            }
        }
    }
    n_cleanup_primes();
}
smooth_ptr
smooth_create_pem(slong * size, pem_srcptr tab, slong len)
{
    slong s, k;
    smooth_ptr new;
    new = flint_malloc(len * sizeof(smooth_struct));
    for (k = 0, s = 0; k < len; k++)
        if (tab[k].m > 1)
            new[s++] = (smooth_struct){ .n = k, .pe = tab[k].pe, .m = tab[k].m };
    *size = s;
    new = flint_realloc(new, s * sizeof(smooth_struct));
    return new;
}

pem_ptr
pem_table_create(slong * size, slong pmax, slong len)
{
    slong i, t, num, num_primes, alloc;
    ulong mod = 1;
    const ulong * prime;
    pem_ptr tab;

    num = n_prime_pi(pmax);
    num_primes = n_prime_pi(len);
    alloc = len;

    prime = n_primes_arr_readonly(num_primes);

    t = 0;
    tab = flint_malloc(alloc * sizeof(pem_struct));
    for (i = 0; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        mod *= p;
        if (mod * p > len)
            break;
        for (pe = p; pe < len; pe *= p)
        {
            ulong m = 1, pem = pe, pemod = pe * mod;
            slong j = i + 1;
            do {
                for (; pem < len; m += mod, pem += pemod)
                    //tab[t++] = (pem_struct){ m, pe};
                    tab[t].m = m, tab[t].pe = pe, t++;
                m = prime[j++];
                pem = pe * m;
            } while (m < mod && pem < len);
        }
    }
    for (; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        for (pe = p; pe < len; pe *= p)
        {
            ulong m = 1, pem = pe;
            slong j = i + 1;
            do {
                tab[t].m = m, tab[t].pe = pe, t++;
                m = prime[j++];
                pem = pe * m;
            } while (pem < len);
        }
    }

    n_cleanup_primes();

    *size = t;
    tab = flint_realloc(tab, t * sizeof(pem_struct));
    return tab;
}

/* variant: ignore trivial decompositions p^e*1 and p*q,
 * and q is prime > sqrt(n)
 *
 * Problem: should reorder.
 */
pem_ptr
coprime_table_create(slong * size, slong pmax, slong len)
{
    slong i, t, num, num_primes, alloc;
    ulong mod;
    const ulong * prime;
    pem_ptr tab;

    num = n_prime_pi(pmax);
    num_primes = n_prime_pi(len);
    alloc = len;

    prime = n_primes_arr_readonly(num_primes);

    t = 0;
    mod = 1;
    tab = flint_malloc(alloc * sizeof(pem_struct));
    for (i = 0; i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        mod *= p;
        for (pe = p; pe < len; pe *= p)
        {
            ulong m = 1 + mod, pem = pe*m, pemod = pe * mod;
            slong j = i + 1;
            do {
                for (; pem < len; m += mod, pem += pemod)
                    tab[t].m = m, tab[t].pe = pe, t++;
                m = prime[j++];
                pem = pe * m;
            } while (m < mod && pem < len);
        }
    }
    for (; 0 && i < num; i++)
    {
        slong p, pe;
        p = prime[i];
        for (pe = p; pe < len; pe *= p)
        {
            ulong m = 1 + mod, pem = pe*m;
            slong j = i + 1;
            for (m = 1 + mod, pem = pe*m; pem < len; m = prime[j++], pem=pe*m)
                tab[t].m = m, tab[t].pe = pe, t++;
        }
    }
    n_cleanup_primes();

    *size = t;
    tab = flint_realloc(tab, t * sizeof(pem_struct));
    return tab;
}
void
pem_embed(pem_ptr tab, slong len, pem_srcptr u, slong size)
{
    slong i;
    for (i = 0; i < len; i++)
        tab[i].pe = tab[i].m = 0;
    for (i = 0; i < size; i++)
        tab[u[i].pe * u[i].m] = u[i];
}
#if 0
/* assume enough room in tab */
void
smooth_merge(pem_ptr tab, pem_srcptr v1, slong len1, pem_srcptr v2, slong len2)
{
    slong k, k1, k2, len = len1 + len2;
    for (k = 0, k1 = 0, k2 = 0; k < len; k++)
    {
        if (v1[k1].n < v2[k2].n)
            tab[k] = v1[k1++];
        else
            tab[k] = v2[k2++];
    }
}
pem_ptr
smooth_create(slong len, slong pmax)
{
    ulong num_primes;
    const ulong * prime;
    pem_ptr tab;

    tab = flint_malloc(len * sizeof(pem_struct)); 
    alt = flint_malloc(len * sizeof(pem_struct));

    num = n_prime_pi(pmax);
    s_index = flint_malloc(num * sizeof(ulong));
    prime = n_primes_arr_readonly(ulong num_primes)

    for (i = 0; i < num; i++)
    {
        ulong p = prime[i];
        new = flint_malloc(s * sizeof(pem_struct));
        for (pe = p; pe < len; pe *= p)
        {
            ulong lim = len / pe;
            for (j = 0; j < i; j++)
            {
                for (k = 0; tab[k].n < lim; k++)
                    new[n++] = (pem_struct){.n = , .pe = pe, .m = tab[k].n };
            }
        }
        alt = flint_malloc(s * sizeof(pem_struct));
        memcpy(alt, tab, s * sizeof(pem_struct));
    }


}
#endif
#if 0
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
pem_table_init_pmax(slong pmax, slong len)
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
    tab2 = flint_malloc(size * sizeof(pem_struct));


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




    
    tab2 = flint_malloc(size * sizeof(pem_struct);
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


    tab = flint_realloc(tab, 2 * size * sizeof(pem_struct));
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
        slong i, size4, size5;
        pem_ptr tab1, tab2, tab3, tab4, tab5;

        len = opt_print;

        tab1 = flint_malloc(len * sizeof(pem_struct));
        tab2 = flint_malloc(len * sizeof(pem_struct));
        tab3 = flint_malloc(len * sizeof(pem_struct));

        pem_table_init_sieve(tab1, len, len);
        pem_table_init_sieve(tab2, n_sqrt(len), len);
        pem_table_init_mod(tab3, n_sqrt(len), len);

        tab4 = pem_table_create(&size4, n_sqrt(len), len);
        tab5 = coprime_table_create(&size5, n_sqrt(len), len);

        flint_printf("size: %ld\n", size4);
        for (i = 0; i < 100 && i < size4; i++)
        {
            slong pe = tab4[i].pe, m = tab4[i].m;
            flint_printf("%2ld: %2ld * %2ld \n", pe*m, pe, m);
        }
        flint_printf("size: %ld\n", size5);
        for (i = 0; i < 100 && i < size5; i++)
        {
            slong pe = tab5[i].pe, m = tab5[i].m;
            flint_printf("%2ld: %2ld * %2ld \n", pe*m, pe, m);
        }
        for (i = 1; i < len; i++)
            flint_printf("%2ld: #  "
                                "%2ld * %2ld  #  "
                                "%2ld * %2ld  #  "
                                "%2ld * %2ld \n",
                    i, tab1[i].pe, tab1[i].m,
                       tab2[i].pe, tab2[i].m,
                       tab3[i].pe, tab3[i].m);

        flint_free(tab1);
        flint_free(tab2);
        flint_free(tab3);
        flint_free(tab4);
        flint_free(tab5);

        return 0;
    }

    len_min = 1UL << opt_min;
    len_max = 1UL << opt_max;
    
    for (len = len_min; len < len_max; len <<= 1)
    {
        timeit_t t0, t1, t2, t3, t32, t4, t5, t6;
        slong i, size;
        pem_ptr tab, queue;
        smooth_ptr squeue;

        /* init empty table */
        timeit_start(t0);
        tab = flint_malloc(len * sizeof(pem_struct));
        for (i = 0; i < len; i++)
            tab[i].pe = i, tab[i].m = 1;
        flint_free(tab);
        timeit_stop(t0);

        /* complete table */
        timeit_start(t1);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_table_init_sieve(tab, len, len);
        flint_free(tab);
        timeit_stop(t1);

        /* only primes p^2 < len */
        timeit_start(t2);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_table_init_sieve(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t2);

        /* mod table */
        timeit_start(t3);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_table_init_mod(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t3);

        /* mod table */
        timeit_start(t3);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_table_init_mod(tab, n_sqrt(len), len);
        flint_free(tab);
        timeit_stop(t3);

        /* same and keep only products */
        timeit_start(t32);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_table_init_mod(tab, n_sqrt(len), len);
        squeue = smooth_create_pem(&size, tab, len);
        flint_free(tab);
        flint_free(squeue);
        timeit_stop(t32);

        /* create mod table */
        timeit_start(t4);
        queue = pem_table_create(&size, n_sqrt(len), len);
        flint_free(queue);
        timeit_stop(t4);

        /* create mod table */
        timeit_start(t5);
        queue = coprime_table_create(&size, n_sqrt(len), len);
        flint_free(queue);
        timeit_stop(t5);

        /* create and reorder mod table */
        timeit_start(t6);
        queue = coprime_table_create(&size, n_sqrt(len), len);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_embed(tab, len, queue, size);
        flint_free(queue);
        flint_free(tab);
        timeit_stop(t6);


        flint_printf("len = %9wd, %s = %3wd, %s = %3wd, %s = %3wd, %s = %3wd, "
                      "%s = %3wd, %s = %3wd, %s = %3wd, %s = %3wd ms\n",
                len, "none", t0->wall,
                     "all", t1->wall,
                     "sqrt", t2->wall,
                     "mod", t3->wall,
                     "mod2", t32->wall,
                     "create", t4->wall,
                     "coprime", t5->wall,
                     "reorder", t6->wall);

    }
}
