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

/* factor table */
typedef struct {
    slong a;
    slong b;
    slong prev;
    slong next;
} coprime_struct;
typedef coprime_struct * coprime_ptr;
typedef const coprime_ptr coprime_srcptr;

struct rough {
    ulong m;
    struct rough * prev;
    struct rough * next;
};
//typedef struct {
//    ulong m;
//    ulong prev;
//    ulong next;
//} rough_struct;
typedef struct rough * rough_ptr;

void
pem_init(pem_ptr tab, slong len)
{
    ulong k;
    for (k = 0; k < len; k++)
        tab[k].pe = 1, tab[k].m = k;
}

/* assume tab has been inited */
void
pem_init_rough_lim(pem_ptr tab, slong lim, slong len)
{
    slong m, len1 = len / 2;
    rough_ptr rough, p1, m1;

    rough = flint_malloc((len1+1) * sizeof(struct rough));

    /* p=2 is done separately */
    //for (m = 3; m < len; m += 2)
    //    for (pe = 2, pem = 2*m; pem < len; pem <<= 1, pe <<= 1)
    //        tab[pem].pe = pe, tab[pem].m = m;
    //for (pe = 2; pe < len; pe >> 1)
    //{
    //    for (m = 3, pem = 3*pe; pem < len; pem += pe, m+=2)
    //        tab[pem].pe = pe, tab[pem].m = m;
    //}

    /* now 2-rough numbers are the odd numbers */
    rough->m = 2;
    rough->prev = NULL;
    rough->next = rough + 1;
    for (m1 = rough + 1, m = 3; m < len; m1++, m += 2)
    {
        m1->m = m;
        m1->prev = m1 - 1;
        m1->next = m1 + 1;
    }
    /* terminate */
    m1->m = len;
    m1->prev = m1 - 1;
    m1->next = NULL;

    for (p1 = rough + 1; p1->m < lim; p1 = p1->next)
    {
        slong p = p1->m, pe;
        /* skip prime powers p^e */
        for (pe = p; pe < len; pe *= p)
        {
            rough_ptr pe1 = rough + (pe>>1);
            pe1->next->prev = pe1->prev;
            pe1->prev->next = pe1->next;
        }
        /* then loop on p-rough numbers */
        for (pe = p; pe < len; pe *= p)
        {
            ulong lim = len / pe;
            /* loop m in p-rough numbers */
            for (m1 = p1->next; m1->m < lim; m1 = m1->next)
            {
                slong pem = pe * m1->m;
                tab[pem].pe = pe;
                tab[pem].m = m1->m;
                /* update links to skip pem */
                rough_ptr pem1 = rough + (pem>>1);
                pem1->next->prev = pem1->prev;
                pem1->prev->next = pem1->next;
            }
        }
    }
    /* if lim = sqrt(len), rough now links primes > lim */
    flint_free(rough);
}

void
pem_init_rough(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_rough_lim(tab, n_sqrt(len), len);
}

void
pem_init_rough_all(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_rough_lim(tab, len-1, len);
}
/* can ignore even numbers */
/* identify pmax-pem numbers in complete table of numbers less than len */
void
pem_init_sieve_range(pem_ptr tab, slong pmin, slong pmax, slong len)
{
    slong p, m;
    n_primes_t iter;
    n_primes_init(iter);
    if (pmin > 2) n_primes_jump_after(iter, pmin-1);
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                    if (tab[pem].pe == 1) tab[pem].m = m, tab[pem].pe = pe;
    }
    n_primes_clear(iter);
}
void
pem_init_sieve_all(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_sieve_range(tab, 2, len, len);
}
void
pem_init_sieve_sqrt(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_sieve_range(tab, 2, n_sqrt(len), len);
}
/* visit p^e m with m p-smooth */
void
pem_init_smooth_lim(pem_ptr tab, slong pmax, slong len)
{
    slong p, m;
    char * smooth;
    n_primes_t iter;
    n_primes_init(iter);

    smooth = flint_malloc(len * sizeof(char));
    memset(smooth, 0, len);

    smooth[1] = 1;
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        slong pe, pem, r;
        for (pe = p; pe < len; pe *= p)
            for (m = 1, pem = pe; pem < len; m++, pem += pe)
                for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                    if (smooth[m])
                    {
                        tab[pem].m = m, tab[pem].pe = pe;
                        smooth[pem] = 1;
                    }
    }
    n_primes_clear(iter);
    flint_free(smooth);
}
/* sieve up to sqrt(n) */
void
pem_init_smooth_sqrt(pem_ptr tab, slong len)
{
    pem_init(tab, len);
    pem_init_smooth_lim(tab, n_sqrt(len), len);
}
/* smooth then usual primes */
void
pem_init_smooth_all(pem_ptr tab, slong len)
{
    n_primes_t iter;
    slong p, pmax = n_sqrt(len);
    pem_init(tab, len);
    pem_init_smooth_lim(tab, pmax, len);
    n_primes_init(iter);
    n_primes_jump_after(iter, pmax);
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        ulong m, pm;
        for (m = 1, pm = p; pm < len; m++, pm += p)
        {
            tab[pm].pe = p;
            tab[pm].m = m;
        }
    }
    n_primes_clear(iter);
 }

void
pem_init_sieve_smooth2(pem_ptr tab, slong len)
{
    slong p, m, pmax = n_sqrt(len);
    char * smooth;
    n_primes_t iter;

    smooth = flint_malloc(len * sizeof(char));
    memset(smooth, 0, len);

    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;
    smooth[1] = 1;

    n_primes_init(iter);
    p = n_primes_next(iter);
    for (ulong pe = 2; pe < len; pe *= 2)
        tab[pe].pe = pe, tab[pe].m = 1, smooth[pe] = 1;
    for (p = n_primes_next(iter); p < pmax; p = n_primes_next(iter))
    {
        ulong m, pm;
        /* We must take care: m must be strictly p-smooth in
         * this loop but we declare p-smooth numbers at the
         * same time.
         * This is not a big issue when we loop on p^e then on
         * m, since the last decomposition written is the good
         * one (p^e maximal).
         * When reversing the loops we must loop backwards on m
         * to solve the problem. */
        ulong e, pe[32], lim[32];
        pe[0] = p;
        /* precompute valid exponents */
        for (e = 0; pe[e] < len; e++)
        {
            lim[e] = len / pe[e];
            pe[e+1] = pe[e] * p;
        }
        lim[e] = 0;
        /* loop */
        for (m = lim[0], pm = p*m; m; m--, pm-=p)
        {
            for (; pm && !smooth[m]; m--, pm-=p)
            /* by definition of smooth m is not divisible by p */
            smooth[pm] = 1;
            tab[pm].m = m;
            tab[pm].pe = p;
            for (e = 1; m < lim[e]; e++)
            {
               ulong pem = pe[e]*m;
               tab[pem].m = m;
               tab[pem].pe = pe[e];
            }
        }
        /* set smooth p^e last */
        for (e = 0; lim[e]; e++)
        {
            smooth[pe[e]] = 1;
        }
    }
    n_primes_clear(iter);
    flint_free(smooth);
}
/* set n=pe*m with smallest p */ 
/* all even numbers, ie 2^e * (1 mod 2)
 * then all 3^e * prime to 2*3 = 6, ie 1, 5 mod 6.
 * then all 5^e * prime to 2*3*5 = 30 ie 1, 7, 11, 13, ...
 * etc.
 * ie shifts for prime p are invertibles mod primorial p#.
 * For p=7 these are 1 and prime numbers,
 * for p=11 we also have composite values...
 *
 * Store numbers in complete table tab of length len. All numbers
 * having a prime factor > pmax can be ignored.
 *
 * FIXME: implementation below is *wrong*. One must loop
 * over p-rough numbers mod p# for residue and not only primes.
 */
void
pem_init_mod(pem_ptr tab, slong pmax, slong len)
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
/* compute only non trivial factorizations k=p^em,
   p smallest prime factor.

   Still wrong, need rough implementation.
*/
void
pem_init_mod_strict(pem_ptr tab, slong len)
{
    slong i, m, num;
    ulong mod;
    ulong pmax = n_sqrt(len);
    const ulong * prime;
    for (m = 0; m < len; m++)
        tab[m].m = tab[m].pe = 0;
    tab[1].m = tab[1].pe = 1;

    num = n_prime_pi(pmax);
    prime = n_primes_arr_readonly(pmax);

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
            slong qmax = FLINT_MIN(mod, lim1);
            ulong pem;
            /* numbers = 1 mod mod */
            for (m = len + 1, pem = pe * m; m <= lim1; m += mod, pem += pemod)
                tab[pem].pe = pe, tab[pem].m = m; 

            if (p <= 7) {
              for (m = prime[j++]; m <= qmax; m = prime[j++])
              //for (m = 1; m <= qmax; m = prime[j++])
              {
                  ulong pem = pe * m;
                  for (; m <= lim1; m += mod, pem += pemod)
                      tab[pem].pe = pe, tab[pem].m = m; 
              }
            } else {
              for (m = prime[j++]; m <= qmax; m += 2)
              {
                  ulong pem = pe * m;
                  for (; m <= lim1; m += mod, pem += pemod)
                      tab[pem].pe = pe, tab[pem].m = m; 
              }
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
            slong pem, dpe = 2*pe;
            slong lim = len / pe;
            for (m = prime[j++], pem = pe*m; m <= lim; m += 2, pem += dpe)
                tab[pem].pe = pe, tab[pem].m = m; 
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



typedef struct {
    char * name;
    void (*func)(pem_ptr, slong);
} smooth_func;

int main(int argc, char * argv[])
{
    slong i;
    int e, opt_min = 15, opt_max = 28;
    slong opt_print = 0;

#define NUM 6
    const smooth_func func[NUM] = {
        (const smooth_func){ "all", &pem_init_sieve_all },
        (const smooth_func){ "sqrt", &pem_init_sieve_sqrt },
        (const smooth_func){ "smoothtab", &pem_init_smooth_sqrt },
        (const smooth_func){ "tab all", &pem_init_smooth_all },
        (const smooth_func){ "rough", &pem_init_rough },
        (const smooth_func){ "rough all", &pem_init_rough_all }
    };

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
        slong i, j, len;
        pem_ptr tab[NUM];
        len = opt_print;

        for (i = 0; i < NUM; i++)
        {
            tab[i] = flint_malloc(len * sizeof(pem_struct));
            (func[i].func)(tab[i], len);
        }
        flint_printf("index");
        for (i = 0; i < NUM; i++)
            flint_printf("  # %9s  ", func[i].name);
        flint_printf("\n");
        for (j = 1; j < len; j++)
        {
            flint_printf("%5ld", j);
            for (i = 0; i < NUM; i++)
            {
                pem_struct pem = tab[i][j];
                if (pem.pe > 1)
                    flint_printf("  #  %3ld * %4ld", pem.pe, pem.m);
                else 
                    flint_printf("  #    . *    .");
            }
            flint_printf("\n");
        }
        for (i = 0; i < NUM; i++)
            flint_free(tab[i]);

        return 0;
    }

    flint_printf("2^k# none");
    for (i = 0; i < NUM; i++)
        flint_printf(" #  %8s ", func[i].name);
    flint_printf("\n");
    for (e = opt_min; e <= opt_max; e++)
    {
        timeit_t t0, t1;
        slong i, len = 1UL<<e;
        double ref;
        pem_ptr tab;

        timeit_start(t0);
        tab = flint_malloc(len * sizeof(pem_struct));
        pem_init(tab, len);
        flint_free(tab);
        timeit_stop(t0);
        flint_printf("%2d # %4wd", e, t0->wall);
        ref = 1. / (1 + t0->wall);
        for (i = 0; i < NUM; i++)
        {
            const smooth_func f = func[i];
            timeit_start(t1);
            tab = flint_malloc(len * sizeof(pem_struct));
            (f.func)(tab, len);
            flint_free(tab);
            timeit_stop(t1);
            flint_printf(" # %4wd [%02.1f]", t1->wall, t1->wall * ref);
        }
        flint_printf("\n");
    }
}
