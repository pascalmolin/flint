/* This file is public domain. Author: Pascal Molin. */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_mat.h"
#include "fft_small.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "dirichlet.h"
#include "profiler.h"

/*
 Fast computation of modular forms coefficients via
 representation as products of Eisenstein series.

 A modular form f in S_2(N) is given as

 f = c_0 E2(N) + sum_{i=1}^{n} c_i E1(chi_i)E1(chi_i^-1)

 with E2(N) the usual Eisenstein series of weight 2 and level N
 and E1(chi) the weight one Eisenstein series of character chi.


 The coefficients of Eisenstein series are computed as lines in
 a nmod matrix (the modulus having enough roots of unity).

 For big lengths this matrix may be compressed to keep only
 prime indices.

 E1 a2 a3 a5 ...
 E2 a2 a3 a5 ...

 A form is a combination of Ei whose coefficients belong
 to some number ring.

 It is given as a matrix whose rows correspond to integral
 basis (and columns to generators).

 When several forms are required the matrices may be
 concatenated.

 A form is then given by a matrix whose lines correspond to an
 integral basis of the Hecke field.

 We choose to transpose the output so that each line corresponds
 to a Fourier coefficient. By default only coefficients of
 prime index are output.
*/

/* data format, one character */
struct mf_eis_space {
    const slong N;       // level
    const slong k;       // weight
                         // TODO: space character
    const slong nchi;    // number of Dirichlet character used
    const slong * chi;   // characters by Conrey index mod N
    
    const slong ord;     // order of root of unity
    const ulong modp;    // fft prime used for expression
    const ulong z;       // root of unity for character

    const slong num;     // number of Eisenstein generators
    const slong * l;     // weight
    const slong * c;     // character index (-1 for Ek)
    const slong * d;     // Bd operator
    const ulong * e0;    // constant terms (two)
    
    const slong rank;    // rank of output basis
    const ulong * basis; // conversion matrix from generators to basis
                         // rank * num, could be nmod_mat

                         // FIXME: skip this, depends on form
    const slong dim;     // number of actual forms
    const slong * deg;   // degree of form
    const char ** poly;  // Hecke polynomial
    const char ** zk;    // integral basis used as string
};
typedef struct mf_eis_space mf_space_t[1];

// [11, [ [2, Mod(1,11)], -3/2; [1, Mod(-1,11), 1, -1], 5/2 ], [-1,2] ]
const struct mf_eis_space mf11 = {
    11, 2,

    2, (const slong[]){ 1, 10 },

    2, 928284166586369, 928284166586368,

    2,
    (const slong[]){ 2, 1 },
    (const slong[]){ 0, 1 },
    (const slong[]){ 1, 1 },
    (const ulong[]){ 0, 0, 464142083293185,464142083293185 },

    1,
    (const ulong[]){ 464142083293183,464142083293187 },

    1,
    (const slong[]){ 1 },
    (const char*[]){ "y-1" },
    (const char*[]){ "[1]" }
};
const struct mf_eis_space mf23 = {
    23, 2,
    3, (const slong[]){ 1,5,22 },
    22, 570715254292481, 165626152232624,
    3,
    (const slong[]){ 2,1,1 }, // weight
    (const slong[]){ 1,3,2 }, // index
    (const slong[]){ 1,1,1 }, // d
    (const ulong[]){ 0,0,1,1,87660967136177,153751443341630 }, // e[0]
    3, // rank
    (const ulong[]){ 73563789228892,64666259459551,216951869898920,205217271741967,195003784695111,556419154499629 },
    1, // forms
    (const slong[]){ 2 }, // hecke degree
    (const char *[]){ "y^2 - y - 1" }, // hecke pol
    (const char *[]){ "[1, y]" } // basis
};
/* [[1, 30, 6], y^2 - y - 1, [Mod(t, t^2 - t + 1), 6], [[3/2], [3/2], [Mod(-t + 2, t^2 - t + 1)], [Mod(t + 1, t^2 - t + 1)]], [-7/2, 3/2, 1/6, 1/2, 4/3, -1]] */
//const struct mf_eis_space mf31 = {
//    31, 2,
//
//    3, (const slong[]){ 1, 30, 6 },
//
//    6, 905789275373569, 644816079264124,
//
//
//    3,
//    (const slong[]){ 2, 1, 1}, // weight l
//    (const slong[]){ 0, 1, 2}, // character index
//    (const slong[]){ 1, 1, 1}, // expansion
//    (const ulong[]){ 0, 0, 452894637686786,452894637686786,260973196109447,644816079264125 },
//
//    2, 
//    (const ulong[]){ 452894637686781,452894637686786,754824396144641,452894637686785,603859516915714,905789275373568 },
//
//    1,
//    (const slong[]){ 1 },
//    (const char*[]){ "y^2-y-1" },
//    (const char*[]){ "[1,y]" }
//};
//  [41, [ [2,Mod(1,41)] ,  24071217223589*y^2 + 251361673798276*y + 21090979747625;
//         [1, Mod(3,41), 1, -1] , 301634451880304*y^2 + 339546933529028*y + 40434702472896;
//         [1, Mod(6,41), 1, -1], 383565831106330*y^2 + 441101456979627*y + 216060380163714;
//         [1, Mod(11,41), 1, -1], 246106924681630*y^2 + 520635285518189*y + 91457017589189 ]
//       , [Mod(227071490884881, 587207928709121), 40], y^3 - y^2 - 3*y + 1 ],
const struct mf_eis_space mf41 = {
    41, 2,
    4, (const slong[]){ 1,3,6,11 },
    40, 587207928709121, 227071490884881,
    4,
    (const slong[]){ 2,1,1,1 }, // weight
    (const slong[]){ 1,2,3,4 }, // index
    (const slong[]){ 1,1,1,1 }, // d
    (const ulong[]){ 1,1,127699192864814,556572729132039,580642649055909,377124033433130,464057956047059,214951021409062 }, // e[0]
    4, // rank
    (const ulong[]){ 69233414194803,275432891021865,24071217223589,56495677524383,53973456700211,301634451880304,395984113667253,237459359376836,383565831106330,583670866952449,179534281490698,246106924681630 },
    1, // forms
    (const slong[]){ 3 }, // hecke degree
    (const char *[]){ "y^3 - y^2 - 3*y + 1" }, // hecke pol
    (const char *[]){ "[1, y, y^2 - y - 2]" } // basis
};
// [43, [ [2, Mod(1, 43)], 1/2*y-1/2;
//        [1, Mod(42, 43), 1, -1], y/2+5/6;
//        [1, Mod(7, 43), 1, -1], -y+2/3
//        ], [Mod(t, t^2 - t + 1), 6], y^2-2 ]
//const struct mf_space_t mf43 = {
//    43, 2,
//    3, (const slong[]){ 1, 42, 7 },
//    2, 928284166586369, 928284166586368,
//
//    3,
//    (const slong[]){ 2, 1 },
//    (const slong[]){ 0, 1 },
//    (const slong[]){ 1, 1 },
//    (const ulong[]){ 0, 0, 464142083293185,464142083293185 },
//
//    1,
//    (const ulong[]){ 464142083293183,464142083293187 }
//};

/* character values mod N */
struct mf_char_ctx {
    ulong q;
    ulong * chivec;
};
typedef struct mf_char_ctx mf_char_ctx_t[1];
typedef struct mf_char_ctx * mf_char_ctx_ptr;

/* Compute the list of character values as powers of z mod p,
 * assume z has exact order ord */
void
dirichlet_chi_vec_nmod(ulong *v, slong nv, const dirichlet_group_t G,
        const dirichlet_char_t chi, ulong ord, ulong z, nmod_t mod)
{
    slong k;
    dirichlet_chi_vec_order(v, G, chi, ord, nv);
    for (k = 0; k < nv; k++)
    {
        if (v[k] == DIRICHLET_CHI_NULL)
            v[k] = 0;
        else
            v[k] = nmod_pow_ui(z, v[k], mod);
    }
}

/* store values chi(k) */
void
mf_char_ctx_init(mf_char_ctx_t ctx, const dirichlet_group_t G, slong a, ulong ord, ulong z, nmod_t mod)
{
    dirichlet_char_t chi;
    dirichlet_char_init(chi, G);
    dirichlet_char_log(chi, G, a);
    ctx->q = G->q;
    ctx->chivec = (ulong *)flint_malloc(G->q*sizeof(ulong));
    dirichlet_chi_vec_nmod(ctx->chivec, G->q, G, chi, ord, z, mod);
    dirichlet_char_clear(chi);
}
void
mf_char_ctx_clear(mf_char_ctx_t ctx)
{
    flint_free(ctx->chivec);
}

/* change chi -> chi^(-1) */
void
mf_char_ctx_dual(mf_char_ctx_t ctx, nmod_t mod)
{
    slong k;
    for (k = 0; k < ctx->q; k++)
        if (ctx->chivec[k])
            ctx->chivec[k] = nmod_inv(ctx->chivec[k], mod);
}

/* precompute coprime decomposition k = p^e*m, p smallest prime */
typedef struct {
    slong pe;
    slong m;
} pem_struct;
typedef pem_struct * pem_ptr;
typedef const pem_ptr pem_srcptr;
struct rough {
    ulong m;
    struct rough * prev;
    struct rough * next;
};
typedef struct rough * rough_ptr;

void
pem_init_rough(pem_ptr tab, slong lim, slong len)
{
    ulong m, pe, pem, len1 = (len-1) / 2;
    rough_ptr rough, p1, m1;

    /* p=2 can be done separately */
    for (m = 3; m <= len1; m += 2)
        for (pe = 2, pem = 2*m; pem < len; pem <<= 1, pe <<= 1)
            tab[pem].pe = pe, tab[pem].m = m;

    /* now need 2-rough (ie odd) numbers up to len / 3 */
    len1 = (len-1) / 3;
    rough = flint_malloc((2 + len1/2) * sizeof(struct rough));
    rough->m = 1;
    rough->prev = NULL;
    rough->next = rough + 1;
    for (m1 = rough + 1, m = 3; m <= len1; m1++, m += 2)
    {
        m1->m = m;
        m1->prev = m1 - 1;
        m1->next = m1 + 1;
    }
    /* terminate */
    m1->m = len;
    m1->prev = m1 - 1;
    m1->next = NULL;

    for (p1 = rough + 1; p1->m <= lim; p1 = p1->next)
    {
        slong p = p1->m, pe;
        /* skip prime powers p^e up to len1 */
        for (pe = p; pe <= len1; pe *= p)
        {
            rough_ptr pe1 = rough + (pe>>1);
            pe1->next->prev = pe1->prev;
            pe1->prev->next = pe1->next;
        }
        /* then loop on p-rough numbers (p>=3) */
        for (pe = p; pe <= len1; pe *= p)
        {
            ulong lim = (len-1) / pe;
            /* loop m in p-rough numbers */
            for (m1 = p1->next; m1->m <= lim; m1 = m1->next)
            {
                slong pem = pe * m1->m;
                tab[pem].pe = pe;
                tab[pem].m = m1->m;
                /* update links to skip pem */
                if (pem >= len1) continue;
                rough_ptr pem1 = rough + (pem>>1);
                pem1->next->prev = pem1->prev;
                pem1->prev->next = pem1->next;
            }
        }
    }
    /* if lim = sqrt(len), rough now links primes > lim */
    flint_free(rough);
}

typedef struct {
    ulong n;
    ulong a;
    ulong b;
} coprime_t;
typedef coprime_t * coprime_ptr;

coprime_ptr
coprime_table_init(slong * size, slong len)
{
    slong k, n;
    coprime_ptr fac;
    pem_ptr tab = flint_malloc(len * sizeof(pem_struct));
    for (k = 0; k < len; k++)
        tab[k].pe = tab[k].m = 0;
    pem_init_rough(tab, n_sqrt(len), len);
    fac = flint_malloc(len * sizeof(coprime_t));
    for (n = 0, k = 1; k < len; k++)
        if(tab[k].m)
            fac[n++] = (coprime_t){ .n = k, .a = tab[k].pe, .b = tab[k].m };
    flint_free(tab);
    *size = n;
    fac = flint_realloc(fac, n * sizeof(coprime_t));
    return fac;
}

/* Eisenstein series: direct algorithm */
void
_nmod_poly_eisenstein_series(nn_ptr z, slong len, slong k, mf_char_ctx_t psi, const coprime_ptr tab, slong size, nmod_t mod)
{
    slong p;
    n_primes_t iter;
    n_primes_init(iter);
    z[1] = 1;
    /* first expand Euler factors */
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong pe, pe1 = 1;
        ulong cp = psi ? psi->chivec[p % psi->q] : 1;
        if (k > 1)
            cp = nmod_mul(cp, nmod_pow_ui(p, k-1, mod), mod);
        z[p] = nmod_add(cp, 1, mod);
        for (pe = p, pe1 = 1; pe < len; pe1 = pe, pe *= p)
            z[pe] = nmod_add(nmod_mul(z[pe1], cp, mod), 1, mod);
    }
    /* then fill composite */
    for (k = 0; k < size; k++)
        z[tab[k].n] = nmod_mul(z[tab[k].a], z[tab[k].b], mod);
}

/* Modular form */
void
nmod_vec_set_primes(nn_ptr a, nn_srcptr g, slong len)
{
    n_primes_t iter;
    slong j, p;

    /* assume a and g large enough */
    n_primes_init(iter);
    for (j = 0, p = n_primes_next(iter); p < len; j++, p = n_primes_next(iter))
        a[j] = g[p];
    n_primes_clear(iter);
}

void
nmod_mat_modular_form_series(nmod_mat_t a, const mf_space_t mf, slong len)
{
    nmod_t mod;
    coprime_ptr tab = NULL;
    dirichlet_group_t G;
    mf_char_ctx_ptr char_ctx;
    nn_ptr g1, g2, g12;
    nmod_mat_t eis, basis;
    mpn_ctx_t fft_ctx;
    const ulong * b;
    slong size, cols, i, j;
    /* init */
    
    nmod_mat_set_mod(a, mf->modp);
    nmod_init(&mod, mf->modp);

    /* sanity checks */
    cols = n_prime_pi(len);
    FLINT_ASSERT(n_prime_pi(len) == nmod_mat_ncols(a));
    FLINT_ASSERT(mf->rank == nmod_mat_nrows(a));
    FLINT_ASSERT(n_trailing_zeros(f.modp-1) > n_clog2(len));
    
    /*
     Critical part: force mpn_mul (_nmod_poly_mul_mid_mpn_ctx)
     to use custom 50 bits prime p = mod.n
     Initialize context fft_ctx accordingly.
     */
    mpn_ctx_init(fft_ctx, mod.n);

    /* tabulate composite */
    tab = coprime_table_init(&size, len);

    /* precompute chars */
    char_ctx = flint_malloc(mf->nchi * sizeof(struct mf_char_ctx));
    dirichlet_group_init(G, mf->N);
    for (i = 0; i < mf->nchi; i++)
        mf_char_ctx_init(char_ctx + i, G, mf->chi[i], mf->ord, mf->z, mod);
    dirichlet_group_clear(G);

    /* compute eisenstein expansions */
    nmod_mat_init(eis, mf->num, cols, mf->modp);
    g1 = _nmod_vec_init(len);
    g2 = _nmod_vec_init(len);
    g12 = _nmod_vec_init(len);
    for (i = 0; i < mf->num; i++)
    {
        slong k = mf->l[i], c = mf->c[i];
        nn_ptr row = nmod_mat_entry_ptr(eis, i, 0);
        mf_char_ctx_ptr psi = (c == -1) ? NULL : char_ctx + mf->c[i];
        if (k < mf->k)
        {
            _nmod_poly_eisenstein_series(g1, len, k, psi, tab, size, mod);
            mf_char_ctx_dual(psi, mod);
            _nmod_poly_eisenstein_series(g2, len, k, psi, tab, size, mod);
            mf_char_ctx_dual(psi, mod);
            g1[0] = nmod_set_ui(mf->e0[2*i], mod);
            g2[0] = nmod_set_ui(mf->e0[2*i+1], mod);
            _nmod_poly_mul_mid_mpn_ctx(g12, 0, len, g1, len, g2, len, mod, fft_ctx);
            g12[0] = 0;
            nmod_vec_set_primes(row, g12, len);
        }
        else
        {
            /* FIXME: need only prime indices */
            _nmod_poly_eisenstein_series(g12, len, k, psi, tab, size, mod);
            nmod_vec_set_primes(row, g12, len);
        }
    }
    _nmod_vec_clear(g1);
    _nmod_vec_clear(g2);
    _nmod_vec_clear(g12);

    flint_free(tab);

    /* convert to basis */
    nmod_mat_init(basis, mf->rank, mf->num, mf->modp);
    for (b = mf->basis, i = 0; i < basis->r; i++)
        for (j = 0; j < basis->c; j++, b++)
            nmod_mat_entry(basis, i, j) = *b;

    nmod_mat_mul(a, basis, eis);

    nmod_mat_clear(basis);
    nmod_mat_clear(eis);
}
void
fmpz_mat_set_transpose_nmod_mat(fmpz_mat_t b, const nmod_mat_t a)
{
    slong i, j;

    for (i = 0; i < a->r; i++)
        for (j = 0; j < a->c; j++)
            fmpz_set_ui_smod(fmpz_mat_entry(b, j, i),
                             nmod_mat_entry(a, i, j), a->mod.n);
}

int usage(int count, const char * fname[])
{
    int i;
    flint_printf("mfcoefs [options] <level> <length>\n");
    flint_printf("where <level> is in");
    for (i = 0; i < count; i++) flint_printf(" %s,", fname[i]);
    flint_printf("\n and <length> is the number of terms to compute\n");
    flint_printf("output coefficients as a matrix, one row per coefficient a_p,");
    flint_printf(" columns indexed by an integral basis of the value field\n");
    flint_printf("options:\n");
    flint_printf(" --raw: raw flint output (matrix size followed by space separated values)\n");
    flint_printf(" --all: all coefficients a_n (default only a_p)\n");
    flint_printf(" --tail <n>: output only last <n> coefficients (implies --all)\n");
    flint_printf(" --time: time each step (implies --tail 0)\n");
    return EXIT_FAILURE;
}

int main(int argc, char * argv[])
{
    slong i, len;

    nmod_mat_t a;
    fmpz_mat_t m;
    slong count = 3, cols;
    const char * mf_name[]       = { "11", "23", "41" };
    const struct mf_eis_space *f = NULL, mf[] = { mf11 , mf23 , mf41 };
    int opt_all = 0, opt_raw = 0, opt_time = 0, opt_smooth = 0, opt_test = 0;
    long opt_tail = -1;

    /* options */
    for (i = 1; i < argc;)
    {
        if (strcmp(argv[i], "--all") == 0)
            opt_all = 1, i++;
        else if (strcmp(argv[i], "--raw") == 0)
            opt_raw = 1, i++;
        else if (strcmp(argv[i], "--tail") == 0 && i + 1 < argc)
            opt_tail = atol(argv[i+1]), i +=2 ;
        else if (strcmp(argv[i], "--time") == 0)
            opt_time = 1, i++;
        else if (strcmp(argv[i], "--smooth") == 0)
            opt_smooth = 1, i++;
        else if (strcmp(argv[i], "--smooth2") == 0)
            opt_smooth = 2, i++;
        else if (strcmp(argv[i], "--test") == 0)
            opt_test = 1, i++;
        else break;
    }

    if (opt_tail >= 0) opt_all = 1;

    /* form and length */
    if (argc == i + 2)
    {
        slong j;
        for (j = 0; j < count; j++)
           if (strcmp(argv[i], mf_name[j]) == 0)
               break;
        f = mf + j;
        len = atol(argv[i+1]);
    }

    if (argc != i + 2 || len < 1 || f == NULL)
        return usage(count, mf_name);

    cols = n_prime_pi(len);
    nmod_mat_init(a, f->rank, cols, f->modp);

    nmod_mat_modular_form_series(a, f, len);

    if (opt_time)
        return 0;

    fmpz_mat_init(m, cols, f->rank);
    fmpz_mat_set_transpose_nmod_mat(m, a);

    if (opt_raw)
        fmpz_mat_print(m);
    else
        fmpz_mat_print_pretty(m);
    flint_printf("\n");

    fmpz_mat_clear(m);
}
