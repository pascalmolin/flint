/* This file is public domain. Author: Pascal Molin. */

#include <stdlib.h>
#include <stdio.h>
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "fft_small.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "dirichlet.h"
#include "profiler.h"

/* fast computation of modular forms coefficients.
   Use precomputed representation as Eisenstein series */

/* modular form expression via Eisenstein series */

struct mf_eis_desc {
    const slong N;
    const slong nchi;
    const slong * chi;   /* nchi values mod N */

    const slong degy;    /* components of final result */
    const slong * poly;  /* hecke field coefficients */

    const slong ord;     /* order of root of unity */
    const slong modp;    /* fft prime used for expression */
    const slong z;       /* root of unity for character */

    const slong denom;   /* denominator to take into account */
    const slong * E0;    /* (1+2*nchi) constant values of E1(chi)*denom */
    const slong * coefs; /* (1+nchi) * degy values mod modp */
};

/* [[1, 10], t - 1, [Mod(-1, t + 1), 2], [[1/2], [1/2]], [-3/2, 5/2]] */
const struct mf_eis_desc f11 = {
    //11, 1, (const slong[]){ 10 },
    //1, (const slong[]) { 0 },
    //1, 1, 1,
    //2,
    //(const slong[]){ 1 , 1},
    //(const slong[]){ -3, 5 }
    11, 1, (const long[]){ 10 },
    1, (const long[]){ -1 },
    2, 992515402498049, 992515402498048,
    2,
    (const long[]){ -496257701249024,-496257701249024 },
    (const long[]){ -3,5 }
};
/* [[1, 30, 6], y^2 - y - 1, [Mod(t, t^2 - t + 1), 6], [[3/2], [3/2], [Mod(-t + 2, t^2 - t + 1)], [Mod(t + 1, t^2 - t + 1)]], [-7/2, 3/2, 1/6, 1/2, 4/3, -1]] */
const struct mf_eis_desc f31 = { // [3/2*y - 7/2, 1/2*y + 1/6, -y + 4/3]
    31, 2, (const slong[]){ 30,6 },
    2, (const slong[]){ -1,-1 },
    6, 583189986803713, 490018832356491,
    6,
    (const long[]){ -291594993401855,-291594993401855,93171154447224,-93171154447221 },
    (const slong[]){ -21,9,1,3,8,-6 }
};
const struct mf_eis_desc f61 = {
    61, 3, (const slong[]){ 11,21,29 },
    3, (const slong[]){ 1,-3,-1 },
    12, 1108668498051073, 1071723438326607,
    6,
    (const long[]){ 92782715567329,-92782715567328,-439364898327501,551040210013227,-483738168588914,372062856903190 },
    (const slong[]){ 0,-3,0,12,3,-4,-111675311685725,-111675311685722,111675311685724,111675311685725,111675311685728,-111675311685726 }
};

/* store character values */
struct mf_eis_ctx {
    ulong q;
    ulong * chivec;
};
typedef struct mf_eis_ctx mf_eis_ctx_t[1];

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

void
mf_eis_ctx_init(mf_eis_ctx_t ctx, const dirichlet_group_t G, slong a, ulong ord, ulong z, nmod_t mod)
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
mf_eis_ctx_clear(mf_eis_ctx_t ctx)
{
    flint_free(ctx->chivec);
}
void
mf_eis_ctx_dual(mf_eis_ctx_t ctx, nmod_t mod)
{
    slong k;
    for (k = 0; k < ctx->q; k++)
        if (ctx->chivec[k])
            ctx->chivec[k] = nmod_inv(ctx->chivec[k], mod);
}

/* use mpn_mul.c, initialize context mpn_ctx_init with my first chosen
 * 50 bits prime, then mul_mid_mpn_ctx */

/* assume deg s.t. p^deg < len */
void
_nmod_poly_set_euler_factor(nn_ptr z, slong len, nn_srcptr c, slong deg, slong p, nmod_t mod)
{
    slong e, pe;
    for (e = 1, pe = p; e <= deg; e++, pe *= p)
    {
        slong m, pem, r;
        for (m = 1, pem = pe; pem < len; m++, pem += pe)
            for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                z[pem] = nmod_mul(z[m], c[e], mod);
    }
}

void
_nmod_poly_set_euler_factor_conj(nn_ptr z1, nn_ptr z2, slong len, nn_srcptr c1, nn_srcptr c2, slong deg, slong p, nmod_t mod)
{
    slong e, pe;
    for (e = 1, pe = p; e <= deg; e++, pe *= p)
    {
        slong m, pem, r;
        for (m = 1, pem = pe; pem < len; m++, pem += pe)
        {
            for(r = 1; r < p && pem < len; m++, pem += pe, r++)
            {
                z1[pem] = nmod_mul(z1[m], c1[e], mod);
                z2[pem] = nmod_mul(z2[m], c2[e], mod);
            }
        }
    }
}

/* euler factor of Eisenstein series, 1/((1-chi(p)*x)*(1-p*x)) */
typedef void _nmod_euler_func_t(nn_ptr fp, slong deg, slong p, void * ctx, nmod_t mod);

/* E_2( Mod(1,N) ) as prod_p ((1-p^(-s))*(1-p^(1-s)))^(-1). Assume N prime. */
void
_nmod_euler_factor_E2_1N(nn_ptr fp, slong deg, slong p, void * ctx, nmod_t mod)
{
    long N = (long)ctx;
    if (p == N)
    {
        ulong Q[2] = { 1, mod.n - 1 };
        _nmod_poly_inv_series(fp, Q, 2, deg+1, mod);
    }
    else
    {
        /* (1-p*x)*(1-x) = p*x^2 + (-p - 1)*x + 1 */
        ulong Q[3] = { 1, mod.n - p - 1, p };
        _nmod_poly_inv_series(fp, Q, 3, deg+1, mod);
    }
}

/* E_1(chi) as prod_p ((1-chi(p)p^(-s))*(1-p^(-s)))^(-1) */
void
_nmod_euler_factor_E1_chi(nn_ptr fp, slong deg, slong p, void * ctx_ptr, nmod_t mod)
{
    struct mf_eis_ctx * ctx = (struct mf_eis_ctx *)ctx_ptr;
    ulong chip = ctx->chivec[p % ctx->q];
    ulong Q[3] = { 1, mod.n - chip - 1, chip }; /* (1-chip*x)*(1-x) */
    _nmod_poly_inv_series(fp, Q, 3, deg+1, mod);
}

void
_nmod_poly_euler_product(nn_ptr z, slong len, _nmod_euler_func_t factor, void * ctx, nmod_t mod)
{
    n_primes_t iter;
    slong p;
    ulong fp[30];
    n_primes_init(iter);
    z[1] = 1;
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong deg = n_flog(len, p);
        factor(fp, deg, p, ctx, mod);
        _nmod_poly_set_euler_factor(z, len, fp, deg, p, mod);
    }
    n_primes_clear(iter);
}

slong
_nmod_poly_set_euler_factor_small(nn_ptr z, slong len, nn_srcptr c, slong deg, slong p, slong * s, slong ns, nmod_t mod)
{
    slong e, pe, ns1 = ns;
    for (e = 1, pe = p; e <= deg; e++, pe *= p)
    {
        slong i;
        for (i = 0; i < ns; i++)
        {
            slong m = s[i], pem = pe*m;
            if (pem > len) continue;
            z[pem] = nmod_mul(z[m], c[e], mod);
            s[ns1++] = pem;
        }
    }
    return ns1;
}

void
_nmod_poly_euler_product_stages(nn_ptr z, slong len, _nmod_euler_func_t factor, void * ctx, nmod_t mod)
{
    n_primes_t iter;
    slong p, B, ns;
    slong * s;
    ulong fp[30];
    n_primes_init(iter);

    /* number of smooth indices stored */
    B = n_cbrt(len);

    s = (slong *)flint_malloc(len * sizeof(slong));
    z[1] = 1;

    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong deg = n_flog(len, p);
        factor(fp, deg, p, ctx, mod);
        if (p < B && ns < len)
            ns = _nmod_poly_set_euler_factor_small(z, len, fp, deg, p, s, ns, mod);
        else
            _nmod_poly_set_euler_factor(z, len, fp, deg, p, mod);
    }
    n_primes_clear(iter);
    flint_free(s);
}

void
_nmod_poly_euler_product_conj(nn_ptr z1, nn_ptr z2, slong len, _nmod_euler_func_t factor, void * ctx, void * ctx2, nmod_t mod)
{
    n_primes_t iter;
    slong p;
    ulong fp1[30], fp2[30];
    n_primes_init(iter);

    z1[1] = z2[1] = 1;
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong deg = n_flog(len, p);
        factor(fp1, deg, p, ctx, mod);
        factor(fp2, deg, p, ctx2, mod);
        _nmod_poly_set_euler_factor_conj(z1, z2, len, fp1, fp2, deg, p, mod);
    }
    n_primes_clear(iter);
}
/* currently only output constant coefficients mod poly */
void
_fmpz_modular_form_expansion(fmpz * a, slong len, const struct mf_eis_desc f)
{
    nmod_t mod;
    dirichlet_group_t G;
    slong k;
    ulong pow2 = n_clog2(len), len2;
    nn_ptr g, g1, g2, g12;
    mpn_ctx_t fft_ctx;

    FLINT_ASSERT(n_trailing_zeros(f.modp-1) > pow2);
    nmod_init(&mod, f.modp);
    mpn_ctx_init(fft_ctx, mod.n);
    len2 = 1UL << pow2;
    len2 = len;

    g = _nmod_vec_init(len2);

    /* initialize with c[0] * E_2(1 mod N) */
    {
        ulong c = nmod_set_si(f.coefs[0], mod);
        g[0] = 0;
        _nmod_poly_euler_product(g, len2, _nmod_euler_factor_E2_1N, (void *)f.N, mod);
        _nmod_vec_scalar_mul_nmod(g, g, len2, c, mod);
    }

    /* then add products c[k] * E1(chi) * E1(chi^-1) */

    dirichlet_group_init(G, f.N);
    g1 = _nmod_vec_init(len2);
    g2 = _nmod_vec_init(len2);
    g12 = _nmod_vec_init(len2);
    for (k = 0; k < f.nchi; k++)
    {
        mf_eis_ctx_t ctx, ctx2;
        ulong c = nmod_set_si(f.coefs[(k+1)*f.degy], mod);

        TIMEIT_ONCE_START
        flint_printf("[ char %ld ]\n",k+1);
        mf_eis_ctx_init(ctx, G, f.chi[k], f.ord, f.z, mod);
        TIMEIT_ONCE_STOP

        TIMEIT_ONCE_START
        flint_printf("##    euler g1 and g2...", k+1);
        mf_eis_ctx_init(ctx2, G, f.chi[k], f.ord, nmod_inv(f.z,mod), mod);
        _nmod_poly_euler_product_conj(g1, g2, len2, _nmod_euler_factor_E1_chi, ctx, ctx2, mod);
        flint_printf("[done]\n");
        TIMEIT_ONCE_STOP

        TIMEIT_ONCE_START
        flint_printf("##    euler g1 with stages...", k+1);
        g1[0] = nmod_set_si(f.E0[2*k], mod);
        _nmod_poly_euler_product_stages(g1, len2, _nmod_euler_factor_E1_chi, ctx, mod);
        flint_printf("[done]\n");
        TIMEIT_ONCE_STOP

        TIMEIT_ONCE_START
        flint_printf("##    euler g2 direct...");
        mf_eis_ctx_dual(ctx, mod);
        g2[0] = nmod_set_si(f.E0[2*k+1], mod);
        _nmod_poly_euler_product(g2, len2, _nmod_euler_factor_E1_chi, ctx, mod);
        flint_printf("[done]\n");
        TIMEIT_ONCE_STOP

        mf_eis_ctx_clear(ctx);

#if 0
        TIMEIT_ONCE_START
        flint_printf("           default product g1 * g2...");
        _nmod_poly_mullow(g12, g1, len2, g2, len2, len2, mod);
        flint_printf("[done]\n");
        TIMEIT_ONCE_STOP
#else
        /* force fft_small to use prime mod.n */
        TIMEIT_ONCE_START
        flint_printf("           hacked product g1 * g2...");
        _nmod_poly_mul_mid_mpn_ctx(g12, 0, len2, g1, len2, g2, len2, mod, fft_ctx);
        flint_printf("[done]\n");
        TIMEIT_ONCE_STOP
#endif
        _nmod_vec_scalar_addmul_nmod(g, g12, len2, c, mod);

    }
    _nmod_vec_clear(g12);
    _nmod_vec_clear(g1);
    _nmod_vec_clear(g2);
    mpn_ctx_clear(fft_ctx);

    g[0] = 0;
    if (f.denom > 1)
    {
        ulong inv = nmod_inv(f.denom, mod);
        _nmod_vec_scalar_mul_nmod(g, g, len2, inv, mod);
    }
    _fmpz_vec_set_nmod_vec(a, g, len, mod);

    _nmod_vec_clear(g);
    dirichlet_group_clear(G);
}

int main(int argc, char* argv[])
{
    slong N, len;
    fmpz * a;

    if (argc == 3)
    {
        N = atol(argv[1]);
        len = atol(argv[2]);
    }

    if (argc != 3 || len < 1)
    {
        flint_printf("Syntax: mfcoefs <level> <length>\n");
        flint_printf("where <length> is the (positive) number of terms to compute\n");
        return EXIT_FAILURE;
    }

    a = _fmpz_vec_init(len);

    if (N == 11)
        _fmpz_modular_form_expansion(a, len, f11);
    else if (N == 31)
        _fmpz_modular_form_expansion(a, len, f31);
    else if (N == 61)
        _fmpz_modular_form_expansion(a, len, f61);
    else
        return EXIT_FAILURE;

    if (len < 1000)
    {
        _fmpz_vec_print(a, len); flint_printf("\n");
    }
    else
    {
        _fmpz_vec_print(a, 100); flint_printf(" [...] \n");
        _fmpz_vec_print(a + len - 101, 100); flint_printf("\n");
    }

    _fmpz_vec_clear(a, len);
}
