/* This file is public domain. Author: Pascal Molin. */

#include <stdlib.h>
#include <stdio.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "dirichlet.h"
#include "profiler.h"

/* fast computation of modular forms coefficients.
   Use precomputed representation as Eisenstein series */

/* modular forms expressions using Eisenstein series */

struct mf_eis_desc {
    const long N;
    const long nchi;
    const long * chi;   /* nchi values mod N */
    const long degy;    /* components of final result */
    const long denom;   /* denominator to take into account */
    const long modp;    /* fft prime used for expression */
    const long ord;     /* order of root of unity */
    const long z;       /* root of unity */
    const long * E0;    /* nchi constant values of E1(chi) */
    const long * coefs; /* (nchi+1) * degy values mod modp */
};

const struct mf_eis_desc f11 = { // [ -3/2, 5/2 ]
    11, 1, (const long[]){ 10 }, 1, 2, 1, (const long[]){ 1 }, (const long[]){ -3, 5 }
};
const struct mf_eis_desc f31 = { // [3/2*y - 7/2, 1/2*y + 1/6, -y + 4/3]
    31, 2, (const long[]){ 30,6 }, 2, 6, 1, (const long[]){ -21, 9, 1, 3, 8, -6 }
};
const struct mf_eis_desc f61 = {
    61, 3, (const long[]){ 11,21,29 }, 3, 1, 580400405544961, (const long[]){ 0, 290200202772480, 0 , 2, 290200202772481, 386933603696640, 272827375838823, 563027578611304, 404306430630298, 307573029706138, 17372826933658, 369560776762983 }
};
const struct mf_eis_desc f131 = {
    131, 10, (const long[]){ 130,2,8,10,14,17,22,23,26,29 }, 10, 1, 750416685957121, (const long[]){ 370734114198171, 729894233470458, 330282884820119, 660359122210055, 480630494708254, 538607622402821, 118999898035140, 105017328772988, 160470152861453, 641017891267604 ,525286275703765, 446906851922888, 648464688534361, 722766519581883, 91830264880274, 133052133306831, 24042270749256, 105833938924033, 440199226054243, 118957193690636 ,264205262756277, 475681407279711, 649915584210309, 397245222261465, 118074290260434, 303734824077910, 273075581915145, 927237608767, 702962797799259, 248149031006752 ,676967494884818, 477065382409982, 156045808251716, 189105787097326, 398182340790368, 571153325189563, 293615647872839, 37697553012801, 559411872262861, 78824329448577 ,106554305132658, 27767371010013, 224746701016500, 241910887508324, 567833672728038, 462031815659143, 70631418062722, 13231990886753, 2216191060509, 534841720380290 ,449737772031356, 562760798893037, 151815537414119, 290614423776698, 456404230983720, 212088781194802, 270711705603041, 327393231253606, 691339992714197, 742841151715342 ,603919873770823, 331977044634075, 16397153063625, 197492068090930, 319258415610478, 314696186479672, 199900541538731, 363877897520820, 228228519741079, 305136940736957 ,385039421300682, 497579470273777, 91951008582168, 525037837304924, 428152110919880, 38090884690076, 529685856752249, 273827378660059, 644577688468547, 171653722468745 ,365885363766238, 247434318594215, 304259918194539, 552679549811891, 534428390077273, 492456707523719, 521416587571554, 510641236447450, 408983029419457, 160159060373857 ,532055837227643, 595765292650530, 219382466374002, 656695470323621, 331778353932609, 49463691438540, 113917403777178, 386518917016142, 474434463116732, 366818586730318 ,55753377972844, 574713055405017, 252901391700437, 469664967465012, 459851008366049, 459596472303613, 211144458793846, 338970471260765, 108869929728107, 410883371371910 }
};

/* use mpn_mul.c, initialize context mpn_ctx_init with my first chosen
 * 50 bits prime, then mul_mid_mpn_ctx */

/* assume deg s.t. p^deg < len */
void
_fmpz_poly_set_euler_si_mul(fmpz *z, slong len, const slong *c, slong deg, slong p)
{
    slong e, pe;
    for (e = 1, pe = p; e <= deg; e++, pe *= p)
    {
        slong m, pem, r;
        for (m = 1, pem = pe; pem < len; m++, pem += pe)
            for(r = 1; r < p && pem < len; m++, pem += pe, r++)
                fmpz_mul_si(z + pem, z+m, c[e]);
    } 
}

/* euler factor of Eisenstein series, 1/((1-chi(p)*x)*(1-p*x)) */
typedef void _fmpz_euler_func_t(fmpz *fp, slong deg, slong p, void * ctx);

/* E_2( Mod(1,N) ) as prod_p ((1-p^(-s))*(1-p^(1-s)))^(-1). Assume N prime. */
void
_fmpz_euler_factor_E2_1N(fmpz *fp, slong deg, slong p, void * ctx)
{
    long N = (long)ctx;
    if (p == N)
    {
        fmpz Q[2] = { 1, -1 };
        _fmpz_poly_inv_series(fp, Q, 2, deg+1);
    }
    else
    {
        /* (1-p*x)*(1-x) = p*x^2 + (-p - 1)*x + 1 */
        fmpz Q[3] = { 1, -p-1, p };
        _fmpz_poly_inv_series(fp, Q, 3, deg+1);
    }
}

/* E_1( (-11/.) ) as prod_p ((1-(-11/p)p^(-s))*(1-p^(-s)))^(-1) */
/* warning: ignore a0 = 1/2 */
void
_fmpz_euler_factor_E1_11(fmpz *fp, slong deg, slong p, void * ctx)
{
    slong chip = n_jacobi(-11, p);
    fmpz Q[3] = { 1, -chip-1, chip }; /* (1-chip*x)*(1-x) */
    _fmpz_poly_inv_series(fp, Q, 3, deg+1);
}

void
_fmpz_vec_linear(fmpz *z, const fmpz **x, slong len, const fmpz * coefs, slong num)
{
    slong j;
    for (j = 0; j < num; j++)
    {
        if (j == 0)
            _fmpz_vec_scalar_mul_fmpz(z, x[j], len, coefs + j);
        else
            _fmpz_vec_scalar_addmul_fmpz(z, x[j], len, coefs + j);
    }
}

void
_fmpz_poly_euler_product(fmpz *z, slong len, _fmpz_euler_func_t factor, void * ctx)
{

    n_primes_t iter;
    slong p;
    n_primes_init(iter);
    fmpz * s;
    s = _fmpz_vec_init(len);
    s[0] = 0; s[1] = 1;
    z[1] = 1;
    for (p = n_primes_next(iter); p < len; p = n_primes_next(iter))
    {
        slong deg = n_flog(len, p);
        fmpz * fp = _fmpz_vec_init(deg);
        factor(fp, deg, p, ctx);
        _fmpz_poly_set_euler_si_mul(z, len, fp, deg, p);
        _fmpz_vec_clear(fp, deg);
    }
    n_primes_clear(iter);
    _fmpz_vec_clear(s, len);
}

int main(int argc, char* argv[])
{
    slong len;

    if (argc == 2)
        len = atol(argv[1]);

    if (argc != 2 || len < 1)
    {
        flint_printf("Syntax: ell11an <integer>\n");
        flint_printf("where <integer> is the (positive) number of terms to compute\n");
        return EXIT_FAILURE;
    }

    /* euler product */
    fmpz *f, *g1, *g2, *g;
    f = _fmpz_vec_init(len);
    g1 = _fmpz_vec_init(len);
    g2 = _fmpz_vec_init(len);
    g = _fmpz_vec_init(len);

    TIMEIT_ONCE_START
    flint_printf("compute f=1*chi as Euler product 1...");
    _fmpz_poly_euler_product(f, len, _fmpz_euler_factor_E1_11, NULL);
    flint_printf("[done]\n");
    TIMEIT_ONCE_STOP

    TIMEIT_ONCE_START
    flint_printf("compute g1 = f^2 + 2f*a0...");
    _fmpz_poly_mullow(g1, f, len, f, len, len);
    _fmpz_poly_add(g1, g1, len, f, len); /* + 2*f*a0 = f */
    flint_printf("[done]\n");
    TIMEIT_ONCE_STOP

    TIMEIT_ONCE_START
    flint_printf("compute g2 = n*chi as Euler product 1...");
    _fmpz_poly_euler_product(g2, len, _fmpz_euler_factor_E2_1N, NULL);
    flint_printf("[done]\n");
    TIMEIT_ONCE_STOP

    TIMEIT_ONCE_START
    const fmpz * G[2] = { g1, g2 };
    fmpz C[2] = { 5, -3 };
    flint_printf("compute g = (5*g1 - 3*g2)/2...");
    _fmpz_vec_linear(g, G, len, C, 2);
    _fmpz_vec_scalar_tdiv_q_2exp(g, g, len, 1);
    flint_printf("[done]\n");
    TIMEIT_ONCE_STOP

    if (len <= 100)
    {
      flint_printf("Series f = 1*chi:\n");
      _fmpz_vec_print(f, len); flint_printf("\n");
      flint_printf("Series g1 = (1*chi)^2:\n");
      _fmpz_vec_print(g1, len); flint_printf("\n");
      flint_printf("Series g2 = n*chi:\n");
      _fmpz_vec_print(g2, len); flint_printf("\n");
      flint_printf("Series g = (5*g1 - 3*g2)/2 = ellan(11a1):\n");
      _fmpz_vec_print(g, len); flint_printf("\n");
    }

    _fmpz_vec_clear(f, len);
    _fmpz_vec_clear(g1, len);
    _fmpz_vec_clear(g2, len);
    _fmpz_vec_clear(g, len);
}
