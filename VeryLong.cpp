#include "VeryLong.h"
#include "VeryLongModular.h"
#include <time.h>
#include "pow.h"
#include <limits.h>
#ifdef WIN32
#define LONG_LONG_MAX LLONG_MAX
#endif
#ifdef __SUNPRO_CC
#define LONG_LONG_MAX INT64_MAX
#endif
#include <cstdlib>
using namespace std;

namespace
{
bool Debug = false;
}

void VeryLong::setDebug()
{
    Debug = true;
}

void VeryLong::clearDebug()
{
    Debug = false;
}

extern "C"
{
#include "lip.h"
}
#define USE_LIP 1

// Stolen from lip.c (c) A. Lenstra
// Converted to GMP by C.T. Card 2002

#if 0
const int NBITS = 30;
const int NBITSH = (NBITS>>1);
const int RADIX = (1<<NBITS);
#endif

#ifdef USE_LIP
long mpz_sqrts(long n)
{
    long a;
    long ndiva;
    long newa;
    mpz_t ln;
    mpz_t rr;
    mpz_t dif;
    mpz_init(ln);
    mpz_init(rr);
    mpz_init(dif);

    if (n <= 0) return (0);
    if (n <= 3) return (1);
    if (n <= 8) return (2);
    if (n >= RADIX)
    {
        long ret;

        mpz_set_ui(ln, n);
        mpz_sqrtrem(rr, dif, ln);
        ret = mpz_get_ui(rr);
        mpz_clear(rr);
        mpz_clear(ln);
        mpz_clear(dif);
        return ret;
    }
    mpz_clear(rr);
    mpz_clear(ln);
    mpz_clear(dif);
    newa = 3L << (2 * (NBITSH - 1));
    a = 1L << NBITSH;
    while (!(n & newa))
    {
        newa >>= 2;
        a >>= 1;
    }
    while (1)
    {
        newa = ((ndiva = n / a) + a) / 2;
        if (newa - ndiva <= 1)
        {
            if (newa * newa <= n)
                return (newa);
            else
                return (ndiva);
        }
        a = newa;
    }
}
#endif

/* for small prime genaration */
#if 0
const int PRIM_BND = 16500;
const int PRIM_UP = ((((PRIM_BND<<1)+1)*((PRIM_BND<<1)+1))-(NBITS<<2));
#endif
#ifdef USE_LIP
static short *lowsieve = 0;
static short *movesieve = 0;
static long pindex = 0;
static long pshift = -1;
static long lastp = 0;

void mpz_pstart()
{
    long i;
    long j;
    long jstep;
    long jstart;
    long ibnd;
    short int *p;

    if (!lowsieve)
    {
        lowsieve = (short int *)calloc((size_t)PRIM_BND, sizeof(short int));
        p = &lowsieve[0];
        for (i = PRIM_BND; i; i--)
            *p++ = 1;
        jstep = 1;
        jstart = -1;
        ibnd = (mpz_sqrts((long) (2 * PRIM_BND + 1)) - 3) / 2;
        for (i = 0; i <= ibnd; i++)
        {
            jstart += 2 * ((jstep += 2) - 1);
            if (lowsieve[i])
                for (j = jstart; j < PRIM_BND; j += jstep)
                    lowsieve[j] = 0;
        }
    }
    lastp = 0;
    pshift = 0;
    pindex = -2;
}

void
mpz_pstart2()
{
    lastp = 0;
    pshift = -1;
}

static void mpz_pshift()
{
    /* auxiliary routine for prime generator */
    long i;
    long j;
    long jstep;
    long jstart;
    long ibound;
    short *p;
    short *pi;

    if (!movesieve)
        movesieve = (short int *)calloc((size_t)PRIM_BND, sizeof(short int));
    pi = &movesieve[0];
    if (!pshift)
    {
        p = &lowsieve[0];
        for (i = PRIM_BND; i; i--)
            *pi++ = *p++;
    }
    else
    {
        for (i = PRIM_BND; i; i--)
            *pi++ = 1;
        jstep = 3;
        ibound = pshift + 2 * PRIM_BND + 1;
        for (i = 0; jstep * jstep <= ibound; i++)
        {
            if (lowsieve[i])
            {
                if (!((jstart = (pshift + 2) / jstep + 1) & 1))
                    jstart++;
                if (jstart <= jstep)
                    jstart = jstep;
                jstart = (jstart * jstep - pshift - 3) / 2;
                for (j = jstart; j < PRIM_BND; j += jstep)
                    movesieve[j] = 0;
            }
            jstep += 2;
        }
    }
}

long mpz_pnext()
{
    if (pshift < 0)
    {
        mpz_pstart();
        return (lastp = 2);
    }
    if (pindex == -2)
    {
        pindex = 0;
        mpz_pshift();
        return (lastp = 3);
    }
    for (;;)
    {
        while ((++pindex) < PRIM_BND)
        {
            if (movesieve[pindex])
                return (lastp = pshift + 2 * pindex + 3);
        }
        if ((pshift += 2 * PRIM_BND) > 2 * PRIM_BND * (2 * PRIM_BND + 1))
        {
            mpz_pstart();
            return (lastp = 2);
        }
        mpz_pshift();
        pindex = -1;
    }
}

long mpz_p()
{
    return (lastp);
}

long mpz_pnextb(long b)
{
    if (b >= PRIM_UP)
        return(0);
    mpz_pstart();
    if (b < (2*PRIM_BND))
    {
        mpz_pstart2();
        while (mpz_pnext() < b);
        return(lastp);
    }
    pshift = ((b-2) / (2*PRIM_BND)) * (2*PRIM_BND);
    pindex = -1;
    mpz_pshift();
    while (mpz_pnext() < b);
    return(lastp);
}
// END Stolen from lip.c (c) A. Lenstra
// END Converted to GMP by C.T. Card 2002
#endif
// static class members
std::vector<long int> VeryLong::Prime_power;
std::vector<long int> VeryLong::Prime_power2;
long int VeryLong::Max_diff = 0;
std::vector<VeryLong> VeryLong::extraPrimes_;

VeryLong exp(const VeryLong& vl1, const VeryLong& vl2)
{
    VeryLong vl = pow<VeryLong, VeryLong>(vl1, vl2);
    return vl;
}
#ifdef USE_LIP
void VeryLong::generate_prime_powers()
{
    long B1 = 1000000;
    long B2 = 20000000;
    long p = 0;
    long p0 = 0;
    mpz_pstart2();

    //cout << "Generating prime powers up to B = " << B1 << " ..." << flush;
    p = mpz_pnext();
    while (p < B1)
    {

        p0 = p;
        long p1;
        while ((p1 = p*p0) > 0 && p1 < B1) p = p1;
        Prime_power.push_back(p);
        p = mpz_pnext();
    }
    //cout << " done (" << Prime_power.size() << ")" << endl;
    //cout << "Generating prime powers for step 2 up to B2 = " << B2 << " ..." << flush;
    long int pp = p0;
    long int diff = 0;
    while (p0 < B2)
    {
        Prime_power2.push_back(p0);
        pp = p0;
        p0 = mpz_pnext();
        diff = p0 - pp;
        if (diff > Max_diff) Max_diff = diff;
    }
    //cout << " done (" << Prime_power2.size() << ")" << endl;
}
#endif

bool VeryLong::factorise_p_minus_1(VeryLong* factor, VeryLong* new_N)
{
    const VeryLong one(1L);
    // Pollards p-1 algorithm, looking for factors of N = pq, when p-1 is smooth
    //cerr << "p-1 method" << endl;
    VeryLongModular x(*this, 2L);
    static int first_time = 1;
    if (first_time)
    {
        generate_prime_powers();
        first_time = 0;
    }

    int done = 0;
    long p;
    int i = 0;
    auto iter = Prime_power.begin();
    while (!done)
    {
        p = *iter;
        x.exp(p); // x = x^p mod N
        // test gcd(x-1, N)
        i++;
        if (i%100 == 0)
        {
            VeryLong y = x.get_very_long();
            y -= one;
            VeryLong g = gcd(y, *this);
            if (g != one && g != *this)
            {
                //cout << "Factor found : " << g << endl;
                *factor = g;
                *new_N = *this / g;
                //cout <<"leaving : " << *new_N << endl;
                return true;
            }
            //cout << "." << flush;
        }
        //if (i%8000 == 0) cout << endl;
        ++iter;
        if (iter == Prime_power.end()) done = 1;
    }
    // now for step 2
    //cerr << " ... step 2" << endl;
    long int diff = 0;
    VeryLongModular y;
    std::vector<VeryLongModular> yVector;
    for (diff = 2; diff <= Max_diff; diff += 2)
    {
        y = x;
        y.exp(diff);
        yVector.push_back(y);
    }

    auto iter2 = Prime_power2.begin();
    long int p0 = *iter2;
    long int p1 = 0;
    ++iter2;
    for (;
            iter2 != Prime_power2.end();
            ++iter2)
    {
        p1 = *iter2;
        diff = p1 - p0;
        y = yVector [ (diff - 2) / 2 ];
        x *= y;
        VeryLong xx = x.get_very_long();
        xx -= one;
        //VeryLong g = gcd(x.get_very_long() - VeryLong(1L), *this);
        VeryLong g = gcd(xx, *this);
        if (g != one && g != *this)
        {
            //cout << "Factor found in step 2: " << g << endl;
            *factor = g;
            *new_N = *this / g;
            //cout <<"leaving : " << *new_N << endl;
            return true;
        }
        p0 = p1;
    }

    return false;
}

// PrimeTable_ stores differences between primes
const int MAX_PRIME = 10000000;
const int PRIME_TABLE_SIZE = 1000000;
long int VeryLong::Max_prime_ = MAX_PRIME;
int VeryLong::PrimeTableSize_ = 0;
unsigned char* VeryLong::PrimeTable_ = 0;
unsigned char* VeryLong::PrimeTablePtr_ = VeryLong::PrimeTable_;
long int VeryLong::CurrentPrime_ = 2L;

void VeryLong::set_max_prime(long int B)
{
    if (B < MAX_PRIME) Max_prime_ = B;
    else Max_prime_ = MAX_PRIME;
}

void VeryLong::clear_prime_table()
{
    if (PrimeTable_)
        delete [] PrimeTable_;
    Max_prime_ = MAX_PRIME;
    PrimeTableSize_ = 0;
    PrimeTable_ = 0;
    PrimeTablePtr_ = VeryLong::PrimeTable_;
    CurrentPrime_ = 2L;
}

#ifdef USE_LIP
void VeryLong::generate_prime_table()
{
    if (PrimeTable_) return;
    PrimeTable_ = new unsigned char[PRIME_TABLE_SIZE];
    mpz_pstart2();
    long int p = mpz_pnext(); // 2
    long int prev_p = p;
    p = mpz_pnext(); // 3
    unsigned char* pt = PrimeTable_;
    while (p < MAX_PRIME)
    {
        int diff = p - prev_p;
        if (diff < 256)
        {
            *pt++ = static_cast<unsigned char>(p - prev_p);
            PrimeTableSize_++;
            prev_p = p;
            p = mpz_pnext();
        }
        else
        {
            cerr << "Prime difference too big : " << p << " - " << prev_p << " = " << diff << endl;
            break;
        }
    }
    PrimeTable_[PrimeTableSize_] = 3;
    //std::cerr << "PrimeTableSize_ = " << PrimeTableSize_ << std::endl;
}
#endif

static unsigned char equiv_class_mod30(long int p)
{
    unsigned char m = static_cast<unsigned char>(p % 30);
    static unsigned char lookup[] =
    {
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 2, 0, 3, 0, 0, 0, 4, 0, 5,
        0, 0, 0, 6, 0, 0, 0, 0, 0, 7
    };
    return lookup[m];
#if 0
    switch (m)
    {
    case 1:
        return 0;
    case 7:
        return 1;
    case 11:
        return 2;
    case 13:
        return 3;
    case 17:
        return 4;
    case 19:
        return 5;
    case 23:
        return 6;
    case 29:
        return 7;
    default:
        return 0;
    }
#endif
}
// Algorithm 8.1.1 (Trial Division)
bool VeryLong::factorise_trial_division(std::vector<VeryLong> * factors,
                                        VeryLong* new_N) const
{
    const VeryLong zero(0L);
    static unsigned char t[8] =
    {
        6, 4, 2, 4, 2, 4, 6, 2
    };
    const long int B = Max_prime_;
    generate_prime_table();

    *new_N = *this;
    // Step 1. [Initialise]
    int i = -1;
    int m = 0;
    VeryLong l = nth_root(2);

    long int d = 0;
    long int p = firstPrime();
    while (1)
    {
        if (i < 0)
        {
            // Step 2. [Next prime]
            m++;
            if (m > PrimeTableSize_ || d > B)
            {
                // Step 5. [Next divisor]
                i = equiv_class_mod30(p);
                d += t[i];
                if (d > B)
                {
                    for (auto& factor: extraPrimes_)
                    {
                        while (*new_N % factor == zero)
                        {
                            *new_N /= factor;
                            factors->push_back(factor);
                        }
                    }
                    return true;
                }
            }
            else
            {
                if (m > 1) p = nextPrime();
                d = p;
            }
        }
        else
        {
            // Step 5. [Next divisor]
            i++;
            if (i > 7) i = 0;
            d += t[i];
            if (d > B)
            {
                for (auto& factor: extraPrimes_)
                {
                    while (*new_N % factor == zero)
                    {
                        *new_N /= factor;
                        factors->push_back(factor);
                    }
                }
                return true;
            }
        }

        // Step 3. [Trial divide]
        VeryLong dd(d);
        const VeryLong zero(0L);
        int factor_found = 0;
        ;
        while (*new_N % dd == zero)
        {
            *new_N /= dd;
            factors->push_back(dd);
            factor_found = 1;
        }
        if (factor_found) l = new_N->nth_root(2);

        // Step 4. [Prime?]
        if (dd >= l) return true;
    }
}

bool VeryLong::factorise_trial_division(std::vector<long int> * factors,
                                        VeryLong* new_N) const
{
    static unsigned char t[8] =
    {
        6, 4, 2, 4, 2, 4, 6, 2
    };
    const long int B = Max_prime_;
    generate_prime_table();

    *new_N = *this;
    // Step 1. [Initialise]
    int i = -1;
    VeryLong l = nth_root(2);

    long int d = 0L;

    int m = 0;
    long int p = firstPrime();
    while (1)
    {
        if (i < 0)
        {
            // Step 2. [Next prime]
            m++;
            if (m > PrimeTableSize_ || d > B)
            {
                // Step 5. [Next divisor]
                i = equiv_class_mod30(p);
                d += t[i];
                if (d > B) return true;
            }
            else
            {
                if (m > 1) p = nextPrime();
                d = p;
            }
        }
        else
        {
            // Step 5. [Next divisor]
            i++;
            if (i > 7) i = 0;
            d += t[i];
            if (d > B) return true;
        }

        // Step 3. [Trial divide]
        int factor_found = 0;
        ;
        while (*new_N % d == 0L)
        {
            *new_N /= d;
            factors->push_back(d);
            factor_found = 1;
        }
        if (factor_found) l = new_N->nth_root(2);

        // Step 4. [Prime?]
        if (l <= d) return true;
    }
}

bool VeryLong::factorise_trial_division(std::vector<long int> * factors,
                                        VeryLong* new_N, long int B) const
{
    static long int t[8] =
    {
        6L, 4L, 2L, 4L, 2L, 4L, 6L, 2L
    };
    //const long int B = Max_prime_;
    generate_prime_table();
    //std::cout << "In factorise_trial_division, N = " << *this << ", B = " << B << std::endl;

    *new_N = *this;
    // Step 1. [Initialise]
    //std::cout << "Step 1" << std::endl;
    int i = -1;
    VeryLong l = nth_root(2);

    long int d = 0L;

    int m = 0;
    long int p = firstPrime();
    while (1)
    {
        if (i < 0)
        {
            // Step 2. [Next prime]
            //std::cout << "Step 2" << std::endl;
            m++;
            if (m > PrimeTableSize_ || d > B)
            {
                // Step 5. [Next divisor]
                //std::cout << "Step 5" << std::endl;

                i = equiv_class_mod30(p);
                d += t[i];
                if (d > B) return true;
            }
            else
            {
                if (m > 1) p = nextPrime();
                d = p;
            }
        }
        else
        {
            // Step 5. [Next divisor]
            //std::cout << "Step 5 (1)" << std::endl;
            i++;
            if (i > 7) i = 0;
            d += t[i];
            if (d > B) return true;
        }

        // Step 3. [Trial divide]
        //std::cout << "Step 3" << std::endl;
        int factor_found = 0;
        ;
        while (*new_N % d == 0L)
        {
            *new_N /= d;
            //std::cout << d << std::endl;
            factors->push_back(d);
            factor_found = 1;
        }
        if (factor_found) l = new_N->nth_root(2);

        // Step 4. [Prime?]
        //if (d % 5L == 0L) std::cout << "Step 4 : l = " << l << ", d = " << d << std::endl;
        if (l <= d) return true;
    }
}

bool VeryLong::factorise_trial_division(VeryLong* factor, VeryLong* new_N) const
{
    //cerr << "trial division method" << endl;
    long int B = Max_prime_;
    // try dividing by all odd numbers up to B
    long int trial_divisor = 2;
    if ((*this) % trial_divisor == 0L)
    {
        *factor = trial_divisor;
        *new_N = (*this) / trial_divisor;
        return true;
    }
    trial_divisor = 3;
    while (trial_divisor < B)
    {
        if ((*this) % trial_divisor == 0L)
        {
            *factor = trial_divisor;
            *new_N = (*this) / trial_divisor;
            return true;
        }
        trial_divisor += 2;
    }
    return false;
}

bool VeryLong::factorise_fermat(VeryLong* factor, VeryLong* new_N)
{
    const VeryLong one(1L);
    // find starting point near sqrt(N)
    //cerr << "Fermat method" << endl;
    VeryLong root;
    VeryLong diff;
    mpz_sqrtrem(root.vl_, diff.vl_, vl_);
    if (diff.is_zero())
    {
        // perfect square!
        *factor = root;
        *new_N = root;
        return true;
    }

    long B = 1000000;
    VeryLong r2 = *this;
    r2 -= diff;
    //cout << "r = " << root << endl;
    //cout << "d = " << diff << endl;
//	cout << "r^2 = " << r2 << endl;
    VeryLong g;
    long epsilon = 1;
    while (epsilon < B)
    {
        r2 -= epsilon;

        if ((g = gcd(r2, *this)) != one)
        {
            *factor = g;
            *new_N = (*this) / g;
            return true;
        }
        epsilon += 2;
        //if ((epsilon - 1)% 1000 == 0) cout << "." <<flush;
    }
    //cout << endl;

    return false;
}

double log10(const VeryLong& vl)
{
    long int power = static_cast<long int>(mpz_sizeinbase(vl.vl_, 2) - 1L);
    const VeryLong two(2L);
    mpz_t powerOf2;
    mpz_init(powerOf2);
    mpz_pow_ui(powerOf2, two.vl_, power);
    mpq_t q1;
    mpq_t q2;
    mpq_init(q1);
    mpq_init(q2);

    mpq_set_z(q1, vl.vl_);
    mpq_set_z(q2, powerOf2);
    mpq_div(q1, q1, q2);
    mpq_canonicalize(q1);

    double res = log10(mpq_get_d(q1)) + power * log10(2.0);

    mpz_clear(powerOf2);
    mpq_clear(q1);
    mpq_clear(q2);

    return res;
}

double ln(const VeryLong& vl)
{
    long int power = static_cast<long int>(mpz_sizeinbase(vl.vl_, 2) - 1L);
    const VeryLong two(2L);
    mpz_t powerOf2;
    mpz_init(powerOf2);
    mpz_pow_ui(powerOf2, two.vl_, power);
    mpq_t q1;
    mpq_t q2;
    mpq_init(q1);
    mpq_init(q2);

    mpq_set_z(q1, vl.vl_);
    mpq_set_z(q2, powerOf2);
    mpq_div(q1, q1, q2);
    mpq_canonicalize(q1);

    double res = log(mpq_get_d(q1)) + power * log(2.0);

    mpz_clear(powerOf2);
    mpq_clear(q1);
    mpq_clear(q2);

    return res;
}

VeryLong VeryLong::random(const VeryLong& a, const VeryLong& b)
{
    // return a random number in range [a,b]
    static int first_time = 1;
    static gmp_randstate_t state;
    if (first_time)
    {
        gmp_randinit_default(state);
        gmp_randseed_ui(state, static_cast<unsigned long int>(time(0)));
        first_time = 0;
    }
    VeryLong B = b;
    B -= a;
    if (b < a) B = -B;
    VeryLong r;
    mpz_urandomm(r.vl_, state, B.vl_);
    if (b < a) r += b;
    else r += a;
    return r;
}

#ifdef FASTVERYLONG
const VeryLong FastVeryLong::max_long_long_(LONG_LONG_MAX);
const long long FastVeryLong::max_long_(LONG_MAX);
#endif
