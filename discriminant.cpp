#include <iostream>
#include "Polynomial.h"
#include "Quotient.h"
#include "VeryLong.h"
#include <vector>
#include <iostream>
#include "pow.h"

//#define LIDIA_RESULTANT 1
#ifdef LIDIA_RESULTANT
// Code taken from LiDIA
// The following two functions were contributed
//  by Roland Dreier (dreier@math.berkeley.edu)

// Calculate resultant of aa and bb via subresultant algorithm.
// Algorithm cribbed verbatim from Algorithm 3.3.7 of H. Cohen's "A
// course in computational algebraic number theory."  (Even the
// variables are named pretty much the same as in his book!)
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//

VeryLong
resultant(const Polynomial <VeryLong> &aa,
          const Polynomial <VeryLong> &bb)
{
    bool debug(false);
    if (std::getenv("FACTOR_VERBOSE_OUTPUT")) debug = true;
    const VeryLong zero(0L);
    const VeryLong one(1L);
    // Return zero if either polynomial is zero.
    if ((aa == Polynomial<VeryLong>(0L)) || (bb == Polynomial<VeryLong>(0L)))
    {
        return zero;
    }
    // otherwise...

    // Initialization steps:

    VeryLong acont, bcont;
    Polynomial <VeryLong> a, b;
    bool neg = false;

    // Maybe skip one reduction by making sure deg(a) >= deg(b).
    if (aa.deg() >= bb.deg())
    {
        a = aa;
        b = bb;
    }
    else
    {
        a = bb;
        b = aa;

        // Possibly change sign!!
        neg = ((a.deg() % 2) && (b.deg() % 2));
    }
    if (debug) std::cout << "::::: resultant of " << a << " and " << b << std::endl;

    acont = a.content();
    bcont = b.content();

    a.make_primitive();
    b.make_primitive();

    VeryLong g = 1L;
    VeryLong h = 1L;

    VeryLong t;
    VeryLong pow_temp;

    //t = exp(VeryLong(acont), VeryLong((long int)b.deg()));
    t = pow<VeryLong, long int>(acont, (long int)b.deg());
    //pow_temp = exp(VeryLong(bcont), VeryLong((long int)a.deg()));
    pow_temp = pow<VeryLong, long int>(bcont, (long int)a.deg());
    t = t * pow_temp;

    // Loop over pseudo-division and reduction steps:

    long int delta;
    Polynomial <VeryLong> r;

    do
    {
        delta = a.deg() - b.deg();

        if ((a.deg() % 2) && (b.deg() % 2))
        {
            neg = !neg;
        }

        if (debug) std::cout << "::::: about to calculate remainder(" << a << "," << b << ")" << std::endl;

        r = remainder(a, b);
        a = b;

        //pow_temp = exp(h, VeryLong(delta));
        pow_temp = pow<VeryLong, long int>(h, delta);
        if (debug) std::cout << "::::: 1. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << std::endl;
        b = r / pow_temp;
        b = b / g;

        if (debug) std::cout << "::::: b = " << b << std::endl;

        if (debug) std::cout << "::::: Top coeff of a = " << a.coefficient(a.deg()) << std::endl;
        g = a.coefficient(a.deg());
        //pow_temp = exp(g, VeryLong(delta--));
        pow_temp = pow<VeryLong, long int>(g, delta--);
        if (debug) std::cout << "::::: 2. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << std::endl;

        if (delta<=0)
        {
            //h = exp(h, VeryLong(-delta));
            h = pow<VeryLong, long int>(h, -delta);
            h = h * pow_temp;

        }
        else
        {
            //h = exp(h, VeryLong(delta));
            h = pow<VeryLong, long int>(h, delta);
            h = pow_temp / h;
        }
        if (debug) std::cout << "::::: 3. r = " << r << ", pow_temp = " << pow_temp << ", g = " << g << ", h = " << h << std::endl;
    }
    while (b.deg() > 0);

    // Finish up:

    VeryLong temp(0L);
    if (b.deg() > 0)
    {
       temp = b.coefficient(b.deg());
    }
    if (debug) std::cout << "::::: temp = " << temp << std::endl;
    //pow_temp = exp(temp, VeryLong((long int)a.deg()));
    pow_temp = pow<VeryLong, long int>(temp, (long int)a.deg());

    delta = a.deg()-1;
    if (delta <= 0)
    {
        //h = exp(h, VeryLong(-delta));
        h = pow<VeryLong, long int>(h, -delta);
        h = h * pow_temp;
    }
    else
    {
        //h = exp(h,VeryLong(delta));
        h = pow<VeryLong, long int>(h,delta);
        h = pow_temp / h;
    }
    if (neg)
        t = -one * t;
    VeryLong res = t*h;
    if (debug) std::cout << "::::: Resultant = " << res << std::endl;
    return res;
}
#else
// Algorithm 3.3.7 [Sub-Resultant]
VeryLong resultant(const Polynomial <VeryLong>& AA, const Polynomial <VeryLong>& BB)
{
    bool debug(false);
    if (std::getenv("FACTOR_VERBOSE_OUTPUT")) debug = true;
    const VeryLong zero(0L);
    const Polynomial<VeryLong> zeropoly(0L);
    const VeryLong one(1L);
    auto A(AA);
    auto B(BB);
    // Step 1. [Initializations and reductions]
    if (A == zeropoly || B == zeropoly)
    {
        return zero;
    }
    VeryLong a = A.content();
    VeryLong b = B.content();
    A /= a;
    B /= b;
    VeryLong g(1L);
    VeryLong h(1L);
    long int s(1L);
    auto t = pow<VeryLong, int>(a, B.deg()) * pow<VeryLong, int>(b, A.deg());
    if (debug) std::cout << "::::: A = " << A << ", B = " << B << ", g = " << g << ", h = " << h << ", s = " << s << std::endl;

    if (A.deg() < B.deg())
    {
        auto tmp = A;
        A = B;
        B = tmp;
        if (A.deg() % 2 && B.deg() % 2)
        {
            s = -1L;
        }
    }
    if (debug) std::cout << "::::: A = " << A << ", B = " << B << ", g = " << g << ", h = " << h << ", s = " << s << std::endl;

    while (B.deg() > 0)
    {
        // Step 2. [Pseudo division]
        auto delta = A.deg() - B.deg();
        if (A.deg() % 2 && B.deg() % 2)
        {
            s = -s;
        }
        Polynomial<VeryLong> Q;
        Polynomial<VeryLong> R;
        pseudo_divide(A, B, Q, R);
        if (debug) std::cout << "::::: Step 2. A = " << A << ", B = " << B << ", Q = " << Q << ", R = " << R << std::endl;

        // Step 3. [Reduce remainder]
        A = B;
        if (delta < zero)
        {
            B = R * pow<VeryLong, int>(h, -delta) / g;
        }
        else
        {
            B = R / (g * pow<VeryLong, int>(h, delta));
        }
        if (debug) std::cout << "::::: Step 3. A = " << A << ", B = " << B << std::endl;

        // Step 4. [Finished]
        g = A.coefficient(A.deg());
        if (delta <= 0)
        {
            h = pow<VeryLong, int>(h, 1 - delta) / pow<VeryLong, int>(g, -delta);
        }
        else 
        {
            h = pow<VeryLong, int>(h, 1 - delta) * pow<VeryLong, int>(g, delta);
        }
        if (debug) std::cout << "::::: Step 4. A = " << A << ", B = " << B << ", g = " << g << ", h = " << h << ", s = " << s << std::endl;
    }
    if (debug) std::cout << "::::: A = " << A << ", B = " << B << ", g = " << g << ", h = " << h << ", s = " << s << std::endl;
    auto l_B = B.coefficient(B.deg());
    if (debug) std::cout << "::::: l_B = " << l_B << std::endl;
    if (A.deg() <= 1)
    {
        h = pow<VeryLong, int>(h, 1 - A.deg()) * pow<VeryLong, int>(l_B, A.deg());
    }
    else
    {
        h = pow<VeryLong, int>(l_B, A.deg()) / pow<VeryLong, int>(h, A.deg() - 1);
    }
    if (debug) std::cout << "::::: h = " << h << ", s = " << s << ", t = " << t << ", s*t*h = " << s * t * h << std::endl;
    return s * t * h;
}
#endif

//
// Calculate discriminant of a polynomial using the formula:
//   disc(a) = (-1)^(d*(d-1)/2) * resultant(a, a') / lead_coeff(a)
// where d = deg(f).
//
// Rather than raising -1 to a power, just look at d mod 4.
// Remember, d*(d-1)/2 is even iff d = 0 or 1 mod 4.
//
// Author      : Roland Dreier (dreier@math.berkeley.edu)
//
VeryLong
discriminant(const Polynomial <VeryLong> &a)
{
//	cout << "discriminant of " << a << endl;
    if (a.deg() <= 0)
    {
        return(VeryLong(0L));
    }

    //VeryLong multiplier = exp(a.coefficient(a.deg()), VeryLong(a.deg()*2L - 3L));
    VeryLong multiplier = pow<VeryLong, long int>(a.coefficient(a.deg()), a.deg()*2L - 3L);
    if ((int(a.deg()) % 4) <= 1)
    {
        return(resultant(a, a.derivative()) / a.coefficient(a.deg()));
    }
    else
    {
        return(-(resultant(a, a.derivative()) / a.coefficient(a.deg())));
    }
}
