#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <cstdlib>
#include "Polynomial.h"
#include "VeryLong.h"
#include "Polynomial.inl"

std::string handler(const std::string& polynomial_string)
{
    Polynomial<VeryLong> poly = Polynomial<VeryLong>::read_polynomial(polynomial_string.c_str());
    std::vector<Polynomial<VeryLong> > factors;
    VeryLong cont;
    Polynomial<VeryLong>::factor(poly, factors, cont);
    //Polynomial<VeryLong>::factor_LLL(poly, factors, cont);
#if 0
    // Check that poly is fully factored
    Polynomial<VeryLong> remaining_part(poly);
    for (auto& factor: factors)
    {
        remaining_part /= factor;
    }
    if (remaining_part.deg() > 0)
    {
        factors.push_back(remaining_part);
    }

    // Check that factors can't be factored more
    bool done = false;
    while (!done)    
    {
        std::vector<Polynomial<VeryLong> > newfactors;
        for (auto& factor: factors)
        {
            std::vector<Polynomial<VeryLong> > subfactors;
            VeryLong subcont;
            Polynomial<VeryLong>::factor(factor, subfactors, subcont);
            newfactors.insert(std::end(newfactors), std::begin(subfactors), std::end(subfactors));
        }
        if (newfactors.size() == factors.size())
        {
            done = true;
        }
        else
        {
            factors = newfactors;
        }
    }
#endif
    // Gather common factors
    std::map<std::string, int> gathered_factors;
    for (auto& factor: factors)
    {
        std::ostringstream oss;
        oss << factor;
        gathered_factors[oss.str()]++;
    }
    std::ostringstream oss;
    if (cont == VeryLong(-1L))
    {
        oss << "-";
    }
    else if (cont != VeryLong(1L))
    {
        oss << cont;
    }
    for (auto& gf: gathered_factors)
    {
        std::string f = gf.first;
        int n = gf.second;
        if (n == 1)
        {
            oss << "(" << f << ")";
        }
        else 
        {
            if (n < 10)
            {
                oss << "(" << f << ")^" << n; 
            }
            else
            {
                oss << "(" << f << ")^" << "{" << n << "}"; 
            }
        }
    }
//    for (auto& factor: factors)
//    {
//        oss << "(" << factor << ")";
//    }
    return oss.str();
}

std::string mathjax_output(const std::string& p)
{
    std::ostringstream oss;
    oss << "<!DOCTYPE html>" << std::endl;
    oss << "<html>" << std::endl;
    oss << "<head>" << std::endl;
    oss << "  <meta charset=\"utf-8\">" << std::endl;
    oss << "  <meta name=\"viewport\" content=\"width=device-width\">" << std::endl;
    oss << "  <title>factorpoly output</title>" << std::endl;
    oss << "  <script type=\"text/javascript\" async" << std::endl;
    oss << "  src=\"https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML\">" << std::endl;
    oss << "</script>" << std::endl;
    oss << "</head>" << std::endl;
    oss << "<body>" << std::endl;
    oss << "<p>" << std::endl;
    oss << "$$" << p << "$$" << std::endl;
    oss << "</p>" << std::endl;
    oss << "</body>" << std::endl;
    oss << "</html>" << std::endl;
    return oss.str();
}

int main(int argc, char* argv[])
{
    std::string input;
    std::getline(std::cin, input);
    if (input.find("polynomial=") == 0)
    {
        input.replace(0, 11, "");
    }
    input.erase(std::remove(input.begin(), input.end(), '\n'), input.end());
    input.erase(std::remove(input.begin(), input.end(), '\r'), input.end());
    input.erase(std::remove(input.begin(), input.end(), '\n'), input.end());
    input.erase(std::remove(input.begin(), input.end(), '\r'), input.end());
    try 
    {
        std::string output = handler(input);
        if (std::getenv("USE_MATHJAX"))
        {
            output = mathjax_output(output);
        }
        std::cout << output << std::endl;
    }
    catch (const std::string& e)
    {
        std::cout << e << std::endl;
    }

    return 0;
}

