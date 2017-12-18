# factorpoly

factorpoly is a standalone single-variate polynomial factorizer, written in C++, based on code taken from [factor-by-gnfs](https://sourceforge.net/projects/factor-by-gnfs/)

# Building
## Standalone executable
Clone the code from Github:
`git clone https://github.com/ChrisCGH/factorpoly.git`
Change directory to factorpoly and build using make:
`cd factorpoly`
`make`
This should create the executable `factorpoly` in the same directory.
### Dependencies
The supplied makefile is written for Linux, and assumes that the standard toolchain for building C++ programs is already installed (e.g. make, gcc, g++). In addition the code makes use of [the GNU MP Bignum library](https://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjDrOfd7JPYAhVUF8AKHQGLBWMQFggnMAA&url=https%3A%2F%2Fgmplib.org%2F&usg=AOvVaw0J-ZxDUBVNbeR6FiUDVsnH), and assumes that the development library for GMP is installed.
Building on other systems is left as an exercise, but the code compiles with g++ and -std=c++14, so should compile with a recent C++ compiler.

## Docker image suitable for OpenFaas
factorpoly can be built as a Docker image suitable for use with [OpenFaas](https://www.openfaas.com/) using the supplied Dockerfile.
The latest image can be found [here](https://hub.docker.com/r/ccard/factorpoly/).

# Running
## Standalone executable
To run the standalone factorpoly simply read a polynomial from stdin:

`$ echo "X^55 + 1" | ./factorpoly`
`(X + 1)(X^10 - X^9 + X^8 - X^7 + X^6 - X^5 + X^4 - X^3 + X^2 - X + 1)(X^4 - X^3 + X^2 - X + 1)(X^40 + X^39 - X^35 - X^34 + X^30 - X^28 - X^25 + X^23 + X^20 + X^17 - X^15 - X^12 + X^10 - X^6 - X^5 + X + 1)`

If the environment variable USE_MATHJAX is set, then the output is an html page which uses mathjax to format the polynomial factors:
```
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width">
  <title>factorpoly output</title>
  <script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.2/MathJax.js?config=TeX-MML-AM_CHTML">
</script>
</head>
<body>
<p>
$$(X + 1)(X^4 - X^3 + X^2 - X + 1)(X^{10} - X^9 + X^8 - X^7 + X^6 - X^5 + X^4 - X^3 + X^2 - X + 1)(X^{40} + X^{39} - X^{35} - X^{34} + X^{30} - X^{28} - X^{25} + X^{23} + X^{20} + X^{17} - X^{15} - X^{12} + X^{10} - X^6 - X^5 + X + 1)$$
</p>
</body>
</html>
```

## OpenFaas function
factorpoly can be run as an OpenFaas function, by deploying the function from the image ccard/factorpoly.
Once deployed, the function can be run e.g. by using curl:
`curl http://localhost:8080/function/factorpoly -X POST --data-binary 'X^55 + 1'`

# References
The algorithm used by factorpoly is an implementation of Algorithm 3.5.7 in "A Course in Computational Algebraic Number Theory" by Henri Cohen (http://www.springer.com/gp/book/9783540556404)
