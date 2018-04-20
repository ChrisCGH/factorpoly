# assume Linux, use gcc, g++
GMPLIBDIR = /usr/local/lib
#COVERAGE_LIBS = -lgcov
#COVERAGE = -fprofile-arcs -ftest-coverage
#COVERAGE_OPT = -pg
COVERAGE_LIBS =
COVERAGE =
COVERAGE_OPT =
LIBS = -L/usr/lib64 -L/usr/local/lib $(COVERAGE_LIBS) -lgmp -lpthread -ltcmalloc_minimal -static-libstdc++
#LIBS = -L/usr/lib64 -L/usr/local/lib $(COVERAGE_LIBS) -lgmp -lpthread -static-libstdc++
INCLUDES = -I/usr/local/include
OPT = -O3
#CCFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH) -std=c++14 -fno-operator-names -Wno-non-template-friend -Wno-uninitialized -DUSING_GCC -Dlinux
CCFLAGS = -Wall -g $(COVERAGE_OPT) $(COVERAGE) $(OPT) $(PROFILE) $(ARCH) -std=c++14 
CFLAGS = -Wall -g $(COVERAGE_OPT) $(COVERAGE) $(OPT) $(PROFILE) $(ARCH)
CC = g++
cc = gcc
OBJEXT = .o
OUTPUTOPTION = -o

.SUFFIXES = $(OBJEXT) .cpp .c
.cpp$(OBJEXT): 
	$(CC) -c $(CCFLAGS) $<

.c$(OBJEXT):
	$(cc) -c $(CFLAGS) $<

%.d: %.cpp
	set -e; g++ -MM $(CPPFLAGS) $< \
       | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; [ -s $@ ] || rm -f $@

%.d: %.c
	set -e; gcc -MM $(CPPFLAGS) $< \
       | sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; [ -s $@ ] || rm -f $@

ALL_PROGRAMS : factorpoly

all : $(ALL_PROGRAMS)

ALL_HEADERS = crt.h discriminant.h gcd.h legendre.h lip.h lippar.h LongModular.h mod.h MPFloat.h mt19937int.h Polynomial.h pow.h Quotient.h VeryLong.h VeryLongModular.h Combinations.h Matrix.h lll.h timings.h
ALL_CPPS = factorpoly.cpp discriminant.cpp LongModular.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp lll.cpp timings.cpp

ALL_CS = lip.c mt19937int.c 
ALL_SRCS = $(ALL_CPPS) $(ALL_CS) $(ALL_HEADERS)
ALL_OBJS = $(ALL_CPPS:.cpp=$(OBJEXT)) $(C_SRCS:.c=$(OBJEXT))

FACTOR_POLYNO_CPP_SRCS = factorpoly.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp LongModular.cpp discriminant.cpp lll.cpp timings.cpp
C_SRCS = mt19937int.c lip.c

include $(FACTOR_POLYNO_CPP_SRCS:.cpp=.d)

FACTOR_POLYNO_OBJS = $(FACTOR_POLYNO_CPP_SRCS:.cpp=$(OBJEXT)) $(C_SRCS:.c=$(OBJEXT))

factorpoly : $(FACTOR_POLYNO_OBJS) 
	$(CC) $(CCFLAGS) $(OUTPUTOPTION) factorpoly $(FACTOR_POLYNO_OBJS) $(LIBS)

clean:
	rm -rf factorpoly *$(OBJEXT) *.d

