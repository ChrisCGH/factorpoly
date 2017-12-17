# assume Linux, use gcc, g++
GMPLIBDIR = /usr/local/lib
LIBS = -L/usr/lib64 -L/usr/local/lib -lgmp -static-libstdc++
INCLUDES = -I/usr/local/include
OPT = -O3
#CCFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH) -std=c++14 -fno-operator-names -Wno-non-template-friend -Wno-uninitialized -DUSING_GCC -Dlinux
CCFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH) -std=c++14 
CFLAGS = -Wall -g $(OPT) $(PROFILE) $(ARCH)
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

ALL_HEADERS = crt.h gcd.h legendre.h lip.h lippar.h LongModular.h mod.h MPFloat.h mt19937int.h Polynomial.h pow.h Quotient.h VeryLong.h VeryLongModular.h
ALL_CPPS = factorpoly.cpp LongModular.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp

ALL_CS = lip.c mt19937int.c 
ALL_SRCS = $(ALL_CPPS) $(ALL_CS) $(ALL_HEADERS)
ALL_OBJS = $(ALL_CPPS:.cpp=$(OBJEXT)) $(C_SRCS:.c=$(OBJEXT))

FACTOR_POLYNO_CPP_SRCS = factorpoly.cpp Polynomial.cpp VeryLong.cpp VeryLongModular.cpp LongModular.cpp
C_SRCS = mt19937int.c lip.c

include $(FACTOR_POLYNO_CPP_SRCS:.cpp=.d)

FACTOR_POLYNO_OBJS = $(FACTOR_POLYNO_CPP_SRCS:.cpp=$(OBJEXT)) $(C_SRCS:.c=$(OBJEXT))

factorpoly : $(FACTOR_POLYNO_OBJS) 
	$(CC) $(CCFLAGS) $(OUTPUTOPTION) factorpoly $(FACTOR_POLYNO_OBJS) $(LIBS)

clean:
	rm -rf factorpoly *$(OBJEXT) *.d

