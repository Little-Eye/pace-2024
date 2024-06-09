/* int128.h - 128-bit integer arithmetic for C++, by Robert Munafo

   (and copied to rhtf/.../int128.h.txt by proj/.../MCS.per)

LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


The license (GPL version 3) for this source code is at:
mrob.com/pub/ries/GPL-3.txt


HOW TO USE THIS OBJECT:

 Single-file projects:
   *- Main.c -*
   #include "int128.c"
   s128 x;
   main() {
     x = 3; x = x * 3;
   }

 Complex projects:
   *- makefile -*
   (Note: the function mk_compile is in rpmlib.pl)
   &mk_compile("/home/munafo/shared/proj/include/int128.c");
   # Compile with flag -DI128_COMPAT_PLUS_EQUAL if on a modern compiler
   &mk_compile("other_module.cxx");
   &mk_compile("math3000.cxx", "int128.o other_module.o");


REVISION HISTORY is in int128.c

*/

/* Compile line should supply -I/Users/munafo/shared/proj/include
   or equivalent */
#include "rpmtypes.h"

#define BITS62_0 0x7FFFFFFFFFFFFFFFLL
#define BIT_52   0x0010000000000000LL
#define BITS51_0 0x000FFFFFFFFFFFFFLL
#define HIGHBIT  0x8000000000000000LL
#define ALLBITS  0xFFFFFFFFFFFFFFFFLL

#define SIGBITS_HI 0x7FFFFFFF00000000LL
#define HI_WORD 0xFFFFFFFF00000000LL
#define LO_WORD 0x00000000FFFFFFFFLL


class s128_o
{
 private:

 public:
  s64 high;
  u64 low;

  s128_o ( );
  s128_o ( int x );
  s128_o ( u64 x );
  s128_o ( s64 x );
  s128_o ( double x );
  s128_o ( const s128_o & );

  operator int ( )
  {
    int result;

    result = low & LO_WORD;
    return result;
  }

  operator unsigned long long ( )
  {
    long long result;

    result = (low & BITS62_0) | (high & HIGHBIT);
    return result;
  }

  operator long long ( )
  {
    long long result;

    return low;
  }

  operator double ( )
  {
    double result;

    if (high >= 0) {
      result = (((double) high) * ((double) 18446744073709551616.0)) + ((double) low);
    } else {
      s64 h; u64 l;
      h = high; l = low;

      h = ~h;
      l = ~l;
      l += 1;
      if (l == 0) {
        h += 1;
      }
      result = -((((double) h) * ((double) 18446744073709551616.0)) + ((double) l));
    }

    return result;
  }

  friend s128_o operator + ( const s128_o &, const s128_o & );
  friend s128_o operator + ( s64, const s128_o & );
  friend s128_o operator + ( const s128_o &, s64 );
  friend s128_o operator + ( const s128_o &, u64 );
  friend s128_o operator + ( const s128_o &, s32 );

  friend s128_o operator - ( const s128_o &, const s128_o & );
  friend s128_o operator - ( const s128_o &, s32 );

  friend s128_o operator - ( const s128_o & );

  friend s128_o s128_shr(s128_o);
  friend s128_o s128_shl(s128_o);

  friend s128_o mult1(s128_o x, s128_o y);
  friend s128_o operator * ( const s128_o &, const s128_o & );
  friend s128_o operator * ( const s128_o &, u64 );
  friend s128_o operator * ( const s128_o &, s32 );

  friend s128_o div1(s128_o, s128_o, s128_o *);
  friend s128_o operator / ( const s128_o &, const s128_o & );
  friend s128_o operator / ( const s128_o &, u64 );

  friend s128_o operator += ( const s128_o &, const s128_o & );
  friend s128_o operator += ( const s128_o &, const s64 & );
  friend s128_o operator -= ( const s128_o &, const s128_o & );
  friend s128_o operator -= ( const s128_o &, const s64 & );
  friend s128_o operator *= ( const s128_o &, const s128_o & );
  friend s128_o operator *= ( const s128_o &, const s64 & );
  friend s128_o operator /= ( const s128_o &, const s128_o & );
  friend s128_o operator /= ( const s128_o &, const s64 & );

  friend int operator < ( const s128_o &, const s128_o & );
  friend int operator < ( const s128_o &, s32 );

  friend int operator <= ( const s128_o &, const s128_o & );
  friend int operator <= ( const s128_o &, s32 );

  friend int operator == ( const s128_o &, const s128_o & );
  friend int operator == ( const s128_o &, s32 );

  friend int operator != ( const s128_o &, const s128_o & );
  friend int operator != ( const s128_o &, s32 );

  friend int operator > ( const s128_o &, const s128_o & );
  friend int operator > ( const s128_o &, s32 );

  friend int operator >= ( const s128_o &, const s128_o & );
};

typedef s128_o s128;

extern u64 s128_u64(s128 x);
extern u32 s128_u32(s128 x);
extern void s128_str(s128 x, char *s);

inline s128_o::s128_o ( ): high(0), low(0) { }

inline s128_o::s128_o ( u64 x ): high(0), low(x) { }

inline s128_o::s128_o ( const s128_o & x ): high(x.high), low(x.low) { }

/* When I originally made this library (sometime before 2004) the
compiler I was using at the time was probably GCC on the PowerPC G4
iMac. For inlined += operator, the following syntax worked:

  inline s128_o::s128_o operator += ( s128_o & lhs, const s128_o & rhs )
  {
    lhs = lhs + rhs;
  }

Sometime since then, this syntax became illegal, and generates the
useless error "'s128_o::s128_o' names the constructor, not the type".
The simplest way to fix it appears to be to just change "s128_o::s128_o"
to "s128_o", thus:

  inline s128_o operator += ( s128_o & lhs, const s128_o & rhs )
  {
    lhs = lhs + rhs;
  }

For more, see stackoverflow.com/questions/13074590

*/

#ifdef I128_COMPAT_PLUS_EQUAL
# define S128O_S128O s128_o
#else
# ifdef I128_COMPAT_PEQ_AMP
#   define S128O_S128O s128_o&
# else
#   define S128O_S128O s128_o::s128_o
# endif
#endif


inline S128O_S128O operator += ( s128_o & lhs, const s128_o & rhs ) 
{
  lhs = lhs + rhs;
  return lhs;
}

inline S128O_S128O operator += ( s128_o & lhs, const s64 & rhs ) 
{
  lhs = lhs + rhs;
  return lhs;
}

inline S128O_S128O operator -= ( s128_o & lhs, const s128_o & rhs ) 
{
  lhs = lhs - rhs;
  return lhs;
}

inline S128O_S128O operator -= ( s128_o & lhs, const s64 & rhs ) 
{
  lhs = lhs - ((s128_o)rhs);
  return lhs;
}

inline S128O_S128O operator *= ( s128_o & lhs, const s128_o & rhs ) 
{
  lhs = lhs * rhs;
  return lhs;
}

inline S128O_S128O operator *= ( s128_o & lhs, const s64 & rhs ) 
{
  lhs = lhs * ((s128_o) rhs);
  return lhs;
}

inline S128O_S128O operator /= ( s128_o & lhs, const s128_o & rhs ) 
{
  lhs = lhs / rhs;
  return lhs;
}

inline S128O_S128O operator /= ( s128_o & lhs, const s64 & rhs ) 
{
  lhs = lhs / ((s128_o) rhs);
  return lhs;
}


/* end of int128.h */
