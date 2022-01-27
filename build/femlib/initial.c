/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 1 "initial.y"

/*****************************************************************************
 * Synopsis
 * --------
 * initial.y: yacc code for a simple function interpreter, which also allows
 * lookup for named tokens (all are double precision).
 *
 * Modelled on "hoc3" in Chapter 8 of "The UNIX Programming Environment"
 * by Kernighan & Pike, Prentice-Hall 1984.  A hash-table data structure is
 * used in place of the linked list of Symbols which they employed.
 * Hashing code from Chapter 6 of "The C Programming Language", 2nd Edn,
 * by Kernighan & Ritchie, Prentice-Hall 1988.
 *
 * Copyright (c) 1994 <--> $Date$, 
 *   Ron Henderson, Hugh Blackburn, Spencer Sherwin
 *
 * Summary
 * -------
 * void    yy_initialize (void);
 * void    yy_help       (void);
 * void    yy_show       (void);
 * int_t   yy_dump       (char*, const int_t);
 *
 * double  yy_interpret  (const char*);
 *
 * void    yy_vec_init   (const char*, const char*);
 * void    yy_vec_interp (const int_t, ...);
 *
 * Notes
 * -----
 * 1. yy_initialize must be called before other routines will work.
 * 2. yy_help prints a summary of available functions on stdout.
 * 3. yy_show prints a summary of installed variables on stdout.
 * 4. yy_dump is similar to yy_show, but prints into a string of given length.
 * 5. yy_interpret is the central routine.  It can be used:
 *    (a), to install a named symbol in internal tables: e.g. "name = value";
 *    (b), to retrieve the value of an installed symbol: e.g. "name";
 *    (c), to evaluate a function string: e.g. "cos(x)*exp(-t)".
 * 6. yy_vec_init is used to set up the interpreter for vector evaluation.
 * 7. yy_vec_interp subsequently used for "vectorized" calls to yy_interpret.
 *
 * Operators
 * ---------
 * Unary:      -
 * Binary:     -, +, *, /, ^ (exponentiation), ~ (atan2), & (hypot), % (fmod)
 * Functions:  sin,  cos,  tan,  abs, floor, ceil, int, heav (Heaviside),
 *             asin, acos, atan, log, log10, exp,  sqrt,
 *             sinh, cosh, tanh, asinh, acosh, atanh,
 *             erf, erfc, gamma, lgamma, sgn,
 *             j0, j1, y0, y1, jn.
 * Procedures: jn, yn, rad, ang, rejn, imjn, jacobi, womcos, womsin.
 *
 * --
 * This file is part of Semtex.
 * 
 * Semtex is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * Semtex is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Semtex (see the file COPYING); if not, write to the Free
 * Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define __USE_MISC
#include <math.h>

#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>

#include <cfemdef.h>
#include <cveclib.h>

#if 1
#define HASHSIZE 199
#else
#define HASHSIZE 37
#endif
#define HASHSEED 31
#define VEC_MAX  32

typedef double (*PFD)( ); /* NB: no arguments -- non-ANSI (on purpose). */

typedef struct symbol {
  char* name;
  short type;
  union {
    double val;			/* -- If VAR.   */
    PFD    ptr;			/* -- If BLTIN. */
  } u;
  struct symbol* next;
} Symbol;

static double
  Sgn    (double),
  Heavi  (double), 
  White  (double),
  Gamma  (double),
  Step   (double,double),
  Jn     (double,double),
  Yn     (double,double),
  Rad    (double,double),
  Ang    (double,double),
  ReJn   (double,double,double),
  ImJn   (double,double,double),
  Jacobi (double,double,double,double),
  Womcos (double,double,double,double,double),
  Womsin (double,double,double,double,double);

static unsigned hash     (const char*);
static Symbol*  lookup   (const char*);
static Symbol*  install  (const char*, const int_t, const double);
static void*    emalloc  (const size_t);

       int      yyparse (void);
static int      yylex   (void);
static void     yyerror (char*);

static double  value;
static Symbol* hashtab[HASHSIZE];
static char    func_string[STR_MAX], *cur_string;
static int_t   nvec = 0;
static Symbol* vs[VEC_MAX];
extern int     errno;

static struct {			    /* -- Built-in functions. */
  char* name;
  short narg;
  PFD   func;
} builtin[] = {
  "cos"   ,  1, cos     ,
  "sin"   ,  1, sin     ,
  "tan"   ,  1, tan     ,
  "exp"   ,  1, exp     ,
  "sinh"  ,  1, sinh    ,
  "cosh"  ,  1, cosh    ,
  "tanh"  ,  1, tanh    ,
  "erf"   ,  1, erf     ,
  "erfc"  ,  1, erfc    ,
  "int"   ,  1, rint    ,
  "abs"   ,  1, fabs    ,
  "floor" ,  1, floor   ,
  "ceil"  ,  1, ceil    ,
  "acos"  ,  1, acos    ,
  "asin"  ,  1, asin    ,
  "atan"  ,  1, atan    ,
  "acosh" ,  1, acosh   ,
  "asinh" ,  1, asinh   ,
  "atanh" ,  1, atanh   ,
  "log"   ,  1, log     ,
  "log10" ,  1, log10   ,
  "sqrt"  ,  1, sqrt    ,
  "heav"  ,  1, Heavi   ,
  "j0"    ,  1, j0      ,
  "j1"    ,  1, j1      ,
  "y0"    ,  1, y0      ,
  "y1"    ,  1, y1      ,
  "gamma" ,  1, Gamma   ,
  "lgamma",  1, lgamma  ,
  "white",   1, White   ,
  "sgn",     1, Sgn     ,

  "jn"    ,  2, Jn      ,
  "yn"    ,  2, Yn      ,
  "rad"   ,  2, Rad     ,
  "ang"   ,  2, Ang     ,
  "step"  ,  2, Step    ,

  "rejn"  ,  3, ReJn    ,
  "imjn"  ,  3, ImJn    ,
  
  "jacobi",  4, Jacobi  ,

  "womcos",  5, Womcos  ,
  "womsin",  5, Womsin  ,

  NULL, 0, NULL
};

#include "defaults.h"


#line 266 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

#include "initial.h"
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_NUMBER = 3,                     /* NUMBER  */
  YYSYMBOL_VAR = 4,                        /* VAR  */
  YYSYMBOL_BLTIN_UNARY = 5,                /* BLTIN_UNARY  */
  YYSYMBOL_BLTIN_BINARY = 6,               /* BLTIN_BINARY  */
  YYSYMBOL_BLTIN_TERNARY = 7,              /* BLTIN_TERNARY  */
  YYSYMBOL_BLTIN_QUATERNARY = 8,           /* BLTIN_QUATERNARY  */
  YYSYMBOL_BLTIN_QUINTERNARY = 9,          /* BLTIN_QUINTERNARY  */
  YYSYMBOL_UNDEF = 10,                     /* UNDEF  */
  YYSYMBOL_11_ = 11,                       /* '='  */
  YYSYMBOL_12_ = 12,                       /* '+'  */
  YYSYMBOL_13_ = 13,                       /* '-'  */
  YYSYMBOL_14_ = 14,                       /* '*'  */
  YYSYMBOL_15_ = 15,                       /* '/'  */
  YYSYMBOL_UNARYMINUS = 16,                /* UNARYMINUS  */
  YYSYMBOL_17_ = 17,                       /* '^'  */
  YYSYMBOL_18_ = 18,                       /* '~'  */
  YYSYMBOL_19_ = 19,                       /* '&'  */
  YYSYMBOL_20_n_ = 20,                     /* '\n'  */
  YYSYMBOL_21_ = 21,                       /* '('  */
  YYSYMBOL_22_ = 22,                       /* ')'  */
  YYSYMBOL_23_ = 23,                       /* ','  */
  YYSYMBOL_24_ = 24,                       /* '%'  */
  YYSYMBOL_YYACCEPT = 25,                  /* $accept  */
  YYSYMBOL_list = 26,                      /* list  */
  YYSYMBOL_asgn = 27,                      /* asgn  */
  YYSYMBOL_expr = 28                       /* expr  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int8 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   317

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  25
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  4
/* YYNRULES -- Number of rules.  */
#define YYNRULES  24
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  74

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   266


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
      20,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,    24,    19,     2,
      21,    22,    14,    12,    23,    13,     2,    15,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    11,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,    17,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,    18,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    16
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_uint8 yyrline[] =
{
       0,   209,   209,   210,   211,   212,   214,   216,   217,   223,
     224,   226,   228,   230,   232,   234,   235,   236,   237,   243,
     244,   245,   246,   247,   248
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "NUMBER", "VAR",
  "BLTIN_UNARY", "BLTIN_BINARY", "BLTIN_TERNARY", "BLTIN_QUATERNARY",
  "BLTIN_QUINTERNARY", "UNDEF", "'='", "'+'", "'-'", "'*'", "'/'",
  "UNARYMINUS", "'^'", "'~'", "'&'", "'\\n'", "'('", "')'", "','", "'%'",
  "$accept", "list", "asgn", "expr", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-21)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-1)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
     -21,    24,   -21,   -21,    -8,   -20,   -17,    -9,    -7,     1,
      59,   -21,    59,     5,    64,    59,    59,    59,    59,    59,
      59,   -21,    29,    77,   -21,    59,    59,    59,    59,    59,
      59,    59,   -21,    59,   285,    90,   103,   116,   129,   142,
     -21,   293,   293,    29,    29,    29,    29,    29,   285,   -21,
      59,    59,    59,    59,   155,   168,   181,   194,   -21,    59,
      59,    59,   207,   220,   233,   -21,    59,    59,   246,   259,
     -21,    59,   272,   -21
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int8 yydefact[] =
{
       2,     0,     1,     7,     8,     0,     0,     0,     0,     0,
       0,     3,     0,     9,     0,     0,     0,     0,     0,     0,
       0,     9,    24,     0,     4,     0,     0,     0,     0,     0,
       0,     0,     5,     0,     6,     0,     0,     0,     0,     0,
      23,    15,    16,    17,    18,    19,    21,    20,    22,    10,
       0,     0,     0,     0,     0,     0,     0,     0,    11,     0,
       0,     0,     0,     0,     0,    12,     0,     0,     0,     0,
      13,     0,     0,    14
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -21,   -21,    10,   -10
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
       0,     1,    21,    14
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int8 yytable[] =
{
      22,    16,    23,    15,    17,    34,    35,    36,    37,    38,
      39,    13,    18,     0,    19,    41,    42,    43,    44,    45,
      46,    47,    20,    48,     2,    24,     0,     3,     4,     5,
       6,     7,     8,     9,     0,     0,     0,    10,     0,     0,
      54,    55,    56,    57,    11,    12,    29,    30,    31,    62,
      63,    64,     0,    33,     0,     0,    68,    69,     0,     0,
       0,    72,     3,     4,     5,     6,     7,     8,     9,     0,
       0,     0,    10,     0,     0,     0,    25,    26,    27,    28,
      12,    29,    30,    31,    32,     0,     0,     0,    33,    25,
      26,    27,    28,     0,    29,    30,    31,     0,     0,    40,
       0,    33,    25,    26,    27,    28,     0,    29,    30,    31,
       0,     0,    49,     0,    33,    25,    26,    27,    28,     0,
      29,    30,    31,     0,     0,     0,    50,    33,    25,    26,
      27,    28,     0,    29,    30,    31,     0,     0,     0,    51,
      33,    25,    26,    27,    28,     0,    29,    30,    31,     0,
       0,     0,    52,    33,    25,    26,    27,    28,     0,    29,
      30,    31,     0,     0,     0,    53,    33,    25,    26,    27,
      28,     0,    29,    30,    31,     0,     0,    58,     0,    33,
      25,    26,    27,    28,     0,    29,    30,    31,     0,     0,
       0,    59,    33,    25,    26,    27,    28,     0,    29,    30,
      31,     0,     0,     0,    60,    33,    25,    26,    27,    28,
       0,    29,    30,    31,     0,     0,     0,    61,    33,    25,
      26,    27,    28,     0,    29,    30,    31,     0,     0,    65,
       0,    33,    25,    26,    27,    28,     0,    29,    30,    31,
       0,     0,     0,    66,    33,    25,    26,    27,    28,     0,
      29,    30,    31,     0,     0,     0,    67,    33,    25,    26,
      27,    28,     0,    29,    30,    31,     0,     0,    70,     0,
      33,    25,    26,    27,    28,     0,    29,    30,    31,     0,
       0,     0,    71,    33,    25,    26,    27,    28,     0,    29,
      30,    31,     0,     0,    73,     0,    33,    25,    26,    27,
      28,     0,    29,    30,    31,     0,     0,    27,    28,    33,
      29,    30,    31,     0,     0,     0,     0,    33
};

static const yytype_int8 yycheck[] =
{
      10,    21,    12,    11,    21,    15,    16,    17,    18,    19,
      20,     1,    21,    -1,    21,    25,    26,    27,    28,    29,
      30,    31,    21,    33,     0,    20,    -1,     3,     4,     5,
       6,     7,     8,     9,    -1,    -1,    -1,    13,    -1,    -1,
      50,    51,    52,    53,    20,    21,    17,    18,    19,    59,
      60,    61,    -1,    24,    -1,    -1,    66,    67,    -1,    -1,
      -1,    71,     3,     4,     5,     6,     7,     8,     9,    -1,
      -1,    -1,    13,    -1,    -1,    -1,    12,    13,    14,    15,
      21,    17,    18,    19,    20,    -1,    -1,    -1,    24,    12,
      13,    14,    15,    -1,    17,    18,    19,    -1,    -1,    22,
      -1,    24,    12,    13,    14,    15,    -1,    17,    18,    19,
      -1,    -1,    22,    -1,    24,    12,    13,    14,    15,    -1,
      17,    18,    19,    -1,    -1,    -1,    23,    24,    12,    13,
      14,    15,    -1,    17,    18,    19,    -1,    -1,    -1,    23,
      24,    12,    13,    14,    15,    -1,    17,    18,    19,    -1,
      -1,    -1,    23,    24,    12,    13,    14,    15,    -1,    17,
      18,    19,    -1,    -1,    -1,    23,    24,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,    -1,    22,    -1,    24,
      12,    13,    14,    15,    -1,    17,    18,    19,    -1,    -1,
      -1,    23,    24,    12,    13,    14,    15,    -1,    17,    18,
      19,    -1,    -1,    -1,    23,    24,    12,    13,    14,    15,
      -1,    17,    18,    19,    -1,    -1,    -1,    23,    24,    12,
      13,    14,    15,    -1,    17,    18,    19,    -1,    -1,    22,
      -1,    24,    12,    13,    14,    15,    -1,    17,    18,    19,
      -1,    -1,    -1,    23,    24,    12,    13,    14,    15,    -1,
      17,    18,    19,    -1,    -1,    -1,    23,    24,    12,    13,
      14,    15,    -1,    17,    18,    19,    -1,    -1,    22,    -1,
      24,    12,    13,    14,    15,    -1,    17,    18,    19,    -1,
      -1,    -1,    23,    24,    12,    13,    14,    15,    -1,    17,
      18,    19,    -1,    -1,    22,    -1,    24,    12,    13,    14,
      15,    -1,    17,    18,    19,    -1,    -1,    14,    15,    24,
      17,    18,    19,    -1,    -1,    -1,    -1,    24
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int8 yystos[] =
{
       0,    26,     0,     3,     4,     5,     6,     7,     8,     9,
      13,    20,    21,    27,    28,    11,    21,    21,    21,    21,
      21,    27,    28,    28,    20,    12,    13,    14,    15,    17,
      18,    19,    20,    24,    28,    28,    28,    28,    28,    28,
      22,    28,    28,    28,    28,    28,    28,    28,    28,    22,
      23,    23,    23,    23,    28,    28,    28,    28,    22,    23,
      23,    23,    28,    28,    28,    22,    23,    23,    28,    28,
      22,    23,    28,    22
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr1[] =
{
       0,    25,    26,    26,    26,    26,    27,    28,    28,    28,
      28,    28,    28,    28,    28,    28,    28,    28,    28,    28,
      28,    28,    28,    28,    28
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     0,     2,     3,     3,     3,     1,     1,     1,
       4,     6,     8,    10,    12,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     2
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 5: /* list: list expr '\n'  */
#line 212 "initial.y"
                             { value = (yyvsp[-1].val); }
#line 1360 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 6: /* asgn: VAR '=' expr  */
#line 214 "initial.y"
                             { (yyval.val)=(yyvsp[-2].sym)->u.val=(yyvsp[0].val); (yyvsp[-2].sym)->type = VAR; }
#line 1366 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 8: /* expr: VAR  */
#line 217 "initial.y"
                             { if ((yyvsp[0].sym)->type == UNDEF) {
				 message ("yyparse: undefined variable ",
					  (yyvsp[0].sym)->name, WARNING);
			       }
			       (yyval.val) = (yyvsp[0].sym)->u.val;
			     }
#line 1377 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 10: /* expr: BLTIN_UNARY '(' expr ')'  */
#line 225 "initial.y"
          { (yyval.val) = (*((yyvsp[-3].sym)->u.ptr))((yyvsp[-1].val)); }
#line 1383 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 11: /* expr: BLTIN_BINARY '(' expr ',' expr ')'  */
#line 227 "initial.y"
          { (yyval.val) = (*((yyvsp[-5].sym)->u.ptr))((yyvsp[-3].val),(yyvsp[-1].val)); }
#line 1389 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 12: /* expr: BLTIN_TERNARY '(' expr ',' expr ',' expr ')'  */
#line 229 "initial.y"
          { (yyval.val) = (*((yyvsp[-7].sym)->u.ptr))((yyvsp[-5].val),(yyvsp[-3].val),(yyvsp[-1].val)); }
#line 1395 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 13: /* expr: BLTIN_QUATERNARY '(' expr ',' expr ',' expr ',' expr ')'  */
#line 231 "initial.y"
          { (yyval.val) = (*((yyvsp[-9].sym)->u.ptr))((yyvsp[-7].val),(yyvsp[-5].val),(yyvsp[-3].val),(yyvsp[-1].val)); }
#line 1401 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 14: /* expr: BLTIN_QUINTERNARY '(' expr ',' expr ',' expr ',' expr ',' expr ')'  */
#line 233 "initial.y"
          { (yyval.val) = (*((yyvsp[-11].sym)->u.ptr))((yyvsp[-9].val),(yyvsp[-7].val),(yyvsp[-5].val),(yyvsp[-3].val),(yyvsp[-1].val)); }
#line 1407 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 15: /* expr: expr '+' expr  */
#line 234 "initial.y"
                             { (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val); }
#line 1413 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 16: /* expr: expr '-' expr  */
#line 235 "initial.y"
                             { (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val); }
#line 1419 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 17: /* expr: expr '*' expr  */
#line 236 "initial.y"
                             { (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val); }
#line 1425 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 18: /* expr: expr '/' expr  */
#line 237 "initial.y"
                             { if ((yyvsp[0].val) == 0.0) {
          message ("yyparse", "division by zero", WARNING); (yyval.val) = 0.0;
	  } else { (yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val); }}
#line 1433 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 19: /* expr: expr '^' expr  */
#line 243 "initial.y"
                             { (yyval.val) = pow   ((yyvsp[-2].val), (yyvsp[0].val)); }
#line 1439 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 20: /* expr: expr '&' expr  */
#line 244 "initial.y"
                             { (yyval.val) = hypot ((yyvsp[-2].val), (yyvsp[0].val)); }
#line 1445 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 21: /* expr: expr '~' expr  */
#line 245 "initial.y"
                             { (yyval.val) = atan2 ((yyvsp[-2].val), (yyvsp[0].val)); }
#line 1451 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 22: /* expr: expr '%' expr  */
#line 246 "initial.y"
                             { (yyval.val) = fmod  ((yyvsp[-2].val), (yyvsp[0].val)); }
#line 1457 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 23: /* expr: '(' expr ')'  */
#line 247 "initial.y"
                             { (yyval.val) = (yyvsp[-1].val); }
#line 1463 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;

  case 24: /* expr: '-' expr  */
#line 248 "initial.y"
                                    { (yyval.val) = -(yyvsp[0].val); }
#line 1469 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"
    break;


#line 1473 "/Users/hmb/develop-git/semtex-xxt/build/femlib/initial.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 250 "initial.y"



void yy_initialize (void)
/* ------------------------------------------------------------------------- *
 * Load lookup tables and symbol table with default values.
 *
 * This routine should be called at start of run-time.
 * ------------------------------------------------------------------------- */
{
  static   int_t   initialized = 0;
  register int_t   i;
  register Symbol* s;

  if (!initialized) {
    for (i = 0; consts[i].name; i++)
      install (consts[i].name, VAR, consts[i].cval);

    for (i = 0; builtin[i].name; i++)
      switch (builtin[i].narg) {
      case 1:
	s = install (builtin[i].name, BLTIN_UNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 2:
	s = install (builtin[i].name, BLTIN_BINARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 3:
	s = install (builtin[i].name, BLTIN_TERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 4:
	s = install (builtin[i].name, BLTIN_QUATERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      case 5:
	s = install (builtin[i].name, BLTIN_QUINTERNARY, 0.0);
	s -> u.ptr = builtin[i].func;
	break;
      default:
	message ("yy_initialize", "never happen", ERROR);
	break;
      }

    initialized = 1;
  }
}


double yy_interpret (const char* s)
/* ------------------------------------------------------------------------- *
 * Given a string, interpret it as a function using yacc-generated yyparse.
 * ------------------------------------------------------------------------- */
{
  if (strlen (s) > STR_MAX)
    message ("yy_interpret: too many characters passed:\n", s, ERROR);
  
  strcat (strcpy (func_string, s), "\n");
  
  cur_string = func_string;
  
  yyparse ();
  return value;
}


void yy_vec_init (const char* names,
		  const char* fn   )
/* ------------------------------------------------------------------------- *
 * Set up the vector parser.
 *
 * names contains a list of variable names  e.g. "x y z",
 * fn    contains a function for evaluation e.g. "sin(x)*cos(y)*exp(z)".
 *
 * Valid separator characters in name are space, tab, comma, (semi-)colon.
 * Function string can contain previously-defined symbols (e.g. PI).
 * ------------------------------------------------------------------------- */
{
  char    routine   [] = "vecInit()";
  char    separator [] = " ,:;\t";
  char    tmp       [STR_MAX];
  char*   p;
  Symbol* s;

  if (strlen (fn) == 0)
    message (routine, "empty function string", ERROR);
  else if (strlen (fn) > STR_MAX)
    message (routine, "too many characters in function string", ERROR);

  nvec = 0;
  strcpy (tmp, names);

  p = strtok (tmp, separator);
  do {
    if (nvec++ > VEC_MAX) message (routine, "too many variables", ERROR); 
    vs[nvec-1] = (s = lookup (p)) ? s : install (p, VAR, 0.0);
  } while (p = strtok (NULL, separator));

  
  strcat (strcpy (func_string, fn), "\n");
}


void yy_vec_interp (const int_t ntot, ...)
/* ------------------------------------------------------------------------- *
 * Vector parser.  Following ntot there should be passed a number of
 * pointers to double (vectors), of which there should be in number the
 * number of variables named previously to vecInit, plus one: the result of
 * continually re-parsing the string "fn" is placed in the last vector, for
 * a total of ntot parsings.
 *
 * To follow on from the previous example, four vectors would be passed,
 * i.e.  vecInterp(ntot, x, y, z, u); the result fn(x,y,z) is placed in u.
 * ------------------------------------------------------------------------- */
{
  char           routine[] = "yy_vec_interp";
  register int_t i, n;
  double*        x[VEC_MAX];
  double*        fx = NULL;
  va_list        ap;
  
  va_start (ap, ntot);
  for (i = 0; i < nvec; i++) {
    x[i] = NULL;
    if (!(x[i] = va_arg (ap, double*)))
	message (routine, "not enough vectors passed..1", ERROR);
  }
  if (!(fx = va_arg (ap, double*)))
    message (routine, "not enough vectors passed..2", ERROR);
  va_end (ap);

  for (n = 0; n < ntot; n++) {
    cur_string = func_string;
    for (i = 0; i < nvec; i++) vs[i]->u.val = x[i][n];
    yyparse ();
    fx[n] = value;
  }
}


void yy_help (void)
/* ------------------------------------------------------------------------- *
 * Print details of callable functions to stderr.
 * ------------------------------------------------------------------------- */
{
  fprintf 
    (stderr, 
     "Unary:      -\n"
     "Binary:     -, +, *, /, ^ (exponentiation), "
     "~ (atan2), & (hypot), %% (fmod)\n"
     "Functions:  sin,  cos,  tan,  asin,  acos,  atan,\n"
     "            sinh, cosh, tanh, asinh, acosh, atanh,\n"
     "            abs, floor, ceil, int,   heav (Heaviside),\n"
     "            log, log10, exp,  sqrt,  white\n"
     "            erf, erfc, gamma, lgamma, sgn,\n"
     "            j0, j1, y0, y1,\n"
     "Procedures: step, jn, yn, rad, ang, rejn,imjn, jacobi, womcos,womsin\n");
}


void yy_show (void)
/* ------------------------------------------------------------------------- *
 * Print details of installed variables to stderr.
 * ------------------------------------------------------------------------- */
{
  register int_t   i;
  register Symbol* sp;

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) 
	fprintf (stderr, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
}


int_t yy_dump (char*         str,
		 const int_t max)
/* ------------------------------------------------------------------------- *
 * Load description of internal variables into string, to length max.
 * If string overflows, return 0, else 1.
 * ------------------------------------------------------------------------- */
{
  register int_t   i, n = 0;
  register Symbol* sp;
  char             buf[FILENAME_MAX];

  for (i = 0; i < HASHSIZE; i++)
    for (sp = hashtab[i]; sp; sp = sp -> next)
      if (sp -> type == VAR) {
	sprintf (buf, "%-15s = %-.17g\n", sp -> name, sp -> u.val);
	if ((n += strlen (buf)) > max - 2)
	  return 0;
	else
	  strcat (str, buf);
      }

  return 1;
}


static int yylex (void)
/* ------------------------------------------------------------------------- *
 * Lexical analysis routine called by yyparse, using string loaded by
 * yy_interpret.
 * ------------------------------------------------------------------------- */
{
  register int_t c;

  while ((c = *cur_string++) == ' ' || c == '\t');

  if (c == EOF) return 0;

  if (c == '.' || isdigit (c)) {
    yylval.val = strtod (--cur_string, &cur_string);
    return NUMBER;
  }

  if (isalpha (c)) {
    register Symbol* s;
    char             sbuf[STR_MAX];
    register char*   p = sbuf;
    do {
      *p++ = c;
    } while ((c = *cur_string++) != EOF && (isalnum (c) || c == '_'));
    cur_string--;
    *p = '\0';
    if ((s = lookup (sbuf)) == NULL) s = install (sbuf, UNDEF, 0.0);
    yylval.sym = s;
    return (s -> type == UNDEF) ? VAR : s -> type;
  }

  return c;
}


static void yyerror (char *s)
/* ------------------------------------------------------------------------- *
 * Handler for yyparse syntax errors.
 * ------------------------------------------------------------------------- */
{
  message ("yyparse", s, WARNING);
}


static unsigned hash (const char* s)
/* ------------------------------------------------------------------------- *
 * Generate hash table index.
 * ------------------------------------------------------------------------- */
{
  register unsigned hashval;

  for (hashval = 0; *s != '\0'; s++) hashval = *s + HASHSEED * hashval;
  
  return hashval % HASHSIZE;
}


static Symbol* lookup (const char* s)
/* ------------------------------------------------------------------------- *
 * Find s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol* sp;
  
  for (sp = hashtab[hash (s)]; sp; sp = sp->next)
    if (strcmp (s, sp->name) == 0) return sp;

  return NULL;
}


static Symbol* install (const char*   s,
			const int_t t,
			const double  d)
/* ------------------------------------------------------------------------- *
 * Install s in symbol hashtable.
 * ------------------------------------------------------------------------- */
{
  register Symbol*  sp;
  register unsigned hashval;

  if (!(sp = lookup (s))) {	/* -- Not found, install in hashtab. */
    sp = (Symbol *) emalloc (sizeof (Symbol));
    if (sp == NULL || (sp -> name = strdup (s)) == NULL) return NULL;
    hashval          = hash (s);
    sp -> next       = hashtab[hashval];
    hashtab[hashval] = sp;
  }

  sp -> type  = t;
  sp -> u.val = d;
  
  return sp;
}


static void *emalloc (const size_t n)
/* ------------------------------------------------------------------------- *
 * Check return from malloc.
 * ------------------------------------------------------------------------- */
{
  void* p;

  if (!(p = (void *) malloc (n))) message ("emalloc", "out of memory", ERROR);
 
  return p;
}

static double Sgn   (double x) { return copysign (1.0, x) ; }
static double Heavi (double x) { return (x >= 0.0) ? 1.0 : 0.0; }
static double White (double x) { return dnormal (0.0, x); }
static double Step  (double x, double a) { return (x >= a) ? 1.0 : 0.0; }

static double Rad (double x, double y) { return hypot (x, y);  }
static double Ang (double x, double y) { return atan2 (y, x);  }
static double Jn  (double i, double x) { return jn((int_t)i, x); }
static double Yn  (double i, double x) { return yn((int_t)i, x); }

static double Jacobi (double z, double n, double alpha, double beta)
/* ------------------------------------------------------------------------- *
 * Return value of the n_th order Jacobi polynomial
 *   P^(alpha,beta)_n(z) alpha > -1, beta > -1 at z.
 * ------------------------------------------------------------------------- */
{
  const double     apb = alpha + beta;
  register int_t i,k;
  double           a1,a2,a3,a4;
  double           poly, polyn1, polyn2;
  
  polyn2 = 1.0;
  polyn1 = 0.5*(alpha - beta + (alpha + beta + 2.0)*z);
  
  for (k = 2; k <= n; ++k){
    a1 =  2.0*k*(k + apb)*(2.0*k + apb - 2.0);
    a2 = (2.0*k + apb - 1.0)*(alpha*alpha - beta*beta);
    a3 = (2.0*k + apb - 2.0)*(2.0*k + apb - 1.0)*(2.0*k + apb);
    a4 =  2.0*(k + alpha - 1.0)*(k + beta - 1.0)*(2.0*k + apb);
    
    a2 /= a1;
    a3 /= a1;
    a4 /= a1;
    
    poly   = (a2 + a3*z)*polyn1 - a4*polyn2;
    polyn2 = polyn1;
    polyn1 = poly  ;
  }

  return poly;
}

/* -- Gamma function, ex netlib: C "standard" routine not so standard. */

double F77NAME(dgamma) (const double*);

double Gamma(const double x) { double X = x; return F77NAME(dgamma) (&X); }

/* -- Complex Bessel function, ex netlib, and functions that use it. */

void F77NAME(zbesj) (const double*, const double*, const double*, 
		     const int_t*, const int_t*, double*, double*, 
		     int_t*, int_t*);

void Zbesj (const double *x, const double *y, const double ord, 
	    const int_t Kode, const int_t n, double *ReJ, 
	    double *ImJ, int_t* nz, int_t* ierr) {
  int_t N = n;
  int_t K = Kode;
  double  order = ord;

  F77NAME(zbesj) (x, y, &order, &K, &N, ReJ, ImJ, nz, ierr);
  if (*ierr) message ("initial.y", "Zbesj", ERROR);
}

static double ReJn (double n, double x,  double y)
{
  int_t nz, ierr;
  double  rej, imj;
  
  Zbesj (&x,&y,n,1,1,&rej,&imj,&nz,&ierr);
  return rej;
}

static double ImJn (double n, double x, double y)
{
  int_t nz, ierr;
  double  rej, imj;
  
  Zbesj (&x,&y,n,1,1,&rej,&imj,&nz,&ierr);
  return imj;
}

static double Womersley (double A    ,
			 double B    ,
			 double r    ,
			 double R    , 
			 double mu   , 
			 double omega,
			 double t    )
/* -------------------------------------------------------------------------
 * Calculate the Womersley solution at r for a pipe of radius R and circular
 * frequency omega.  The solution is set so that the area-average flow speed
 *
 *                u_avg(t) = A cos(omega t) +  B sin(omega t)
 * ------------------------------------------------------------------------- */
{
  double x,y;

  if (r > R) message ("Womersley", "r > R", ERROR);

  if (omega == 0) /* Return Poiseuille flow with mean of 1. */
    return 2*(1-r*r/R/R);
  else {
    int_t  ierr,nz;
    double cr,ci,J0r,J0i,rej,imj,re,im,fac;
    double isqrt2 = 1.0/sqrt(2.0);
    static double R_str, omega_str,mu_str;
    static double Jr,Ji,alpha,j0r,j0i, isqrt;

    /* 
       For repeated calls with same parameters, store those independent of r.
    */

    if ((R != R_str) || (omega != omega_str) || (mu != mu_str)) {
      double retmp[2],imtmp[2];
      alpha = R*sqrt(omega/mu);

      re  = -alpha*isqrt2;
      im  =  alpha*isqrt2;
      Zbesj(&re,&im,0,1,2,retmp,imtmp,&nz,&ierr);
      j0r = retmp[0]; j0i = imtmp[0];
      rej = retmp[1]; imj = imtmp[1];

      fac = 1/(j0r*j0r+j0i*j0i);
      Jr = 1+2*fac/alpha*((rej*j0r+imj*j0i)*isqrt2-(imj*j0r - rej*j0i)*isqrt2);
      Ji = 2*fac/alpha*((rej*j0r+imj*j0i)*isqrt2 + (imj*j0r - rej*j0i)*isqrt2);

      R_str = R; omega_str = omega; mu_str = mu;
    }

    /* setup cr, ci from pre-stored value of Jr & Ji */

    fac = 1/(Jr*Jr + Ji*Ji);
    cr  =  (A*Jr - B*Ji)*fac;
    ci  = -(A*Ji + B*Jr)*fac;
    
    /* setup J0r, J0i */

    re  = -alpha*isqrt2*r/R;
    im  =  alpha*isqrt2*r/R;
    Zbesj(&re,&im,0,1,1,&rej,&imj,&nz,&ierr);
    fac = 1/(j0r*j0r+j0i*j0i);
    J0r = 1-fac*(rej*j0r+imj*j0i);
    J0i = -fac*(imj*j0r-rej*j0i);

    return (cr*J0r - ci*J0i)*cos(omega*t) - (ci*J0r + cr*J0i)*sin(omega*t);
  }
}

static double Womsin (double r, double R, double mu, double omega, double t) {
  if (omega == 0) return 0.0; /* No sin term for zeroth mode. */
  
  return Womersley (0.0, 1.0, r, R, mu, omega, t);
}

static double Womcos (double r, double R, double mu, double omega, double t) {
  return Womersley (1.0, 0.0, r, R, mu, omega, t);
}

/* ------------------------------------------------------------------------- *
 * Error checking is done for cases where we can have illegal input
 * values (typically negative), otherwise we accept exception returns.
 *
 * 1/6/2001: After much trouble with underflow errors, I've turned these off.
 * We just accept whatever comes back (typically a nan).
 *
 * example: static double Log (double x) { return errcheck (log (x), "log"); }
 * ------------------------------------------------------------------------- */


static double errcheck (const double d,
			const char*  s)
/* ------------------------------------------------------------------------- *
 * Check result of math library call. We don't use this anymore!
 * ------------------------------------------------------------------------- */
{
  if (errno == EDOM) {
    errno = 0;
    message ("errcheck: argument out of domain in call to", s, ERROR);
  } else if (errno == ERANGE) {
    errno = 0;
    message ("errcheck: result out of range in call to", s, ERROR);
  }

  return d;
}

#undef HASHSIZE
#undef HASHSEED
#undef VEC_MAX
