#include <stdlib.h>
#ifndef lint
#ifdef __unused
__unused
#endif
static char const 
yyrcsid[] = "$FreeBSD: src/usr.bin/yacc/skeleton.c,v 1.37 2003/02/12 18:03:55 davidc Exp $";
#endif
#define YYBYACC 1
#define YYMAJOR 1
#define YYMINOR 9
#define YYLEX yylex()
#define YYEMPTY -1
#define yyclearin (yychar=(YYEMPTY))
#define yyerrok (yyerrflag=0)
#define YYRECOVERING() (yyerrflag!=0)
#if defined(__cplusplus) || __STDC__
static int yygrowstack(void);
#else
static int yygrowstack();
#endif
#define YYPREFIX "yy"
#line 7 "src/calc.y"
#include <sim/sad_f2c.h>
#include <math.h>
#define	YYMAXDEPTH	256
#define	YYSTYPE		real8
extern YYSTYPE yyval_cb;	/* yyval copy back*/
int  yylex(void);
void yyerror(const char *);
void yyunget(void);
#line 33 "src/calc.c"
#define YYERRCODE 256
#define NUM 257
#define ID 258
#define UNIT 259
#define EOL 260
#define INVALID 261
#define FNSQRT 262
#define FNEXP 263
#define FNLN 264
#define FNLOG 265
#define FNSIN 266
#define FNCOS 267
#define FNTAN 268
#define FNATAN 269
#define FNATAN2 270
#define FNSINH 271
#define FNCOSH 272
#define FNTANH 273
#define NEG 274
const short yylhs[] = {                                        -1,
    1,    1,    1,    2,    2,    0,    0,    0,    3,    3,
    4,    4,    5,    5,    5,    5,    5,    5,    5,    6,
    6,    6,    6,    6,    7,    7,    7,    8,    8,    8,
    8,    8,    8,    8,    9,    9,    9,    9,    9,    9,
    9,    9,    9,    9,    9,    9,
};
const short yylen[] = {                                         2,
    1,    1,    1,    0,    1,    1,    2,    1,    1,    3,
    1,    3,    1,    3,    3,    3,    4,    4,    4,    1,
    4,    4,    2,    2,    1,    4,    4,    1,    2,    1,
    3,    4,    2,    1,    4,    4,    4,    4,    4,    4,
    4,    4,    4,    4,    4,    6,
};
const short yydefred[] = {                                      0,
    0,   30,    6,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,   34,   29,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,   33,    0,    3,    1,    2,    0,    7,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   31,    0,    0,    0,    0,    0,    0,    0,    0,    5,
    0,    0,    0,    0,    0,   35,   36,   37,   38,   39,
   40,   41,   45,    0,   42,   43,   44,    0,    0,    0,
    0,    0,    0,    0,    0,    0,   46,
};
const short yydgoto[] = {                                      20,
   49,   81,   21,   22,   23,   24,   25,   26,   27,
};
const short yysindex[] = {                                    345,
 -234,    0,    0,  -14,  -13,  -11,  -10,   -6,   -5,   -4,
   -2,   15,   17,   18,   19,  268,  268,  268,  381,    0,
  -39,    2,  -29,  -21,  -19,  -53,    0,    0,  381,  381,
  381,  381,  381,  381,  381,  381,  381,  381,  381,  381,
  -19,  -19,    0,  -41,    0,    0,    0,  381,    0,  381,
  322,  362,  381, -211, -211, -211, -211, -211,  -38,  -37,
  -36,  -35,  -34,  -33,  -31,  -30,  -32,  -24,  -23,  -22,
    0,    2,  -29,  381,  381,  -21,  381,  -21,  -21,    0,
  268,  268,  268,  268,  268,    0,    0,    0,    0,    0,
    0,    0,    0,  381,    0,    0,    0,  -21,  -21,  -21,
  -19,  -19,  -53,  -53,  -53,  -20,    0,
};
const short yyrindex[] = {                                      0,
    1,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   64,  153,  307,  192,  103,    9,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
  113,  145,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,  401,  401,  401,  401,  401,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,  202,  311,    0,    0,  222,    0,  232,  258,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,  262,  271,  298,
  155,  180,   37,   67,   77,    0,    0,
};
const short yygindex[] = {                                      0,
    0,  -42,  463,   24,   16,  323,  100,   47,    0,
};
#define YYTABLESIZE 674
const short yytable[] = {                                      71,
   28,   47,   86,   87,   88,   89,   90,   91,   25,   92,
   93,   94,   82,   83,   84,   85,   95,   96,   97,   46,
  107,   55,   56,   54,   28,   29,   30,   57,   31,   32,
   51,   53,   52,   33,   34,   35,   26,   36,   28,   50,
   58,   28,   28,   28,   28,   28,   25,   28,   80,   25,
   25,   25,   25,   25,   37,   25,   38,   39,   40,   28,
   28,   28,   28,    8,   43,   73,   27,   25,   25,   25,
   25,   72,    0,    0,   26,    0,   32,   26,   26,   26,
   26,   26,   48,   26,   48,   48,   48,   48,   48,   48,
   48,   48,   48,   48,   28,   26,   26,   26,   26,   48,
   48,   48,   20,   48,   27,    0,    0,   27,   27,   27,
   27,   27,   23,   27,   32,   41,   42,   32,   32,   32,
   32,   32,    0,   32,   28,   27,   27,   27,   27,  103,
  104,  105,   25,    0,    0,   32,   32,   32,   32,    0,
   20,    0,    0,   20,   24,   20,   20,   20,    0,    0,
   23,    0,    9,   23,   22,   23,   23,   23,    0,    0,
   26,   20,   20,   20,   20,    0,    0,    0,    0,    0,
    0,   23,   23,   23,   23,    0,    0,    0,    0,   21,
  101,  102,   24,    0,    0,   24,    0,   24,   24,   24,
   27,   13,   22,    9,    0,   22,    9,   22,   22,   22,
   32,   10,    0,   24,   24,   24,   24,    0,    0,    0,
    0,    9,    0,   22,   22,   22,   22,   21,   45,    0,
   21,   14,   21,   21,   21,    0,   20,    0,    0,   13,
    0,   15,   13,    0,    0,   13,   23,    0,   21,   21,
   21,   21,   10,    0,    0,   10,    0,    0,    0,    0,
   13,   13,   13,   13,    0,    0,    0,   16,   28,   14,
   10,   19,   14,    0,    0,   14,   25,    0,   24,   15,
   17,    0,   15,    0,    0,   15,    9,    0,   22,    0,
   14,   14,   14,   14,    0,    0,    0,    0,    0,    0,
   15,   15,   15,   15,   26,   16,    0,   18,   16,   19,
    0,   16,   19,   21,    0,   19,   11,   19,   17,    0,
   12,   17,    0,    0,   17,   13,   16,   16,   16,   16,
   19,   19,   19,   19,   27,   10,    0,    0,    0,   17,
   17,   17,   17,    0,   32,   18,    0,    0,   18,    0,
    0,   18,    0,    0,   11,   14,    0,   11,   12,    0,
   11,   12,    0,    0,   12,   15,   18,   18,   18,   18,
   20,   19,    0,    0,   17,   11,   16,    0,    0,   12,
   23,    0,    0,   76,   78,   79,    0,    0,    0,    0,
    0,   16,   75,   74,   19,   19,    0,   17,    0,   16,
    0,    0,    0,   18,   17,    0,   98,   99,    0,  100,
    0,   19,   24,    0,   17,    0,   16,    0,    0,    0,
    9,    0,   22,    0,    0,    0,    0,    0,    0,    0,
   19,   18,   77,   17,    0,   16,    0,    0,    0,    0,
   11,    0,    0,    0,   12,    0,    0,   21,    0,    0,
    4,    0,    0,    0,    0,    0,    0,   18,    0,   13,
    0,    0,    0,    0,    0,    0,    0,    0,    0,   10,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
   18,    0,    0,    0,    0,    0,    0,    0,    0,   14,
    0,   44,    0,    0,    0,    0,    0,   18,    0,   15,
    0,   59,   60,   61,   62,   63,   64,   65,   66,   67,
   68,   69,   70,    0,    0,    0,   18,    0,    0,    0,
    0,    0,    0,    0,    0,   16,    0,    0,    0,   19,
    0,    0,    0,    0,    1,    2,    4,    0,   17,    4,
    5,    6,    7,    8,    9,   10,   11,   12,   13,   14,
   15,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,   18,  106,    0,    0,    0,
    0,    0,    0,    0,   11,    0,    0,    0,   12,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    1,    2,
    0,    0,    0,    4,    5,    6,    7,    8,    9,   10,
   11,   12,   13,   14,   15,    0,    0,    0,    0,    0,
    0,    1,    2,    3,    0,    0,    4,    5,    6,    7,
    8,    9,   10,   11,   12,   13,   14,   15,    1,    2,
    0,    0,    0,    4,    5,    6,    7,    8,    9,   10,
   11,   12,   13,   14,   15,    0,    0,    1,    2,    0,
    0,    0,    4,    5,    6,    7,    8,    9,   10,   11,
   12,   13,   14,   15,    0,    0,    0,    4,    4,    0,
    0,    0,    4,    4,    4,    4,    4,    4,    4,    4,
    4,    4,    4,    4,
};
const short yycheck[] = {                                      41,
    0,   41,   41,   41,   41,   41,   41,   41,    0,   41,
   41,   44,   55,   56,   57,   58,   41,   41,   41,   59,
   41,   43,   42,   45,  259,   40,   40,   47,   40,   40,
   60,   61,   62,   40,   40,   40,    0,   40,   38,   38,
   94,   41,   42,   43,   44,   45,   38,   47,  260,   41,
   42,   43,   44,   45,   40,   47,   40,   40,   40,   59,
   60,   61,   62,    0,   18,   50,    0,   59,   60,   61,
   62,   48,   -1,   -1,   38,   -1,    0,   41,   42,   43,
   44,   45,  124,   47,  124,  124,  124,  124,  124,  124,
  124,  124,  124,  124,   94,   59,   60,   61,   62,  124,
  124,  124,    0,  124,   38,   -1,   -1,   41,   42,   43,
   44,   45,    0,   47,   38,   16,   17,   41,   42,   43,
   44,   45,   -1,   47,  124,   59,   60,   61,   62,   83,
   84,   85,  124,   -1,   -1,   59,   60,   61,   62,   -1,
   38,   -1,   -1,   41,    0,   43,   44,   45,   -1,   -1,
   38,   -1,    0,   41,    0,   43,   44,   45,   -1,   -1,
  124,   59,   60,   61,   62,   -1,   -1,   -1,   -1,   -1,
   -1,   59,   60,   61,   62,   -1,   -1,   -1,   -1,    0,
   81,   82,   38,   -1,   -1,   41,   -1,   43,   44,   45,
  124,    0,   38,   41,   -1,   41,   44,   43,   44,   45,
  124,    0,   -1,   59,   60,   61,   62,   -1,   -1,   -1,
   -1,   59,   -1,   59,   60,   61,   62,   38,  258,   -1,
   41,    0,   43,   44,   45,   -1,  124,   -1,   -1,   38,
   -1,    0,   41,   -1,   -1,   44,  124,   -1,   59,   60,
   61,   62,   41,   -1,   -1,   44,   -1,   -1,   -1,   -1,
   59,   60,   61,   62,   -1,   -1,   -1,    0,  258,   38,
   59,    0,   41,   -1,   -1,   44,  258,   -1,  124,   38,
    0,   -1,   41,   -1,   -1,   44,  124,   -1,  124,   -1,
   59,   60,   61,   62,   -1,   -1,   -1,   -1,   -1,   -1,
   59,   60,   61,   62,  258,   38,   -1,    0,   41,   38,
   -1,   44,   41,  124,   -1,   44,    0,   40,   38,   -1,
    0,   41,   -1,   -1,   44,  124,   59,   60,   61,   62,
   59,   60,   61,   62,  258,  124,   -1,   -1,   -1,   59,
   60,   61,   62,   -1,  258,   38,   -1,   -1,   41,   -1,
   -1,   44,   -1,   -1,   38,  124,   -1,   41,   38,   -1,
   44,   41,   -1,   -1,   44,  124,   59,   60,   61,   62,
  258,   40,   -1,   -1,   43,   59,   45,   -1,   -1,   59,
  258,   -1,   -1,   51,   52,   53,   -1,   -1,   -1,   -1,
   -1,  124,   61,   62,   40,  124,   -1,   43,   -1,   45,
   -1,   -1,   -1,  126,  124,   -1,   74,   75,   -1,   77,
   -1,   40,  258,   -1,   43,   -1,   45,   -1,   -1,   -1,
  258,   -1,  258,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   40,  124,   61,   43,   -1,   45,   -1,   -1,   -1,   -1,
  124,   -1,   -1,   -1,  124,   -1,   -1,  258,   -1,   -1,
   40,   -1,   -1,   -1,   -1,   -1,   -1,  126,   -1,  258,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  258,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
  126,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  258,
   -1,   19,   -1,   -1,   -1,   -1,   -1,  126,   -1,  258,
   -1,   29,   30,   31,   32,   33,   34,   35,   36,   37,
   38,   39,   40,   -1,   -1,   -1,  126,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,  258,   -1,   -1,   -1,  258,
   -1,   -1,   -1,   -1,  257,  258,  126,   -1,  258,  262,
  263,  264,  265,  266,  267,  268,  269,  270,  271,  272,
  273,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,   -1,  258,   94,   -1,   -1,   -1,
   -1,   -1,   -1,   -1,  258,   -1,   -1,   -1,  258,   -1,
   -1,   -1,   -1,   -1,   -1,   -1,   -1,   -1,  257,  258,
   -1,   -1,   -1,  262,  263,  264,  265,  266,  267,  268,
  269,  270,  271,  272,  273,   -1,   -1,   -1,   -1,   -1,
   -1,  257,  258,  259,   -1,   -1,  262,  263,  264,  265,
  266,  267,  268,  269,  270,  271,  272,  273,  257,  258,
   -1,   -1,   -1,  262,  263,  264,  265,  266,  267,  268,
  269,  270,  271,  272,  273,   -1,   -1,  257,  258,   -1,
   -1,   -1,  262,  263,  264,  265,  266,  267,  268,  269,
  270,  271,  272,  273,   -1,   -1,   -1,  257,  258,   -1,
   -1,   -1,  262,  263,  264,  265,  266,  267,  268,  269,
  270,  271,  272,  273,
};
#define YYFINAL 20
#ifndef YYDEBUG
#define YYDEBUG 0
#endif
#define YYMAXTOKEN 274
#if YYDEBUG
const char * const yyname[] = {
"end-of-file",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,"'&'",0,"'('","')'","'*'","'+'","','","'-'",0,"'/'",0,0,0,0,0,0,0,0,0,0,
0,"';'","'<'","'='","'>'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,"'^'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"'|'",0,
"'~'",0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,"NUM","ID","UNIT","EOL","INVALID","FNSQRT","FNEXP",
"FNLN","FNLOG","FNSIN","FNCOS","FNTAN","FNATAN","FNATAN2","FNSINH","FNCOSH",
"FNTANH","NEG",
};
const char * const yyrule[] = {
"$accept : expr",
"end_of_expr : ';'",
"end_of_expr : ')'",
"end_of_expr : ID",
"opt_eol :",
"opt_eol : EOL",
"expr : UNIT",
"expr : or_expr end_of_expr",
"expr : or_expr",
"or_expr : and_expr",
"or_expr : or_expr '|' and_expr",
"and_expr : cmp_expr",
"and_expr : and_expr '&' cmp_expr",
"cmp_expr : a_expr",
"cmp_expr : cmp_expr '<' a_expr",
"cmp_expr : cmp_expr '>' a_expr",
"cmp_expr : cmp_expr '=' a_expr",
"cmp_expr : cmp_expr '<' '=' a_expr",
"cmp_expr : cmp_expr '>' '=' a_expr",
"cmp_expr : cmp_expr '<' '>' a_expr",
"a_expr : term",
"a_expr : a_expr '+' opt_eol term",
"a_expr : a_expr '-' opt_eol term",
"a_expr : '-' term",
"a_expr : '+' term",
"term : factor",
"term : term '*' opt_eol factor",
"term : term '/' opt_eol factor",
"factor : NUM",
"factor : NUM UNIT",
"factor : ID",
"factor : '(' or_expr ')'",
"factor : factor '^' opt_eol factor",
"factor : '~' factor",
"factor : func_call",
"func_call : FNSQRT '(' or_expr ')'",
"func_call : FNEXP '(' or_expr ')'",
"func_call : FNLN '(' or_expr ')'",
"func_call : FNLOG '(' or_expr ')'",
"func_call : FNSIN '(' or_expr ')'",
"func_call : FNCOS '(' or_expr ')'",
"func_call : FNTAN '(' or_expr ')'",
"func_call : FNSINH '(' or_expr ')'",
"func_call : FNCOSH '(' or_expr ')'",
"func_call : FNTANH '(' or_expr ')'",
"func_call : FNATAN '(' or_expr ')'",
"func_call : FNATAN2 '(' or_expr ',' or_expr ')'",
};
#endif
#ifndef YYSTYPE
typedef int YYSTYPE;
#endif
#if YYDEBUG
#include <stdio.h>
#endif
#ifdef YYSTACKSIZE
#undef YYMAXDEPTH
#define YYMAXDEPTH YYSTACKSIZE
#else
#ifdef YYMAXDEPTH
#define YYSTACKSIZE YYMAXDEPTH
#else
#define YYSTACKSIZE 10000
#define YYMAXDEPTH 10000
#endif
#endif
#define YYINITSTACKSIZE 200
int yydebug;
int yynerrs;
int yyerrflag;
int yychar;
short *yyssp;
YYSTYPE *yyvsp;
YYSTYPE yyval;
YYSTYPE yylval;
short *yyss;
short *yysslim;
YYSTYPE *yyvs;
int yystacksize;
/* allocate initial stack or double stack size, up to YYMAXDEPTH */
static int yygrowstack()
{
    int newsize, i;
    short *newss;
    YYSTYPE *newvs;

    if ((newsize = yystacksize) == 0)
        newsize = YYINITSTACKSIZE;
    else if (newsize >= YYMAXDEPTH)
        return -1;
    else if ((newsize *= 2) > YYMAXDEPTH)
        newsize = YYMAXDEPTH;
    i = yyssp - yyss;
    newss = yyss ? (short *)realloc(yyss, newsize * sizeof *newss) :
      (short *)malloc(newsize * sizeof *newss);
    if (newss == NULL)
        return -1;
    yyss = newss;
    yyssp = newss + i;
    newvs = yyvs ? (YYSTYPE *)realloc(yyvs, newsize * sizeof *newvs) :
      (YYSTYPE *)malloc(newsize * sizeof *newvs);
    if (newvs == NULL)
        return -1;
    yyvs = newvs;
    yyvsp = newvs + i;
    yystacksize = newsize;
    yysslim = yyss + newsize - 1;
    return 0;
}

#define YYABORT goto yyabort
#define YYREJECT goto yyabort
#define YYACCEPT goto yyaccept
#define YYERROR goto yyerrlab

#ifndef YYPARSE_PARAM
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG void
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG
#define YYPARSE_PARAM_DECL
#endif	/* ANSI-C/C++ */
#else	/* YYPARSE_PARAM */
#ifndef YYPARSE_PARAM_TYPE
#define YYPARSE_PARAM_TYPE void *
#endif
#if defined(__cplusplus) || __STDC__
#define YYPARSE_PARAM_ARG YYPARSE_PARAM_TYPE YYPARSE_PARAM
#define YYPARSE_PARAM_DECL
#else	/* ! ANSI-C/C++ */
#define YYPARSE_PARAM_ARG YYPARSE_PARAM
#define YYPARSE_PARAM_DECL YYPARSE_PARAM_TYPE YYPARSE_PARAM;
#endif	/* ANSI-C/C++ */
#endif	/* ! YYPARSE_PARAM */

int
yyparse (YYPARSE_PARAM_ARG)
    YYPARSE_PARAM_DECL
{
    int yym, yyn, yystate;
#if YYDEBUG
    const char *yys;

    if ((yys = getenv("YYDEBUG")))
    {
        yyn = *yys;
        if (yyn >= '0' && yyn <= '9')
            yydebug = yyn - '0';
    }
#endif

    yynerrs = 0;
    yyerrflag = 0;
    yychar = (-1);

    if (yyss == NULL && yygrowstack()) goto yyoverflow;
    yyssp = yyss;
    yyvsp = yyvs;
    *yyssp = yystate = 0;

yyloop:
    if ((yyn = yydefred[yystate])) goto yyreduce;
    if (yychar < 0)
    {
        if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, reading %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
    }
    if ((yyn = yysindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: state %d, shifting to state %d\n",
                    YYPREFIX, yystate, yytable[yyn]);
#endif
        if (yyssp >= yysslim && yygrowstack())
        {
            goto yyoverflow;
        }
        *++yyssp = yystate = yytable[yyn];
        *++yyvsp = yylval;
        yychar = (-1);
        if (yyerrflag > 0)  --yyerrflag;
        goto yyloop;
    }
    if ((yyn = yyrindex[yystate]) && (yyn += yychar) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yychar)
    {
        yyn = yytable[yyn];
        goto yyreduce;
    }
    if (yyerrflag) goto yyinrecovery;
#if defined(lint) || defined(__GNUC__)
    goto yynewerror;
#endif
yynewerror:
    yyerror("syntax error");
#if defined(lint) || defined(__GNUC__)
    goto yyerrlab;
#endif
yyerrlab:
    ++yynerrs;
yyinrecovery:
    if (yyerrflag < 3)
    {
        yyerrflag = 3;
        for (;;)
        {
            if ((yyn = yysindex[*yyssp]) && (yyn += YYERRCODE) >= 0 &&
                    yyn <= YYTABLESIZE && yycheck[yyn] == YYERRCODE)
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: state %d, error recovery shifting\
 to state %d\n", YYPREFIX, *yyssp, yytable[yyn]);
#endif
                if (yyssp >= yysslim && yygrowstack())
                {
                    goto yyoverflow;
                }
                *++yyssp = yystate = yytable[yyn];
                *++yyvsp = yylval;
                goto yyloop;
            }
            else
            {
#if YYDEBUG
                if (yydebug)
                    printf("%sdebug: error recovery discarding state %d\n",
                            YYPREFIX, *yyssp);
#endif
                if (yyssp <= yyss) goto yyabort;
                --yyssp;
                --yyvsp;
            }
        }
    }
    else
    {
        if (yychar == 0) goto yyabort;
#if YYDEBUG
        if (yydebug)
        {
            yys = 0;
            if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
            if (!yys) yys = "illegal-symbol";
            printf("%sdebug: state %d, error recovery discards token %d (%s)\n",
                    YYPREFIX, yystate, yychar, yys);
        }
#endif
        yychar = (-1);
        goto yyloop;
    }
yyreduce:
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: state %d, reducing by rule %d (%s)\n",
                YYPREFIX, yystate, yyn, yyrule[yyn]);
#endif
    yym = yylen[yyn];
    yyval = yyvsp[1-yym];
    switch (yyn)
    {
case 6:
#line 30 "src/calc.y"
{yyerror("Unexpected UNIT");}
break;
case 7:
#line 31 "src/calc.y"
{yyval = yyvsp[-1]; yyunget(); yyval_cb = yyval; YYACCEPT;}
break;
case 8:
#line 32 "src/calc.y"
{yyval = yyvsp[0]; yyunget(); yyval_cb = yyval;}
break;
case 10:
#line 36 "src/calc.y"
{yyval = (yyvsp[-2] != 0 || yyvsp[0] != 0) ? 1 : 0;}
break;
case 12:
#line 40 "src/calc.y"
{yyval = (yyvsp[-2] != 0 && yyvsp[0] != 0) ? 1 : 0;}
break;
case 14:
#line 44 "src/calc.y"
{yyval = ( yyvsp[-2] <  yyvsp[0]) ? 1 : 0;}
break;
case 15:
#line 45 "src/calc.y"
{yyval = ( yyvsp[-2] >  yyvsp[0]) ? 1 : 0;}
break;
case 16:
#line 46 "src/calc.y"
{yyval = ( yyvsp[-2] == yyvsp[0]) ? 1 : 0;}
break;
case 17:
#line 47 "src/calc.y"
{yyval = ( yyvsp[-3] <= yyvsp[0]) ? 1 : 0;}
break;
case 18:
#line 48 "src/calc.y"
{yyval = ( yyvsp[-3] >= yyvsp[0]) ? 1 : 0;}
break;
case 19:
#line 49 "src/calc.y"
{yyval = ( yyvsp[-3] != yyvsp[0]) ? 1 : 0;}
break;
case 20:
#line 52 "src/calc.y"
{yyval = yyvsp[0];}
break;
case 21:
#line 53 "src/calc.y"
{yyval = yyvsp[-3] + yyvsp[0];}
break;
case 22:
#line 54 "src/calc.y"
{yyval = yyvsp[-3] - yyvsp[0];}
break;
case 23:
#line 55 "src/calc.y"
{yyval = - yyvsp[0];}
break;
case 24:
#line 56 "src/calc.y"
{yyval =   yyvsp[0];}
break;
case 25:
#line 59 "src/calc.y"
{yyval = yyvsp[0];}
break;
case 26:
#line 60 "src/calc.y"
{yyval = yyvsp[-3] * yyvsp[0];}
break;
case 27:
#line 61 "src/calc.y"
{yyval = yyvsp[-3] / yyvsp[0];}
break;
case 28:
#line 64 "src/calc.y"
{yyval = yyvsp[0];}
break;
case 29:
#line 65 "src/calc.y"
{yyval = yyvsp[-1] * yyvsp[0];}
break;
case 30:
#line 66 "src/calc.y"
{yyval = yyvsp[0];}
break;
case 31:
#line 67 "src/calc.y"
{yyval = yyvsp[-1];}
break;
case 32:
#line 68 "src/calc.y"
{yyval = pow(yyvsp[-3], yyvsp[0]);}
break;
case 33:
#line 69 "src/calc.y"
{yyval = (yyvsp[0] == 0) ? 1 : 0;}
break;
case 35:
#line 73 "src/calc.y"
{yyval = sqrt(yyvsp[-1]);}
break;
case 36:
#line 74 "src/calc.y"
{yyval = exp(yyvsp[-1]);}
break;
case 37:
#line 75 "src/calc.y"
{yyval = log(yyvsp[-1]);}
break;
case 38:
#line 76 "src/calc.y"
{yyval = log10(yyvsp[-1]);}
break;
case 39:
#line 77 "src/calc.y"
{yyval = sin(yyvsp[-1]);}
break;
case 40:
#line 78 "src/calc.y"
{yyval = cos(yyvsp[-1]);}
break;
case 41:
#line 79 "src/calc.y"
{yyval = tan(yyvsp[-1]);}
break;
case 42:
#line 80 "src/calc.y"
{yyval = sinh(yyvsp[-1]);}
break;
case 43:
#line 81 "src/calc.y"
{yyval = cosh(yyvsp[-1]);}
break;
case 44:
#line 82 "src/calc.y"
{yyval = tanh(yyvsp[-1]);}
break;
case 45:
#line 83 "src/calc.y"
{yyval = atan(yyvsp[-1]);}
break;
case 46:
#line 84 "src/calc.y"
{yyval = atan2(yyvsp[-3], yyvsp[-1]);}
break;
#line 693 "src/calc.c"
    }
    yyssp -= yym;
    yystate = *yyssp;
    yyvsp -= yym;
    yym = yylhs[yyn];
    if (yystate == 0 && yym == 0)
    {
#if YYDEBUG
        if (yydebug)
            printf("%sdebug: after reduction, shifting from state 0 to\
 state %d\n", YYPREFIX, YYFINAL);
#endif
        yystate = YYFINAL;
        *++yyssp = YYFINAL;
        *++yyvsp = yyval;
        if (yychar < 0)
        {
            if ((yychar = yylex()) < 0) yychar = 0;
#if YYDEBUG
            if (yydebug)
            {
                yys = 0;
                if (yychar <= YYMAXTOKEN) yys = yyname[yychar];
                if (!yys) yys = "illegal-symbol";
                printf("%sdebug: state %d, reading %d (%s)\n",
                        YYPREFIX, YYFINAL, yychar, yys);
            }
#endif
        }
        if (yychar == 0) goto yyaccept;
        goto yyloop;
    }
    if ((yyn = yygindex[yym]) && (yyn += yystate) >= 0 &&
            yyn <= YYTABLESIZE && yycheck[yyn] == yystate)
        yystate = yytable[yyn];
    else
        yystate = yydgoto[yym];
#if YYDEBUG
    if (yydebug)
        printf("%sdebug: after reduction, shifting from state %d \
to state %d\n", YYPREFIX, *yyssp, yystate);
#endif
    if (yyssp >= yysslim && yygrowstack())
    {
        goto yyoverflow;
    }
    *++yyssp = yystate;
    *++yyvsp = yyval;
    goto yyloop;
yyoverflow:
    yyerror("yacc stack overflow");
yyabort:
    return (1);
yyaccept:
    return (0);
}
