
/*  A Bison parser, made from calc.y  */

#define	NUM	258
#define	ID	259
#define	UNIT	260
#define	EOL	261
#define	EOF	262
#define	FNSQRT	263
#define	FNEXP	264
#define	FNLN	265
#define	FNLOG	266
#define	FNSIN	267
#define	FNCOS	268
#define	FNTAN	269
#define	FNATAN	270
#define	FNATAN2	271
#define	FNSINH	272
#define	FNCOSH	273
#define	FNTANH	274
#define	NEG	275

#line 3 "calc.y"

# define YYSTYPE REAL*8
# define YYMAXDEPTH 256 
#ifndef NULL
  #define NULL 0L
#endif

#ifndef YYLTYPE
typedef
  struct yyltype
    {
      int timestamp;
      int first_line;
      int first_column;
      int last_line;
      int last_column;
      char *text;
   }
  yyltype;

#define YYLTYPE yyltype
#endif

#ifndef YYSTYPE
#define YYSTYPE int
#endif
#include <stdio.h>

#ifndef __STDC__
#define const
#endif



#define	YYFINAL		104
#define	YYFLAG		-32768
#define	YYNTBASE	35

#define YYTRANSLATE(x) ((unsigned)(x) <= 275 ? yytranslate[x] : 44)

static const char yytranslate[] = {     0,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,    28,     2,    32,
    33,    22,    21,    34,    20,     2,    23,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,    29,
    31,    30,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,    25,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,    27,     2,    26,     2,     2,     2,     2,
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
     2,     2,     2,     2,     2,     1,     2,     3,     4,     5,
     6,     7,     8,     9,    10,    11,    12,    13,    14,    15,
    16,    17,    18,    19,    24
};

static const short yyprhs[] = {     0,
     0,     1,     3,     5,     7,    11,    13,    17,    19,    23,
    27,    31,    36,    41,    46,    48,    53,    58,    61,    64,
    66,    71,    76,    78,    81,    83,    87,    92,    95,    97,
   102,   107,   112,   117,   122,   127,   132,   137,   142,   147,
   152
};

static const short yyrhs[] = {    -1,
     6,     0,    37,     0,    38,     0,    37,    27,    38,     0,
    39,     0,    38,    28,    39,     0,    40,     0,    39,    29,
    40,     0,    39,    30,    40,     0,    39,    31,    40,     0,
    39,    29,    31,    40,     0,    39,    30,    31,    40,     0,
    39,    29,    30,    40,     0,    41,     0,    40,    21,    35,
    41,     0,    40,    20,    35,    41,     0,    20,    41,     0,
    21,    41,     0,    42,     0,    41,    22,    35,    42,     0,
    41,    23,    35,    42,     0,     3,     0,     3,     5,     0,
     4,     0,    32,    37,    33,     0,    42,    25,    35,    42,
     0,    26,    42,     0,    43,     0,     8,    32,    37,    33,
     0,     9,    32,    37,    33,     0,    10,    32,    37,    33,
     0,    11,    32,    37,    33,     0,    12,    32,    37,    33,
     0,    13,    32,    37,    33,     0,    14,    32,    37,    33,
     0,    17,    32,    37,    33,     0,    18,    32,    37,    33,
     0,    19,    32,    37,    33,     0,    15,    32,    37,    33,
     0,    16,    32,    37,    34,    37,    33,     0
};

#if YYDEBUG != 0
static const short yyrline[] = { 0,
    21,    21,    23,    26,    27,    30,    31,    34,    35,    36,
    37,    38,    39,    40,    43,    44,    46,    48,    50,    53,
    54,    56,    59,    60,    61,    62,    63,    64,    65,    67,
    70,    72,    74,    76,    78,    80,    82,    84,    86,    88,
    90
};

static const char * const yytname[] = {   "$",
"error","$illegal.","NUM","ID","UNIT","EOL","EOF","FNSQRT","FNEXP","FNLN",
"FNLOG","FNSIN","FNCOS","FNTAN","FNATAN","FNATAN2","FNSINH","FNCOSH","FNTANH","'-'",
"'+'","'*'","'/'","NEG","'^'","'~'","'|'","'&'","'<'","'>'",
"'='","'('","')'","','","opt_eol","expr","or_expr","and_expr","cmp_expr","a_expr",
"term","factor","func_call",""
};
#endif

static const short yyr1[] = {     0,
    35,    35,    36,    37,    37,    38,    38,    39,    39,    39,
    39,    39,    39,    39,    40,    40,    40,    40,    40,    41,
    41,    41,    42,    42,    42,    42,    42,    42,    42,    43,
    43,    43,    43,    43,    43,    43,    43,    43,    43,    43,
    43
};

static const short yyr2[] = {     0,
     0,     1,     1,     1,     3,     1,     3,     1,     3,     3,
     3,     4,     4,     4,     1,     4,     4,     2,     2,     1,
     4,     4,     1,     2,     1,     3,     4,     2,     1,     4,
     4,     4,     4,     4,     4,     4,     4,     4,     4,     4,
     6
};

static const short yydefact[] = {     0,
    23,    25,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     3,     4,
     6,     8,    15,    20,    29,    24,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,    18,    19,
    28,     0,     0,     0,     0,     0,     0,     1,     1,     1,
     1,     1,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,    26,     5,     7,     0,     0,     9,
     0,    10,    11,     2,     0,     0,     0,     0,     0,    30,
    31,    32,    33,    34,    35,    36,    40,     0,    37,    38,
    39,    14,    12,    13,    17,    16,    21,    22,    27,     0,
    41,     0,     0,     0
};

static const short yydefgoto[] = {    75,
   102,    19,    20,    21,    22,    23,    24,    25
};

static const short yypact[] = {    88,
    -2,-32768,     9,    23,    36,    37,    53,    54,    56,    58,
    83,   117,   118,   126,   113,   113,   113,    88,    29,    17,
    -6,   -14,    20,    33,-32768,-32768,    88,    88,    88,    88,
    88,    88,    88,    88,    88,    88,    88,    88,    20,    20,
-32768,   -25,    88,    88,    18,    63,    88,   112,   112,   112,
   112,   112,    19,    24,    26,    27,    60,    86,   108,   110,
    13,   111,   115,   119,-32768,    17,    -6,    88,    88,   -14,
    88,   -14,   -14,-32768,   113,   113,   113,   113,   113,-32768,
-32768,-32768,-32768,-32768,-32768,-32768,-32768,    88,-32768,-32768,
-32768,   -14,   -14,   -14,    20,    20,    33,    33,    33,   120,
-32768,   140,   151,-32768
};

static const short yypgoto[] = {   105,
-32768,   -18,   116,   121,    65,   -11,   -16,-32768
};


#define	YYLAST		165


static const short yytable[] = {    42,
    41,    43,    26,    39,    40,    48,    49,    65,    53,    54,
    55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
     1,     2,    45,    46,    47,     3,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    43,
    27,    50,    51,    17,    44,    43,    88,    68,    69,    18,
    43,    80,    43,    43,    28,    43,    81,    52,    82,    83,
    97,    98,    99,    95,    96,     1,     2,    29,    30,   100,
     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
    13,    14,    15,    16,    31,    32,    43,    33,    17,    34,
     1,     2,    84,    71,    18,     3,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    70,
    72,    73,    43,    17,    35,     1,     2,    74,    85,    18,
     3,     4,     5,     6,     7,     8,     9,    10,    11,    12,
    13,    14,    92,    93,    43,    94,    43,    43,    17,   103,
    86,    43,    87,    89,    18,    43,    43,    90,    36,    37,
   104,    91,   101,    76,    77,    78,    79,    38,    66,     0,
     0,     0,     0,     0,    67
};

static const short yycheck[] = {    18,
    17,    27,     5,    15,    16,    20,    21,    33,    27,    28,
    29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
     3,     4,    29,    30,    31,     8,     9,    10,    11,    12,
    13,    14,    15,    16,    17,    18,    19,    20,    21,    27,
    32,    22,    23,    26,    28,    27,    34,    30,    31,    32,
    27,    33,    27,    27,    32,    27,    33,    25,    33,    33,
    77,    78,    79,    75,    76,     3,     4,    32,    32,    88,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    20,    21,    32,    32,    27,    32,    26,    32,
     3,     4,    33,    31,    32,     8,     9,    10,    11,    12,
    13,    14,    15,    16,    17,    18,    19,    20,    21,    45,
    46,    47,    27,    26,    32,     3,     4,     6,    33,    32,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,    19,    68,    69,    27,    71,    27,    27,    26,     0,
    33,    27,    33,    33,    32,    27,    27,    33,    32,    32,
     0,    33,    33,    49,    50,    51,    52,    32,    43,    -1,
    -1,    -1,    -1,    -1,    44
};
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/local/lib/bison.simple"

/* Skeleton output parser for bison,
   Copyright (C) 1984, 1989, 1990 Bob Corbett and Richard Stallman

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.  */


#ifndef alloca
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* Not GNU C.  */
#if (!defined (__STDC__) && defined (sparc)) || defined (__sparc__)
#include <alloca.h>
#else /* Not sparc */
#ifdef MSDOS
#include <malloc.h>
#endif /* MSDOS */
#endif /* Not sparc.  */
#endif /* Not GNU C.  */
#endif /* alloca not defined.  */

/* This is the parser code that is written into each bison parser
  when the %semantic_parser declaration is not specified in the grammar.
  It was written by Richard Stallman by simplifying the hairy parser
  used when %semantic_parser is specified.  */

/* Note: there must be only one dollar sign in this file.
   It is replaced by the list of actions, each action
   as one case of the switch.  */

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		-2
#define YYEOF		0
#define YYACCEPT	return(0)
#define YYABORT 	return(1)
#define YYERROR		goto yyerrlab1
/* Like YYERROR except do call yyerror.
   This remains here temporarily to ease the
   transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL		goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(token, value) \
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    { yychar = (token), yylval = (value);			\
      yychar1 = YYTRANSLATE (yychar);				\
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    { yyerror ("syntax error: cannot back up"); YYERROR; }	\
while (0)

#define YYTERROR	1
#define YYERRCODE	256

#ifndef YYPURE
#define YYLEX		yylex()
#endif

#ifdef YYPURE
#ifdef YYLSP_NEEDED
#define YYLEX		yylex(&yylval, &yylloc)
#else
#define YYLEX		yylex(&yylval)
#endif
#endif

/* If nonreentrant, generate the variables here */

#ifndef YYPURE

int	yychar;			/*  the lookahead symbol		*/
YYSTYPE	yylval;			/*  the semantic value of the		*/
				/*  lookahead symbol			*/

#ifdef YYLSP_NEEDED
YYLTYPE yylloc;			/*  location data for the lookahead	*/
				/*  symbol				*/
#endif

int yynerrs;			/*  number of parse errors so far       */
#endif  /* not YYPURE */

#if YYDEBUG != 0
int yydebug;			/*  nonzero means print parse trace	*/
/* Since this is uninitialized, it does not stop multiple parsers
   from coexisting.  */
#endif

/*  YYINITDEPTH indicates the initial size of the parser's stacks	*/

#ifndef	YYINITDEPTH
#define YYINITDEPTH 200
#endif

/*  YYMAXDEPTH is the maximum size the stacks can grow to
    (effective only if the built-in stack extension method is used).  */

#if YYMAXDEPTH == 0
#undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
#define YYMAXDEPTH 10000
#endif

#ifndef __cplusplus

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_bcopy (from, to, count)
     char *from;
     char *to;
     int count;
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#else /* __cplusplus */

/* This is the most reliable way to avoid incompatibilities
   in available built-in functions on various systems.  */
static void
__yy_bcopy (char *from, char *to, int count)
{
  register char *f = from;
  register char *t = to;
  register int i = count;

  while (i-- > 0)
    *t++ = *f++;
}

#endif

#line 160 "/usr/local/lib/bison.simple"
int
yyparse()
{
  register int yystate;
  register int yyn;
  register short *yyssp;
  register YYSTYPE *yyvsp;
  int yyerrstatus;	/*  number of tokens to shift before error messages enabled */
  int yychar1;		/*  lookahead token as an internal (translated) token number */

  short	yyssa[YYINITDEPTH];	/*  the state stack			*/
  YYSTYPE yyvsa[YYINITDEPTH];	/*  the semantic value stack		*/

  short *yyss = yyssa;		/*  refer to the stacks thru separate pointers */
  YYSTYPE *yyvs = yyvsa;	/*  to allow yyoverflow to reallocate them elsewhere */

#ifdef YYLSP_NEEDED
  YYLTYPE *yyls = yylsa;
  YYLTYPE *yylsp;
  YYLTYPE yylsa[YYINITDEPTH];	/*  the location stack			*/

#define YYPOPSTACK   (yyvsp--, yysp--, yylsp--)
#else
#define YYPOPSTACK   (yyvsp--, yysp--)
#endif

  int yystacksize = YYINITDEPTH;

#ifdef YYPURE
  int yychar;
  YYSTYPE yylval;
  int yynerrs;
#ifdef YYLSP_NEEDED
  YYLTYPE yylloc;
#endif
#endif

  YYSTYPE yyval;		/*  the variable used to return		*/
				/*  semantic values from the action	*/
				/*  routines				*/

  int yylen;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Starting parse\n");
#endif

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.  */

  yyssp = yyss - 1;
  yyvsp = yyvs;
#ifdef YYLSP_NEEDED
  yylsp = yyls;
#endif

/* Push a new state, which is found in  yystate  .  */
/* In all cases, when you get here, the value and location stacks
   have just been pushed. so pushing a state here evens the stacks.  */
yynewstate:

  *++yyssp = yystate;

  if (yyssp >= yyss + yystacksize - 1)
    {
      /* Give user a chance to reallocate the stack */
      /* Use copies of these so that the &'s don't force the real ones into memory. */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;
#ifdef YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
#endif

      /* Get the current used size of the three stacks, in elements.  */
      int size = yyssp - yyss + 1;

#ifdef yyoverflow
      /* Each stack pointer address is followed by the size of
	 the data in use in that stack, in bytes.  */
      yyoverflow("parser stack overflow",
		 &yyss1, size * sizeof (*yyssp),
		 &yyvs1, size * sizeof (*yyvsp),
#ifdef YYLSP_NEEDED
		 &yyls1, size * sizeof (*yylsp),
#endif
		 &yystacksize);

      yyss = yyss1; yyvs = yyvs1;
#ifdef YYLSP_NEEDED
      yyls = yyls1;
#endif
#else /* no yyoverflow */
      /* Extend the stack our own way.  */
      if (yystacksize >= YYMAXDEPTH)
	{
	  yyerror("parser stack overflow");
	  return 2;
	}
      yystacksize *= 2;
      if (yystacksize > YYMAXDEPTH)
	yystacksize = YYMAXDEPTH;
      yyss = (short *) alloca (yystacksize * sizeof (*yyssp));
      __yy_bcopy ((char *)yyss1, (char *)yyss, size * sizeof (*yyssp));
      yyvs = (YYSTYPE *) alloca (yystacksize * sizeof (*yyvsp));
      __yy_bcopy ((char *)yyvs1, (char *)yyvs, size * sizeof (*yyvsp));
#ifdef YYLSP_NEEDED
      yyls = (YYLTYPE *) alloca (yystacksize * sizeof (*yylsp));
      __yy_bcopy ((char *)yyls1, (char *)yyls, size * sizeof (*yylsp));
#endif
#endif /* no yyoverflow */

      yyssp = yyss + size - 1;
      yyvsp = yyvs + size - 1;
#ifdef YYLSP_NEEDED
      yylsp = yyls + size - 1;
#endif

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Stack size increased to %d\n", yystacksize);
#endif

      if (yyssp >= yyss + yystacksize - 1)
	YYABORT;
    }

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Entering state %d\n", yystate);
#endif

 yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

  /* First try to decide what to do without reference to lookahead token.  */

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* yychar is either YYEMPTY or YYEOF
     or a valid token in external form.  */

  if (yychar == YYEMPTY)
    {
#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Reading a token: ");
#endif
      yychar = YYLEX;
    }

  /* Convert token to internal form (in yychar1) for indexing tables with */

  if (yychar <= 0)		/* This means end of input. */
    {
      yychar1 = 0;
      yychar = YYEOF;		/* Don't call YYLEX any more */

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Now at end of input.\n");
#endif
    }
  else
    {
      yychar1 = YYTRANSLATE(yychar);

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Next token is %d (%s)\n", yychar, yytname[yychar1]);
#endif
    }

  yyn += yychar1;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
    goto yydefault;

  yyn = yytable[yyn];

  /* yyn is what to do for this token type in this state.
     Negative => reduce, -yyn is rule number.
     Positive => shift, yyn is new state.
       New state is final state => don't bother to shift,
       just return success.
     0, or most negative number => error.  */

  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrlab;

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Shift the lookahead token.  */

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting token %d (%s), ", yychar, yytname[yychar1]);
#endif

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  /* count tokens shifted since error; after three, turn off error status.  */
  if (yyerrstatus) yyerrstatus--;

  yystate = yyn;
  goto yynewstate;

/* Do the default action for the current state.  */
yydefault:

  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;

/* Do a reduction.  yyn is the number of a rule to reduce with.  */
yyreduce:
  yylen = yyr2[yyn];
  yyval = yyvsp[1-yylen]; /* implement default value of the action */

#if YYDEBUG != 0
  if (yydebug)
    {
      int i;

      fprintf (stderr, "Reducing via rule %d (line %d), ",
	       yyn, yyrline[yyn]);

      /* Print the symboles being reduced, and their result.  */
      for (i = yyprhs[yyn]; yyrhs[i] > 0; i++)
	fprintf (stderr, "%s ", yytname[yyrhs[i]]);
      fprintf (stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif


  switch (yyn) {

case 3:
#line 24 "calc.y"
{yyval=yyvsp[0]; call yyunget();;
    break;}
case 5:
#line 28 "calc.y"
{yyval=0d0; if(yyvsp[-2] .ne. 0 .or. yyvsp[0] .ne. 0)  yyval=1d0;;
    break;}
case 7:
#line 32 "calc.y"
{yyval=0d0; if(yyvsp[-2] .ne. 0 .and. yyvsp[0] .ne. 0)  yyval=1d0;;
    break;}
case 9:
#line 36 "calc.y"
{yyval=0d0; if(yyvsp[-2] .lt. yyvsp[0])  yyval=1d0;;
    break;}
case 10:
#line 37 "calc.y"
{yyval=0d0; if(yyvsp[-2] .gt. yyvsp[0])  yyval=1d0;;
    break;}
case 11:
#line 38 "calc.y"
{yyval=0d0; if(yyvsp[-2] .eq. yyvsp[0])  yyval=1d0;;
    break;}
case 12:
#line 39 "calc.y"
{yyval=0d0; if(yyvsp[-3] .le. yyvsp[0])  yyval=1d0;;
    break;}
case 13:
#line 40 "calc.y"
{yyval=0d0; if(yyvsp[-3] .ge. yyvsp[0])  yyval=1d0;;
    break;}
case 14:
#line 41 "calc.y"
{yyval=0d0; if(yyvsp[-3] .ne. yyvsp[0])  yyval=1d0;;
    break;}
case 15:
#line 43 "calc.y"
{yyval=yyvsp[0];;
    break;}
case 16:
#line 45 "calc.y"
{yyval=yyvsp[-3] + yyvsp[0];;
    break;}
case 17:
#line 47 "calc.y"
{yyval=yyvsp[-3] - yyvsp[0];;
    break;}
case 18:
#line 49 "calc.y"
{yyval= - yyvsp[0];;
    break;}
case 19:
#line 51 "calc.y"
{yyval= yyvsp[0];;
    break;}
case 20:
#line 53 "calc.y"
{yyval=yyvsp[0];;
    break;}
case 21:
#line 55 "calc.y"
{yyval=yyvsp[-3] * yyvsp[0];;
    break;}
case 22:
#line 57 "calc.y"
{yyval=yyvsp[-3] / yyvsp[0];;
    break;}
case 23:
#line 59 "calc.y"
{yyval=yyvsp[0];;
    break;}
case 24:
#line 60 "calc.y"
{yyval=yyvsp[-1] * yyvsp[0];;
    break;}
case 25:
#line 61 "calc.y"
{yyval=yyvsp[0];;
    break;}
case 26:
#line 62 "calc.y"
{yyval=yyvsp[-1];;
    break;}
case 27:
#line 63 "calc.y"
{yyval=yyvsp[-3] ** yyvsp[0];;
    break;}
case 28:
#line 64 "calc.y"
{ yyval=0d0; if (yyvsp[0] .eq. 0d0) yyval=1d0;;
    break;}
case 30:
#line 69 "calc.y"
{yyval=sqrt(yyvsp[-1]);;
    break;}
case 31:
#line 71 "calc.y"
{yyval=exp(yyvsp[-1]);;
    break;}
case 32:
#line 73 "calc.y"
{yyval=log(yyvsp[-1]);;
    break;}
case 33:
#line 75 "calc.y"
{yyval=log10(yyvsp[-1]);;
    break;}
case 34:
#line 77 "calc.y"
{yyval=sin(yyvsp[-1]);;
    break;}
case 35:
#line 79 "calc.y"
{yyval=cos(yyvsp[-1]);;
    break;}
case 36:
#line 81 "calc.y"
{yyval=tan(yyvsp[-1]);;
    break;}
case 37:
#line 83 "calc.y"
{yyval=sinh(yyvsp[-1]);;
    break;}
case 38:
#line 85 "calc.y"
{yyval=cosh(yyvsp[-1]);;
    break;}
case 39:
#line 87 "calc.y"
{yyval=tanh(yyvsp[-1]);;
    break;}
case 40:
#line 89 "calc.y"
{yyval=datan(yyvsp[-1]);;
    break;}
case 41:
#line 91 "calc.y"
{yyval=datan2(yyvsp[-3],yyvsp[-1]);;
    break;}
}
   /* the action file gets copied in in place of this dollarsign */
#line 423 "/usr/local/lib/bison.simple"

  yyvsp -= yylen;
  yyssp -= yylen;
#ifdef YYLSP_NEEDED
  yylsp -= yylen;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

  *++yyvsp = yyval;

#ifdef YYLSP_NEEDED
  yylsp++;
  if (yylen == 0)
    {
      yylsp->first_line = yylloc.first_line;
      yylsp->first_column = yylloc.first_column;
      yylsp->last_line = (yylsp-1)->last_line;
      yylsp->last_column = (yylsp-1)->last_column;
      yylsp->text = 0;
    }
  else
    {
      yylsp->last_line = (yylsp+yylen-1)->last_line;
      yylsp->last_column = (yylsp+yylen-1)->last_column;
    }
#endif

  /* Now "shift" the result of the reduction.
     Determine what state that goes to,
     based on the state we popped back to
     and the rule number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
  if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTBASE];

  goto yynewstate;

yyerrlab:   /* here on detecting error */

  if (! yyerrstatus)
    /* If not already recovering from an error, report this error.  */
    {
      ++yynerrs;

#ifdef YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (yyn > YYFLAG && yyn < YYLAST)
	{
	  int size = 0;
	  char *msg;
	  int x, count;

	  count = 0;
	  for (x = 0; x < (sizeof(yytname) / sizeof(char *)); x++)
	    if (yycheck[x + yyn] == x)
	      size += strlen(yytname[x]) + 15, count++;
	  msg = (char *) xmalloc(size + 15);
	  strcpy(msg, "parse error");

	  if (count < 5)
	    {
	      count = 0;
	      for (x = 0; x < (sizeof(yytname) / sizeof(char *)); x++)
		if (yycheck[x + yyn] == x)
		  {
		    strcat(msg, count == 0 ? ", expecting `" : " or `");
		    strcat(msg, yytname[x]);
		    strcat(msg, "'");
		    count++;
		  }
	    }
	  yyerror(msg);
	  free(msg);
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror("parse error");
    }

yyerrlab1:   /* here on error raised explicitly by an action */

  if (yyerrstatus == 3)
    {
      /* if just tried and failed to reuse lookahead token after an error, discard it.  */

      /* return failure if at end of input */
      if (yychar == YYEOF)
	YYABORT;

#if YYDEBUG != 0
      if (yydebug)
	fprintf(stderr, "Discarding token %d (%s).\n", yychar, yytname[yychar1]);
#endif

      yychar = YYEMPTY;
    }

  /* Else will try to reuse lookahead token
     after shifting the error token.  */

  yyerrstatus = 3;		/* Each real token shifted decrements this */

  goto yyerrhandle;

yyerrdefault:  /* current state does not do anything special for the error token. */

#if 0
  /* This is wrong; only states that explicitly want error tokens
     should shift them.  */
  yyn = yydefact[yystate];  /* If its default is to accept any token, ok.  Otherwise pop it.*/
  if (yyn) goto yydefault;
#endif

yyerrpop:   /* pop the current state because it cannot handle the error token */

  if (yyssp == yyss) YYABORT;
  yyvsp--;
  yystate = *--yyssp;
#ifdef YYLSP_NEEDED
  yylsp--;
#endif

#if YYDEBUG != 0
  if (yydebug)
    {
      short *ssp1 = yyss - 1;
      fprintf (stderr, "Error: state stack now");
      while (ssp1 != yyssp)
	fprintf (stderr, " %d", *++ssp1);
      fprintf (stderr, "\n");
    }
#endif

yyerrhandle:

  yyn = yypact[yystate];
  if (yyn == YYFLAG)
    goto yyerrdefault;

  yyn += YYTERROR;
  if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
    goto yyerrdefault;

  yyn = yytable[yyn];
  if (yyn < 0)
    {
      if (yyn == YYFLAG)
	goto yyerrpop;
      yyn = -yyn;
      goto yyreduce;
    }
  else if (yyn == 0)
    goto yyerrpop;

  if (yyn == YYFINAL)
    YYACCEPT;

#if YYDEBUG != 0
  if (yydebug)
    fprintf(stderr, "Shifting error token, ");
#endif

  *++yyvsp = yylval;
#ifdef YYLSP_NEEDED
  *++yylsp = yylloc;
#endif

  yystate = yyn;
  goto yynewstate;
}
#line 93 "calc.y"

