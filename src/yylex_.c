#include <sim/sad_api.h>
#include <sim/MACTTYP.h>
#include <sim/MACCODE.h>
#include <sim/MACCBK.h>
#include <sim/sad_memory.h>

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* external yacc parser symbol */
#define	YYSTYPE		real8
#include <calc.h>
extern YYSTYPE yylval;
extern int yyparse(void);

YYSTYPE yyval_cb;		/* yyval copy back */

/* function name table */
typedef struct {
  const char *name;
  int token;
} yyfunc_tbl_t;

static yyfunc_tbl_t yyfunc_tbl[] = {
  {"SIN",	FNSIN},
  {"COS",	FNCOS},
  {"TAN",	FNTAN},
  {"ATAN",	FNATAN},
  {"ATAN2",	FNATAN2},
  {"SINH",	FNSINH},
  {"COSH",	FNCOSH},
  {"TANH",	FNTANH},
  {"SQRT",	FNSQRT},
  {"EXP",	FNEXP},
  {"LOG",	FNLOG},
  {"LN",	FNLN},
  {NULL,	0}
};

/* Token reader state */
static struct {
  int type;
  int len;
  char token[1024];
} yylex_state;

/* Next token reader for yacc parser */
#ifdef DEBUG_YYLEX
static int debug_level = 0;

static int yylex0(void);

int yylex(void) {
  int yyresult;

  if(debug_level & 0x0010) {
    fprintf(stderr, "yylex: ");
    fprintf(stderr, !(yylex_state.len > 0) ? "  new " : "      ");
  }
  yyresult = yylex0();
  if(debug_level & 0x0010) {
    fprintf(stderr, "token type=%d [", yylex_state.type);
    for(int i = 0; i < abs(yylex_state.len); i++) fprintf(stderr, "%c", yylex_state.token[i]);
    fprintf(stderr, "] -> ");
    fprintf(stderr, " token=%d, yylval=%f\n", yyresult, yylval);
  }
  return yyresult;
}

static int yylex0(void) {
#else
int yylex(void) {
#endif
  integer4 len, type, ival, hash;
  real8 rval;
  char *token;

  token = yylex_state.token;
  type  = yylex_state.type;
  len   = yylex_state.len;
  if(!(len > 0)) {
    len = 0;
    gettok_(token, &len, &type, &rval, &ival, sizeof(yylex_state.token));
    yylex_state.type = type;
  }
  yylex_state.len  = -len; /* Mark token is already got */
  yylval = 0; /* Clear token value */
  switch(type) {
  case ttypID:
    hash = hsrchz_(token, len);
    switch(idtype(hash)) {
    case icGLR:
      yylval = rgetgl1_(token, len);
      return ID;

    case icGLI:
      yylval = igetgl_(token, len);
      return ID;

    default:
      capita_(token, len);
      hash = hsrch_(token, len);
      switch(idtype(hash)) {
      case icUNIT:
	yylval = rlist(idval(hash));
	return UNIT;

      default:
	for(yyfunc_tbl_t *entry = yyfunc_tbl; entry->name != NULL; entry++)
	  if(strlen(entry->name) == len && strncmp(entry->name, token, len) == 0)
	    return entry->token;
      }
    }
    return EOF;

  case ttypNM:
    yylval = atof_(token, &len, len);
    return NUM;

  case ttypDM:
    return token[0];

  case ttypEF:
    return EOF;

  default:
    return EOF;
  }
}

/* Parser support routines */
void yyerror(const char *message) {
  fprintf(stderr, "LALR(1) parser: %s\n", message);
}

void yyunget(void) {
#ifdef DEBUG_YYLEX
  if(debug_level & 0x0004) {
    fprintf(stderr, "yyunget: len=%d type=%d [", yylex_state.len, yylex_state.type);
    for(int i = 0; i < abs(yylex_state.len); i++) fprintf(stderr, "%c", yylex_state.token[i]);
    fprintf(stderr, "]\n");
  }
#endif
  if(yylex_state.len < 0)
    yylex_state.len = -yylex_state.len;
}

/* Fortran access APIs */
integer4 yypushtoken_(const_character token, integer4 *len, integer4 *type, ftnlen tlen) {
  if(*len > sizeof(yylex_state.token)) {
    fprintf(stderr, "yypushtoken(): Too large token(%d > %ubytes)\n", *len, (unsigned int)sizeof(yylex_state.token));
    abort();
  }

  yylex_state.type = *type;
  yylex_state.len  = *len;
  memcpy(yylex_state.token, token, *len);

#ifdef DEBUG_YYLEX
  char *debug_switch = getenv("DEBUG_YYLEX");

  if(debug_switch != NULL) {
    debug_level = strtol(debug_switch, NULL, 10);
  }

  if(debug_level & 0x0001) {
    fprintf(stderr, "yypushtoken: token type=%d [", *type);
    for(int i = 0; i < *len; i++) fprintf(stderr, "%c", token[i]);
    fprintf(stderr, "]\n");
  }
#endif
  return 0;
}

integer4 yypoptoken_(character token, integer4 *len, integer4 *type, ftnlen tlen) {
  int copy_len = abs(yylex_state.len);

  if(copy_len > tlen) {
    fprintf(stderr, "yypushtoken(): Too large token(%d > %ubytes)\n", copy_len, tlen);
    abort();
  }

  *type = yylex_state.type;
  *len  = copy_len;
  memcpy(token, yylex_state.token, copy_len);
  yylex_state.len = -copy_len;

#ifdef DEBUG_YYLEX
  if(debug_level & 0x0002) {
    fprintf(stderr, "yypoptoken:  token type=%d [", *type);
    for(int i = 0; i < *len; i++) fprintf(stderr, "%c", token[i]);
    fprintf(stderr, "]\n");
  }
#endif
  return 0;
}

integer4 yyparse_(integer4 *yyreturn_, real8 *yyval_) {
  int status = yyparse();
  *yyval_ = yyval_cb;
#ifdef DEBUG_YYLEX
  if(debug_level & 0x0008) {
    fprintf(stderr, "yyparse: yyval=%f\n", *yyval_);
  }
#endif
  return status;
}

/* End of File */
