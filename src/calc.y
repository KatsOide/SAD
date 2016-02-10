/* grammer rule for rdexpr routine for SAD1 */
/*
 * Action routines are translated from Fortran to C		by A.Morita(2008/03/20)
 * Grammer is Modified to accept `;', `)' and ID termination	by A.Morita(2008/03/20)
 */
%{
#include <sim/sad_f2c.h>
#include <math.h>
#define	YYMAXDEPTH	256
#define	YYSTYPE		real8
extern YYSTYPE yyval_cb;	/* yyval copy back*/
int  yylex(void);
void yyerror(const char *);
void yyunget(void);
%}
%token NUM ID UNIT EOL INVALID
%token FNSQRT FNEXP FNLN FNLOG
%token FNSIN FNCOS FNTAN FNATAN FNATAN2
%token FNSINH FNCOSH FNTANH 
%left '-' '+'
%left '*' '/' 
%left NEG
%right '^'
%right '~'
%start expr
%%
end_of_expr: ';' | ')' | ID;
opt_eol: /* empty */ | EOL;
expr:
UNIT			{yyerror("Unexpected UNIT");}
| or_expr end_of_expr	{$$ = $1; yyunget(); yyval_cb = $$; YYACCEPT;}
| or_expr		{$$ = $1; yyunget(); yyval_cb = $$;}
;
or_expr:
  and_expr
| or_expr '|' and_expr		{$$ = ($1 != 0 || $3 != 0) ? 1 : 0;}
;
and_expr:
  cmp_expr
| and_expr '&' cmp_expr		{$$ = ($1 != 0 && $3 != 0) ? 1 : 0;}
;
cmp_expr :
  a_expr
| cmp_expr '<' a_expr		{$$ = ( $1 <  $3) ? 1 : 0;}
| cmp_expr '>' a_expr		{$$ = ( $1 >  $3) ? 1 : 0;}
| cmp_expr '=' a_expr		{$$ = ( $1 == $3) ? 1 : 0;}
| cmp_expr '<' '=' a_expr	{$$ = ( $1 <= $4) ? 1 : 0;}
| cmp_expr '>' '=' a_expr	{$$ = ( $1 >= $4) ? 1 : 0;}
| cmp_expr '<' '>' a_expr	{$$ = ( $1 != $4) ? 1 : 0;}
;
a_expr:
  term				{$$ = $1;}
| a_expr '+' opt_eol term	{$$ = $1 + $4;}
| a_expr '-' opt_eol term	{$$ = $1 - $4;}
| '-' term      %prec NEG	{$$ = - $2;}
| '+' term      %prec NEG	{$$ =   $2;}
;
term:
  factor			{$$ = $1;}
| term '*' opt_eol factor	{$$ = $1 * $4;}
| term '/' opt_eol factor	{$$ = $1 / $4;}
;
factor:
  NUM				{$$ = $1;}
| NUM UNIT			{$$ = $1 * $2;}
| ID				{$$ = $1;}
| '(' or_expr ')'		{$$ = $2;}
| factor '^' opt_eol factor	{$$ = pow($1, $4);}
| '~' factor			{$$ = ($2 == 0) ? 1 : 0;}
| func_call
;
func_call:
  FNSQRT  '(' or_expr ')'		{$$ = sqrt($3);}
| FNEXP   '(' or_expr ')'		{$$ = exp($3);}
| FNLN    '(' or_expr ')'		{$$ = log($3);}
| FNLOG   '(' or_expr ')'		{$$ = log10($3);}
| FNSIN   '(' or_expr ')'		{$$ = sin($3);}
| FNCOS   '(' or_expr ')'		{$$ = cos($3);}
| FNTAN   '(' or_expr ')'		{$$ = tan($3);}
| FNSINH  '(' or_expr ')'		{$$ = sinh($3);}
| FNCOSH  '(' or_expr ')'		{$$ = cosh($3);}
| FNTANH  '(' or_expr ')'		{$$ = tanh($3);}
| FNATAN  '(' or_expr ')'		{$$ = atan($3);}
| FNATAN2 '(' or_expr ',' or_expr  ')'	{$$ = atan2($3, $5);}
; 
%%
