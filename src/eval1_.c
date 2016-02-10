#include <sim/sad_f2c.h>

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <sysexits.h>

real8 eval1_(const_character string, const integer4 *l,
	     const integer4 *i1, integer4 *m) {
  const int len = *l;
  int i = *i1 - 1;
  size_t pos_exp = 0;

  *m = 0;
  bool hexadecimal = false;
  int (*digitQ)(int) = isdigit;
  int i_begin = i;

  if(isdigit(string[i])) {
    if(string[i] == '0' &&
       i + 1 < len && (string[i + 1] == 'x' || string[i + 1] == 'X')) {
      if((i + 2 < len && isxdigit(string[i + 2])) ||
	 (i + 3 < len && string[i + 2] =='.' && isxdigit(string[i + 3]))) {
	hexadecimal = true;
	digitQ = isxdigit;
	i += 2;
	/* Hexdecimal */
      }
    }
    /* Decimal */
  } else if(string[i] == '.' &&
	    i + 1 < len && isdigit(string[i + 1])) {
    /* '.' prefixed Decimal */
  } else {
    return 0;	/* Error */
  }

  /* Check digits before `.' */
  while(i < len && digitQ(string[i])) { ++i; }

  /* Check `.' and tailing digits */
  if(i < len && string[i] == '.') {
    if(i + 1 < len && string[i + 1] != '.') {
      i += 1;
      while(i < len && digitQ(string[i])) { ++i; }
    } else if(i + 1 == len) {
      i += 1;
    }
  }

  /* Check exponent */
  if(i < len &&
     (hexadecimal
      ? (string[i] == 'p' || string[i] == 'P')
      : (string[i] == 'e' || string[i] == 'E' ||
	 string[i] == 'd' || string[i] == 'D'))) {
    pos_exp = i - i_begin;
    if(i + 1 < len && isdigit(string[i + 1])) {
      i += 2;
      while(i < len && isdigit(string[i])) { ++i; }
    } else if(i + 2 < len && (string[i + 1] == '+' || string[i + 1] == '-') &&
	      isdigit(string[i + 2])) {
      i += 3;
      while(i < len && isdigit(string[i])) { ++i; }
    }
  }

  *m = i - i_begin;
  if(*m > 0) {
    double r;
    char buffer[*m + 1], *endptr;

    memcpy(buffer, string + i_begin, *m);
    buffer[*m] = '\0';
    if(pos_exp > 0 &&
       (buffer[pos_exp] == 'd' || buffer[pos_exp] == 'D')) {
      /* Replace `d' exponent prefix with `e' */
      buffer[pos_exp] = 'e';
    }
    r = strtod(buffer, &endptr);
    if(endptr != buffer + *m) {
      fprintf(stderr, "eval1_: buffer=\"%s\" rest=\"%s\"\n",
	      buffer, endptr);
      exit(EX_SOFTWARE);
    }
    return r;
  } else {
    return 0;
  }
}

/* End of File */
