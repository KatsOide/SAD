c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      real*8 EPS
      parameter (EPS = 1.0d-7)
      character*1  LCURL,RCURL,LPAR,RPAR,COMMA,SEMIC
      parameter (LCURL='{',RCURL='}',LPAR='(',RPAR=')',COMMA=',')
      parameter (SEMIC=';')
      character*1  EQCHR,MINUS,PLUS,STAR
      parameter (EQCHR='=',MINUS='-',PLUS='+',STAR='*')
