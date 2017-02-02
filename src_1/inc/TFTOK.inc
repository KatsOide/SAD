c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$      character*28 oper
c$$$      character*24 oper1
c$$$      character*15 oper2
c$$$      character*4 oper3
c$$$      character*11 opern
c$$$      parameter (oper=' +-*/=(){}[],;@#:^&<>|~.?"''')
c$$$      parameter (oper1='+-*/=(){}[],;@#:^&<>|~.?')
c$$$      parameter (oper2='=><&|/.+-:#@[*)')
c$$$      parameter (oper3='=>@.')
c$$$      parameter (opern='0123456789.')
c$$$      integer*4 levelop(0:255),ichorder(0:255)
c$$$      common /tok/ levelop,ichorder
