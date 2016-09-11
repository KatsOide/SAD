c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c$$$c
c$$$c     Update History:	
c$$$c     2016/01/12:	from PDG 2014 https://www.google.co.jp/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0ahUKEwjHw4unsKLKAhUEoJQKHaiPCLgQFggjMAE&url=http%3A%2F%2Fpdg.lbl.gov%2F2014%2Freviews%2Frpp2014-rev-phys-constants.pdf&usg=AFQjCNGXw6APvgkpoTdIRLeC-XV651ZbJg&sig2=-wjwSn5pfl6Juvok0hFqow
c$$$c     2008/03/27:	from CODATA 2006
c$$$c     			http://physics.nist.gov/cuu/Constants/Table/allascii.txt
c$$$c
c$$$c     2003/07/15:	from Particle Data Book
c$$$c
c$$$
c$$$c     Including math constant for permeability of vacuum
c$$$      include 'inc/MACMATH.inc'
c$$$
c$$$      real*8 cveloc,mu0,ep0
c$$$      real*8 plank,plankr,elemch,echarg,finest
c$$$      real*8 elmass,prmass,elradi,prradi
c$$$
c$$$c     Speed of light in vacuum:		c
c$$$c     Definition:	299 792 458. m/sec
c$$$      parameter (cveloc = 299792458d0)
c$$$
c$$$c     Permeability of vacuum:		\mu_0
c$$$c     Definition:	4 Pi x 10^-7 N A^-2
c$$$      parameter (mu0 = pi * 4.d-7)
c$$$
c$$$c     Permittivity of vacuum:		\epsilon_0 = 1 / (\mu_0 c^2)
c$$$c     Delivered from defintion
c$$$      parameter (ep0 = 1.d0 / (mu0 * cveloc**2))
c$$$
c$$$c     Planck's constant:		h
c$$$c     CODATA 2006:	6.626 068 96(33)   x 10^-34 Js
c$$$c     PDG 2014::	6.626 069 57(29)  x 10^-34 Js
c$$$      parameter (plank  = 6.62606957d-34)
c$$$
c$$$c     Dirac's constant:			hbar = h / (2 Pi)
c$$$c     PDB 2003:		1.054 571 596(82)  x 10^-34 Js
c$$$c     CODATA 2006:	1.054 571 628(53)  x 10^-34 Js
c$$$c     PDG 2014: 	1.054 571 726(47)  x 10^-34 Js
c$$$      parameter (plankr = 1.054571726d-34)
c$$$
c$$$c     Elementary charge:		e
c$$$c     PDB 2003:		1.602 176 462(63)  x 10^-19 C
c$$$c     CODATA 2006:	1.602 176 487(40)  x 10^-19 C
c$$$c     PDG 2014: 	1.602 176 565(35)  x 10^-19 C
c$$$      parameter (elemch = 1.602176565d-19)
c$$$
c$$$c     Electron charge in elementary charge unit
c$$$      parameter (echarg = 1.0d0)
c$$$
c$$$c     Fine-structure constant		\alpha = \mu_0 e^2 c / (2 h)
c$$$c     PDB 2003:		1 / 137.035 999 76(50)
c$$$c     CODATA 2006:	1 / 137.035 999 679(94)
c$$$c     PDG 2014: 	1 / 137.035 999 074(44)
c$$$      parameter (finest = 1.d0 / 137.035999074d0)
c$$$
c$$$c     Electron mass energy equivalent in eV:	m_e c^2 / e
c$$$c     PDB 2003:		.510 998 902(21) MeV
c$$$c     CODATA 2006:	.510 998 910(13) MeV
c$$$c     PDG 2014: 	.510 998 928(11) MeV
c$$$      parameter (elmass =   0.510998928d6)
c$$$
c$$$c     Proton mass energy equivalent in eV:	m_p c^2 / e
c$$$c     PDB 2003:		938. 271 998(38) MeV
c$$$c     CODATA 2006:	938. 272 013(23) MeV
c$$$c     PDG 2014: 	938. 272 046(21) MeV
c$$$      parameter (prmass = 938.272046d6)
c$$$
c$$$c     Classical electron radius:	r_e = e^2 / (4Pi \epsilon_0 m_e c^2)
c$$$c     					    = (\alpha hbar c) / (m_e c^2)
c$$$c     					    = (e^2 c^2 \mu_0 / 4Pi) / (m_e c^2)
c$$$c     					    = e c^2 * 10^-7 * (e / (m_e c^2))
c$$$c     PDB 2003:		2.817 940 285(??)  x 10^-15 m
c$$$c     CODATA 2006:	2.817 940 2894(58) x 10^-15 m
c$$$c     PDG 2014: 	2.817 940 3267(27) x 10^-15 m
c$$$c      parameter (elradi = finest * plankr * cveloc / (elemch * elmass))
c$$$c      parameter (elradi = elemch * cveloc**2 * 1.d-7 / elmass)
c$$$       parameter (elradi = 2.8179403267d-15)
c$$$
c$$$c     Classical proton radius:		r_p = e^2 / (4Pi \epsilon_0 m_p c^2)
c$$$c     					    = (\alpha hbar c) / (m_p c^2)
c$$$c     					    = (e^2 c^2 \mu_0 / 4Pi) / (m_p c^2)
c$$$c     					    = e c^2 * 10^-7 * (e / (m_p c^2))
c$$$c      parameter (prradi = finest * plankr * cveloc / (elemch * prmass))
c$$$       parameter (prradi = elemch * cveloc**2 * 1.d-7 / prmass)
