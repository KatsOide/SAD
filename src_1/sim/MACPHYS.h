#ifndef _MACPHYS_H_
#define _MACPHYS_H_

/* sim for inc/MACPHYS.inc */
#include <sim/sad_f2c.h>

// Update History:	
// 2008/03/27:	from CODATA 2006
//		http://physics.nist.gov/cuu/Constants/Table/allascii.txt
//
// 2003/07/15:	from Particle Data Book

// Local math constant(Pi) for definition of permeability of vacuum
#define	M_pi			3.14159265358979323846

// Physical constant macros

// Speed of light in vacuum:	c
// Definition:	299 792 458. m/sec
#define	P_SpeedOfLight		299792458.

// Permeability of vacuum:	\mu_0
// Definition:	4 Pi x 10^-7 N A^-2
#define	P_Permeability		(M_pi * 4e-7)

// Permittivity of vacuum:	\epsilon_0 = 1 / (\mu_0 c^2)
// Delivered from defintion
#define	P_Permittivity		(1 / (M_Permeability * P_SpeedOfLight * P_SpeedOfLight))

// Planck's constant:		h
// CODATA 2006:	6.626 068 96(33)   x 10^-34 Js
#define	P_PlanckConstant	6.62606896e-34;

// Dirac's constant:		hbar = h / (2 Pi)
// PDB 2003:	1.054 571 596(82)  x 10^-34 Js
// CODATA 2006:	1.054 571 628(53)  x 10^-34 Js
#define	P_DiracConstant		1.054571628e-34

// Elementary charge:		e
// PDB 2003:	1.602 176 462(63)  x 10^-19 C
// CODATA 2006:	1.602 176 487(40)  x 10^-19 C
#define	P_ElementaryCharge	1.602176487e-19

// Fine-structure constant	\alpha = \mu_0 e^2 c / (2 h)
// PDB 2003:	1 / 137.035 999 76(50)
// CODATA 2006:	1 / 137.035 999 679(94)
#define	P_FineStructureConstant	(1 / 137.035999679)

// Electron mass energy equivalent in eV:	m_e c^2 / e
// PDB 2003:	.510 998 902(21) MeV
// CODATA 2006:	.510 998 910(13) MeV
#define	P_ElectronMassEV	  0.510998910e6

// Proton mass energy equivalent in eV:		m_p c^2 / e
// PDB 2003:	938. 271 998(38) MeV
// CODATA 2006:	938. 272 013(23) MeV
#define	P_ProtonMassEV		938.272013e6

// Classical electron radius:	r_e = e^2 / (4Pi \epsilon_0 m_e c^2)
//				    = (\alpha hbar c) / (m_e c^2)
//				    = (e^2 c^2 \mu_0 / 4Pi) / (m_e c^2)
//				    = e c^2 * 10^-7 * (e / (m_e c^2))
// PDB 2003:	2.817 940 285(??)  x 10^-15 m
// CODATA 2006:	2.817 940 2894(58) x 10^-15 m
//#define	P_ElectronRadius	((P_FineStructureConstant * P_DiracConstant * P_SpeedOfLight) / (P_ElementaryCharge * P_ElectronMassEV))
//#define	P_ElectronRadius	(P_ElementaryCharge * P_SpeedOfLight * P_SpeedOfLight * 1e-7 / P_ElectronMassEV)
#define	P_ElectronRadius	2.8179402894e-15

// Classical proton radius:	r_p = e^2 / (4Pi \epsilon_0 m_p c^2)
//				    = (\alpha hbar c) / (m_p c^2)
//				    = (e^2 c^2 \mu_0 / 4Pi) / (m_p c^2)
//				    = e c^2 * 10^-7 * (e / (m_p c^2))
//#define	P_ProtonRadius		((P_FineStructureConstant * P_DiracConstant * P_SpeedOfLight) / (P_ElementaryCharge * P_ProtonMassEV))
#define	P_ProtonRadius		(P_ElementaryCharge * P_SpeedOfLight * P_SpeedOfLight * 1e-7 / P_ProtonMassEV)

// MACPHYS.inc compatibility variables
// These physical constat variables are defined as ``const rael8''
// in order to refer as Fortran Real*8 value

// c velocity
const real8 cveloc = P_SpeedOfLight;

// Planck constant over 2Pi
const real8 plank  = P_PlanckConstant;

// Planck constant over 2Pi
const real8 plankr = P_DiracConstant;

// Electron charge in elementary charge unit
const real8 echarg = 1.;

// Fine-structure constant
const real8 finest = P_FineStructureConstant;

// Electron mass in eV/c^2 unit
const real8 elmass = P_ElectronMassEV;

// Proton mass in eV/c^2 unit
const real8 prmass = P_ProtonMassEV;

// Classical radius of electron
const real8 elradi = P_ElectronRadius;

// Classical radius of proton
const real8 prradi = P_ProtonRadius;

#undef	M_pi

#endif /* _MACPHYS_H_ */
