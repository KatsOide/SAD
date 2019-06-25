#ifndef _MACMATH_H_
#define _MACMATH_H_

/* sim for inc/MACMATH.inc */
#include <sim/sad_f2c.h>

// m_* math constat is defined as ``const rael8''
// in order to refer as Fortran Real*8 value

// Napier's constant
//	\lim_{n\rightarrow\infty}(1 + \frac{1}{n})^n
const real8 m_e        = 2.7182818284590452354;		// Exp[1]

// Euler-Mascheroni constant
//	\lim_{n\rightarrow\infty}(\sum_{k=1}^{n}\frac{1}{k} - \ln{n})
const real8 m_euler    = 0.57721566490153286061;	// Euler's gamma

const real8 m_log2e    = 1.4426950408889634074;		// Log[2,  Exp[1]]

const real8 m_log10e   = 0.43429448190325182765;	// Log[10, Exp[1]]

const real8 m_ln2      = 0.69314718055994530942;	// Log[Exp[1],  2]

const real8 m_ln10     = 2.30258509299404568402;	// Log[Exp[1], 10]

const real8 m_pi       = 3.14159265358979323846;	//     Pi

const real8 m_2pi      = 6.28318530717958647693;	// 2 * Pi

const real8 m_4pi      = 12.5663706143591729539;	// 4 * Pi

const real8 m_pi_2     = 1.57079632679489661923;	//     Pi / 2

const real8 m_pi_4     = 0.78539816339744830962;	//     Pi / 4

const real8 m_1_pi     = 0.31830988618379067154;	// 1 / Pi

const real8 m_2_pi     = 0.63661977236758134308;	// 2 / Pi

const real8 m_4_pi     = 1.27323954473516268615;	// 4 / Pi

const real8 m_sqrtpi   = 1.77245385090551602730;	//     Sqrt[Pi]

const real8 m_2_sqrtpi = 1.12837916709551257390;	// 2 / Sqrt[Pi]

const real8 m_sqrt2    = 1.41421356237309504880;	//     Sqrt[2]

const real8 m_1_sqrt2  = 0.70710678118654752440;	//     Sqrt[1 / 2]

const real8 m_sqrt3    = 1.73205080756887729353;	//     Sqrt[3]

const real8 m_1_sqrt3  = 0.57735026918962576451;	//     Sqrt[1 / 3]

//     Standard alias for SAD sources
#define	napier	m_e
#define	euler	m_euler
#define	pi	m_pi
#define	pi2	m_2pi
#define	pi4	m_4pi
#define	hpi	m_pi_2
#define	qpi	m_pi_4

#endif /* _MACMATH_H_ */
