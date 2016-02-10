c     Don't confuse, Emacs. This is -*- fortran -*- mode!
c     Math constant
      real*8 m_e, m_euler,
     $     m_log2e, m_log10e, m_ln2, m_ln10,
     $     m_pi, m_2pi, m_4pi, m_pi_2, m_pi_4,
     $     m_1_pi, m_2_pi, m_4_pi,
     $     m_sqrtpi, m_2_sqrtpi, m_sqrt2, m_1_sqrt2,
     $     m_sqrt3, m_1_sqrt3

c     Napier's constant(Exp[1])
c     \lim_{n\rightarrow\infty}(1 + \frac{1}{n})^n
      parameter(m_e        = 2.7182818284590452354d0)

c     Euler-Mascheroni constant(Euler's gamma)
c     \lim_{n\rightarrow\infty}(\sum_{k=1}^{n}\frac{1}{k} - \ln{n})
      parameter(m_euler    = 0.57721566490153286061d0)

c     Log[2,  Exp[1]]
      parameter(m_log2e    = 1.4426950408889634074d0)

c     Log[10, Exp[1]]
      parameter(m_log10e   = 0.43429448190325182765d0)

c     Log[Exp[1],  2]
      parameter(m_ln2      = 0.69314718055994530942d0)

c     Log[Exp[1], 10]
      parameter(m_ln10     = 2.30258509299404568402d0)

c         Pi
      parameter(m_pi       = 3.14159265358979323846d0)

c     2 * Pi
      parameter(m_2pi      = 6.28318530717958647693d0)

c     4 * Pi
      parameter(m_4pi      = 12.5663706143591729539d0)

c         Pi / 2
      parameter(m_pi_2     = 1.57079632679489661923d0)

c         Pi / 4
      parameter(m_pi_4     = 0.78539816339744830962d0)

c     1 / Pi
      parameter(m_1_pi     = 0.31830988618379067154d0)

c     2 / Pi
      parameter(m_2_pi     = 0.63661977236758134308d0)

c     4 / Pi
      parameter(m_4_pi     = 1.27323954473516268615d0)

c         Sqrt[Pi]
      parameter(m_sqrtpi   = 1.77245385090551602730d0)

c     2 / Sqrt[Pi]
      parameter(m_2_sqrtpi = 1.12837916709551257390d0)

c     Sqrt[2]
      parameter(m_sqrt2    = 1.41421356237309504880d0)

c     Sqrt[1/2]
      parameter(m_1_sqrt2  = 0.70710678118654752440d0)

c     Sqrt[3]
      parameter(m_sqrt3    = 1.73205080756887729353d0)

c     Sqrt[1/3]
      parameter(m_1_sqrt3  = 0.57735026918962576451d0)

c     Standard alias for SAD sources
      real*8 napier, euler, pi, pi2, pi4, hpi, qpi
      parameter(napier = m_e, euler = m_euler)
      parameter(pi  = m_pi)
      parameter(pi2 = m_2pi,  pi4 = m_4pi)
      parameter(hpi = m_pi_2, qpi = m_pi_4)

c     End of File
