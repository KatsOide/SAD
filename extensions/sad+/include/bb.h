#ifndef IP_H
#define IP_H
#include <element.h>
#include <matrix.h>

class IP : public Element
{
   double x_angle,cosx,sinx,tanx;
   double N_particle;
   double sigx,sigy,sigz;
   double xi_x,xi_y;
   int n_slice;
   double *np_slice,*z_slice;
   matrix* Beam_envelope;
   matrix* Benv_headon;
   matrix* Benv_slice;
   double cod[6];
   double bcen_fac[5],tilt_angle,cos_t,sin_t;
   double strong_beam_energy,weak_beam_energy;
   double luminosity,dLum;
   char lmon;
public:
   IP(char*);
   IP(double*);
   ~IP(void);
   void IP_reset(double*);

   void print(void);

   double Luminosity(void) {return luminosity*1.e-4;}
   double dLuminosity(void) {return dLum*1.e-4;}
   int FindIP(void) {return 1;}
   void LuminosityMonitorStart(void) {lmon=1; luminosity=0.;}
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   void beam_set(Beam& BEAM) {weak_beam_energy=BEAM.Energy();}
};


#endif
