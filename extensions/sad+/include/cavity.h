
class Cavity : public Element
{
   int harm;
   double phi,freq,voltage;
   double w,sinp,cosp,VE;

public:
   Cavity(double* K);
   Cavity(char*);
   void print(void)
   { Element::print(); 
      cout<< ", voltage=" << voltage << " V, freq=" << freq 
	 << ", phi=" << phi << ")\n" << "      VE=" << VE << "\n";}

   void beam_set(Beam& Beam) {
      VE=voltage/Beam.Energy()*1.e-9; }
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   
};
