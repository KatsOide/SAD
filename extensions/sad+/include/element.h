#ifndef ELEMENT_H
#define ELEMENT_H

#include <iostream>
using std::istream;
using std::ostream;
using std::cout;
#include <cstring>
using std::strcpy;

#include <phys_const.h>
#include <map_double.h>
#include <map_da.h>
#include <map_p_da.h>
#include <track.h>


#define n_data_bend 3
#define n_data_quad 1
#define N_slice_max 100
#define Bres_max 1024

double get_parm(char*,const char*);

class Beam
{
   enum p_type { ELECTRON=1, POSITRON, PROTON, ANTI_PROTON};
   p_type particle;
   double mass,charge;
   double energy;
   double emx,emy,emz;
   double Nparticle;
public:
  matrix* Beam_envelope;
   char b_res[Bres_max];
   //public:
   Beam() {particle=ELECTRON; mass=m_e; charge=-e;
	   energy=0.; emx=0.; emy=0.; emz=0.; Nparticle=0.; 
	   Beam_envelope=new matrix(6,6);
	 }
   Beam(const Beam& BEAM);
   Beam(const double*);
   Beam(const char*);
   ~Beam(void) {delete Beam_envelope;}

   double Energy(void) { return energy;}
   double Emittance_x(void) {return emx;}
   double Emittance_y(void) {return emy;}
   double Emittance_z(void) {return emz;}
   void SetEnergy(double E) { energy=E;}
   void SetEmittance(double ex,double ey,double ez)
     { emx=ex; emy=ey; emz=ez;}
   void SetNparticle(double n) { Nparticle=n;}
   double N_particle(void) {return Nparticle;}
   void ReadFile(const char*);
   void set_beam_envelope(matrix*);
   friend matrix BeamEnvelope(Beam& B) {return *(B.Beam_envelope);}

   void print(void);
   friend ostream& operator<<(ostream&, Beam&);
};

class EMIT
{
 public:
  double emx,emy,emz;
  double cod[6];
  double beig[6];
  matrix* Linear_map;
  matrix* PtoN;
  matrix* Beam_envelope;
  matrix* Damping_matrix;
  matrix* Diffusion_matrix;
  matrix* DifMeig;
  matrix* DifMeig_i_D;
  EMIT() {
    Linear_map=new matrix(6,6);
    PtoN=new matrix(6,6);
    Beam_envelope=new matrix(6,6);
    Damping_matrix=new matrix(6,6);
    Diffusion_matrix=new matrix(6,6);
    DifMeig=new matrix(6,6);
    DifMeig_i_D=new matrix(6,6);
    for(int i=0;i<6;i++) { cod[i]=0.;}
  }
  EMIT(const EMIT&);
  ~EMIT(void) {delete Beam_envelope; delete PtoN; delete Damping_matrix; 
	      delete Diffusion_matrix; delete Linear_map;
	      delete DifMeig; delete DifMeig_i_D; }
  void set_beam_envelope(matrix*);
  void set_EMIT(void);
  void set_EMIT_N(void);
  void DampRateReset(double,double,double,double,double,double);
  void SynchrotronRadiation(map_double&);
  void SynchrotronRadiation(pBeam&);

  friend matrix LinearMap(EMIT& B) {return *(B.Linear_map);}
  friend matrix Phys2Norm(EMIT& B) {return *(B.PtoN);}
  friend matrix BeamEnvelope(EMIT& B) {return *(B.Beam_envelope);}
  friend matrix DampingMatrix(EMIT& B) {return *(B.Damping_matrix);}
  friend matrix DiffusionMatrix(EMIT& B) {return *(B.Diffusion_matrix);}
  friend ostream& operator<<(ostream& s, EMIT&);
};


class Element
{
protected:
//   enum eltype { Drift=1, Bend, Quad, CBend, Thin };
//   eltype type;
   char name[16];
   double length;
public:

   // Constructor
   
   Element(void) {}
   Element(double* K) { length=K[0];}
   Element(double L) { length=L;}
   Element(char* s);

   Element(const Element&);
   Element& operator=(const Element&);
   double get_length(void) { return length;}
   virtual void print(void) 
     { cout << name << "(L= " << length ; }
   friend ostream& operator<<(ostream& s, Element&);
   friend ostream& operator<<(ostream& s, Element*);
   void SetName(const char* s) {strcpy(name,s);}

   virtual void Mapping(map_double& x)=0;
   virtual void Mapping(map_da& x)=0;
   virtual void Mapping(map_p_da& x) =0;
   virtual void Mapping(pBeam& x) =0;
   virtual void set_ent_edge(void) {;}
   virtual void set_exit_edge(void) {;}
   virtual void beam_set(Beam&) {;}
   virtual int FindIP(void) {return 0;}
   virtual double Luminosity(void) {return 0.;}
   virtual double dLuminosity(void) {return 0.;}
   virtual void LuminosityMonitorStart(void) {}
   virtual void IP_reset(double*) {}
   virtual void Xmon(double*) {;}
   virtual void Xmon(void) {;} 
};

void DriftMapping(double,map_double& x);
void DriftMapping(double,map_da& x);
void DriftMapping(double,map_p_da& x);
void DriftMapping(double,pBeam& x);

class Drift : public Element
{
public:
   Drift(double* K) 
        : Element(K) {;}
   Drift(double L) 
        : Element(L) {;}
   Drift(char* s) : Element(s) {;}
   void print(void) { Element::print(); cout << ")\n";}

   void Mapping(map_double& x) {DriftMapping(length,x);}
   void Mapping(map_da& x) {DriftMapping(length,x);}
   void Mapping(map_p_da& x) {DriftMapping(length,x);}
   void Mapping(pBeam& x) {DriftMapping(length,x);}
};

/*class CBend : public Element
{

public:
   CBend(double len,double teta,double K,double e1=0.,double e2=0.) 
          : Element(Element::CBend,len,teta,K,e1,e2)
   {
   //  cout<< "initialize CBend\n"; }
};*/



class Line_element
{
public:
   Element* M;
   double dx,dy,dtheta,cost,sint;
   double DR_space;
   double entrance,exit;
   Line_element* next;
   
   Line_element(void) {;}
   Line_element(char*);
   void print(void);
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   void beam_set(Beam& Beam) { M->beam_set(Beam);}
   friend ostream& operator<<(ostream& s, Line_element&);
   friend ostream& operator<<(ostream& s, Line_element*);
//   Line_element(void) { 
};

class Line_list
{
   Line_element* first;
   Line_element* last;
   double length;
public:
   Line_list(void) {first=last=new Line_element; length=0.;}
   Line_list(const char*);
   void clear() {last=0;}
   Line_list& insert(Line_element* a);
   Line_list& append(Line_element* a);
   void AppendDrift(double);
   void ReadFile(const char*);
   friend int ReadLine(const char*,int,char*,Line_list&,double);
   friend Line_element* ReadLineList(const char*,int,char*,Line_list&,double);
   friend int input(Line_list&);
   void print(void);
   void Mapping(map_double&);
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
   friend istream& operator>>(istream&, Line_element*);
   friend int input(int, char**, Line_list&);
   friend ostream& operator<<(ostream& s, Line_list&);
   void beam_set(Beam&);
   double Length() { return length;}
   friend void Print(Line_list& X) { X.print();}
   friend double Length(Line_list& X) { return X.length;}
   friend void OperatingParameterSet
     (const char* bfile,Beam& BEAM,EMIT& SAD,Line_list& ACC);
   int Luminosity(double*);
   double Luminosity(int n=1);
   void LuminosityMonitorStart(void);
   void IP_reset(double*);
};

class Accelerator: public Line_list
{
   double* rmu;
public:
   Accelerator(void): Line_list() {;}
   Accelerator(const char* mfile) : Line_list(mfile) {} 
   void Mapping(map_double& x) {Line_list::Mapping(x);}
   void Mapping(map_da&);
   void Mapping(map_p_da&);
   void Mapping(pBeam&);
};

void ReadBeam(char*,Beam& BEAM);
void ReadMathList(char* s,matrix& M);
int ReadMathList(char* s,double* M);

#endif
