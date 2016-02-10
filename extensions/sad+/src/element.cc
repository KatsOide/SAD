#include <iostream>
using std::ostream;
using std::ios;
using std::cerr;
#include <fstream>
using std::ifstream;
#include <cstdlib>
using std::atof;
using std::exit;
#include <cstring>
using std::strchr;
using std::strcmp;
using std::strcspn;
using std::strlen;
using std::strncpy;
using std::strstr;
#include <cctype>
using std::iscntrl;
#include <cmath>
using std::fabs;
using std::sin;
using std::cos;

#include <element.h>
#include <lin_map.h>
#include <bend.h>
#include <quad.h>
#include <sext.h>
#include <thin.h>
#include <ph_rot.h>
#include <bb.h>
#include <cavity.h>
#include <wiggler.h>
#include <solenoid.h>
#include <impedance.h>
#include <linelem.h>

#define Npmax 100

void SetRandum(int);

enum { DRIFT=1, BEND, QUAD, CBEND, SEXT, OCTU, THIN, CAV, INTP, PHROT,
       WIGGLER, SOLENOID, LBEND, IMPEDANCE, CHROM, EXIT };

Element::Element(const Element& x)
{
   length=x.length;
}
Element& Element::operator=(const Element& x)
{
   length=x.length;
   return *this;
}


//-------------------------------------------------------------------------
//  mapping of Line_list
//-------------------------------------------------------------------------

void Line_list::Mapping(map_da& x)
{
   Line_element* ptr=first;
   while(ptr!=last) {
      ptr->Mapping(x);
      ptr=ptr->next;
   }
}

void Line_list::Mapping(map_double& x)
{
  Line_element* ptr=first;
  int i=0;

  while(ptr!=last) {
    ptr->Mapping(x);
    if(!is_survive(x)) {
      cout << "Particle loss\n";
      break;
    }
    i++;
//
//  For Debug
//
//    if(i==529) {
//      cout<< "529\n";
//    }
//  cout.width(4); cout << "**** "<< i << "  ";
//  ptr->print();
    cout.setf(ios::scientific,ios::floatfield);
    cout.precision(15);
//  cout << x;
    cout.flush();
    cout.unsetf(ios::scientific);
    cout.unsetf(ios::floatfield);
    cout.precision(6);
    
    ptr=ptr->next;
  }
}

void Line_list::Mapping(pBeam& x)
{
  static int flag=0; 
  Line_element* ptr=first;
  int i=0;
  while(ptr!=last) {
    ptr->Mapping(x);
    flag++;
    if(flag%5==0) is_survive(x);
    i++;
//      {ptr->print();cout << x2.np << x2;}
//    if(i==528) {
//      cout<< "528\n";
//    }
//      ptr->print();
//      cout << x2;
    ptr=ptr->next;
  }
}

int Line_list::Luminosity(double* Lum)
{
  Line_element* ptr=first;
  int i=0;
  while(ptr!=last) {
    if(ptr->M->FindIP()) {
      Lum[i]=ptr->M->Luminosity();
      i++;
    }
    ptr=ptr->next;
  }
  return i;
}

double Line_list::Luminosity(int n)
{
  Line_element* ptr=first;
  int i=0;
  while(ptr!=last) {
    if(ptr->M->FindIP()) {
      i++;
      if(n==i) return ptr->M->Luminosity();
    }
    ptr=ptr->next;
  }
  cout << "There is no " << n << "-th IP\n";
  return 0.0;
}
  
void Line_list::LuminosityMonitorStart(void)
{
  Line_element* ptr=first;
  int i=0;
  while(ptr!=last) {
    if(ptr->M->FindIP()) {
      i++;
      ptr->M->LuminosityMonitorStart();
      cout << "\n ****** Luminosity monitor start at " << i 
	<< "-th IP ***************\n\n";
    }
    ptr=ptr->next;
  }
}

void Line_list::IP_reset(double* K)
{
  Line_element* ptr=first;

  while(ptr!=last) {
    if(ptr->M->FindIP()) {
      ptr->M->IP_reset(K);
      cout << " ############ IP was reset ###############################\n";
      ptr->M->print();
    }
    ptr=ptr->next;
  }
}
  
void Accelerator::Mapping(map_da& x)
{
   Line_list::Mapping(x);
}


void Accelerator::Mapping(pBeam& x)
{
   Line_list::Mapping(x);
}

//-------------------------------------------------------------------------
//  mapping of Line_element
//-------------------------------------------------------------------------

void Line_element::Mapping(map_da& x)
{
  if(DR_space!=0.) DriftMapping(DR_space,x);

  if(dx!=0.) x[0]-=dx;
  if(dy!=0.) x[1]-=dy;
  if(dtheta!=0.) { 
    da zi=x[0];
    x[0]=zi*cost-x[1]*sint;
    x[1]=zi*sint+x[1]*cost;
    zi=x[3];
    x[3]=zi*cost-x[4]*sint;
    x[4]=zi*sint+x[4]*cost;
  }
  M->Mapping(x);
  if(dtheta!=0.) {
    da zi=x[0];
    x[0]=zi*cost+x[1]*sint;
    x[1]=x[1]*cost-zi*sint;
    zi=x[3];
    x[3]=zi*cost+x[4]*sint;
    x[4]=x[4]*cost-zi*sint;
  }
  if(dx!=0.) x[0]+=dx;
  if(dy!=0.) x[1]+=dy;
}

void Line_element::Mapping(map_double& x)
{
  if(DR_space!=0.) DriftMapping(DR_space,x);

  if(dx!=0.) x[0]-=dx;
  if(dy!=0.) x[1]-=dy;
  if(dtheta!=0.) { 
    double zi=x[0];	// ########### difference between this and da.
    x[0]=zi*cost-x[1]*sint;
    x[1]=zi*sint+x[1]*cost;
    zi=x[3];
    x[3]=zi*cost-x[4]*sint;
    x[4]=zi*sint+x[4]*cost;
  }
  M->Mapping(x);
  if(dtheta!=0.) {
    double zi=x[0];	// ###########
    x[0]=zi*cost+x[1]*sint;
    x[1]=x[1]*cost-zi*sint;
    zi=x[3];
    x[3]=zi*cost+x[4]*sint;
    x[4]=x[4]*cost-zi*sint;
  }
  if(dx!=0.) x[0]+=dx;
  if(dy!=0.) x[1]+=dy;
}

void Line_element::Mapping(pBeam& x)
{
  int i;

  if(DR_space!=0.) DriftMapping(DR_space,x);

  if(dx!=0.) {for(i=0;i<x.np;i++) x.x[i]-=dx;}
  if(dy!=0.) {for(i=0;i<x.np;i++) x.y[i]-=dy;}
  if(dtheta!=0.) {
    for(i=0;i<x.np;i++) {
      double zi=x.x[i];	// ########### difference between this and da.
      x.x[i]=zi*cost-x.y[i]*sint;
      x.y[i]=zi*sint+x.y[i]*cost;
      zi=x.px[i];
      x.px[i]=zi*cost-x.py[i]*sint;
      x.py[i]=zi*sint+x.py[i]*cost;
    }
  }
  M->Mapping(x);
  if(dtheta!=0.) {
    for(i=0;i<x.np;i++) {
      double zi=x.x[i];	// ###########
      x.x[i]=zi*cost+x.y[i]*sint;
      x.y[i]=x.y[i]*cost-zi*sint;
      zi=x.px[i];
      x.px[i]=zi*cost+x.py[i]*sint;
      x.py[i]=x.py[i]*cost-zi*sint;
    }
  }
  if(dx!=0.) { for(i=0;i<x.np;i++) x.x[i]+=dx; }
  if(dy!=0.) { for(i=0;i<x.np;i++) x.y[i]+=dy; }
}
//-------------------------------------------------------------------------

void Line_list::print(void)
{
   cout << "\n####### Machine resource print out ############################\n";
   Line_element* ptr=first;
   while(ptr!=last) {
      ptr->print();
      ptr=ptr->next;
   }
   cout<< " Total length = " << length << "m\n\n";
}

/*
void Print(Line_list& ACC)
{
  ACC.print();
}*/

ostream& operator<<(ostream& s,Element& M)
{
   s<< M.name << "   L = " << M.length << '\n';
   return s;
}
ostream& operator<<(ostream& s,Element* M)
{
   s<< M->name << " (L = " << M->length << ")";
   return s;
}

ostream& operator<<(ostream& s,Line_list& ACC)
{
   s << "\n####### Machine resource print out ############################\n";
   Line_element* ptr=ACC.first;
   while(ptr!=ACC.last) {
      s << ptr;
      ptr=ptr->next;
   }
   s << " Total length = " << ACC.length << "m\n\n";
   return s;
}

void Line_element::print(void)
{
   M->print();
   cout.width(8);
   cout << '(' << DR_space << " )";
   cout.width(8); cout <<  entrance << "  "; 
   cout.width(8); cout <<  exit << "\n";
   if(dx!=0. || dy!=0. || dtheta!=0.) {
     cout.width(8);
     cout << "       : " << dx << "  ";
     cout.width(8); cout << dy << "  ";
     cout.width(8); cout << dtheta << "\n";
   }
}

ostream& operator<<(ostream& s, Line_element& LE)
{
   s<< LE.M;
   s.width(8);
   s << '(' << LE.DR_space << " )";
   s.width(8); s << LE.entrance << "  ";
   s.width(8); s << LE.exit << "\n";
   if(LE.dx!=0. || LE.dy!=0. || LE.dtheta!=0.) {
     s.width(8);
     s << "      : " << LE.dx << "  ";
     s.width(8); s << LE.dy << "  ";
     s.width(8); s << LE.dtheta << "\n";
   }
   return s;
}

ostream& operator<<(ostream& s, Line_element* LE)
{
   s<< LE->M;
   if(LE->dx!=0. || LE->dy!=0. || LE->dtheta!=0.)
   s << "   dx = " << LE->dx << "  dy = " << LE->dy
      << "  dtheat = " << LE->dtheta << ";\n";
   s << " : " << LE->entrance << " \n";
   return s;
}


void Line_list::beam_set(Beam& BEAM)
{
  //   cout << " Beam parameter set in IP element \n";
  // cout << " Beam parameter is as follow,\n";
  // BEAM.print();
  // cout <<"\n\n";
   Line_element* ptr=first;
   while(ptr!=last) {
      ptr->beam_set(BEAM);
      ptr=ptr->next;
   }
}


void get_name(const char* s, char* p)
{
  size_t i;
   for(i=0;i<strlen(s);i++) {
     if(s[i]!=' ') break;
   }
   int n=strcspn(s+i,"=( ");
   if(n==0) {p[0]='\0'; return;}
   if(n>=15) {cout << "Element name is too long >15"; exit(1);}
   strncpy(p,s+i,n);
   p[n]='\0';
}

double get_parm(char* s,const char* parm)
{
   char *p,*p2;
   int ibc,ifd;
   int lpar=strlen(parm);
   if((p=strchr(s,'('))==NULL) p=s;
   while((p=strstr(p,parm))!=NULL) {
      if(p==s) ibc=1; 
      else if(*(p-1)==' ' || *(p-1)==',' || *(p-1)=='(') ibc=1;
      else ibc=0;
      if(*(p+lpar)==' ' || *(p+lpar)=='=') ifd=1; else ifd=0;
      p2=strchr(p,'=')+1;
      if(ibc && ifd) return atof(p2);
      p=p2;
   }
   return 0.;
}

Element::Element(char* s)
{
   get_name(s,name);
   length=get_parm(s,"L");
}


Line_element::Line_element(char* s)
{
   dx=get_parm(s,"dx");
   dy=get_parm(s,"dy");
   dtheta=get_parm(s,"dtheta");
   DR_space=get_parm(s,"drift");
   cost=cos(dtheta);
   sint=sin(dtheta);
}


int el_type_number(char* s)
{
   if(!strcmp(s,"Drift")) return DRIFT;
   if(!strcmp(s,"Bend")) return BEND;
   if(!strcmp(s,"Quad")) return QUAD;
   if(!strcmp(s,"Sext")) return SEXT;
   if(!strcmp(s,"Octu")) return OCTU;
   if(!strcmp(s,"Thin")) return THIN;
   if(!strcmp(s,"Cavity")) return CAV;
   if(!strcmp(s,"IP")) return INTP;
   if(!strcmp(s,"Ph_rot")) return PHROT;
   if(!strcmp(s,"Wiggler")) return WIGGLER;
   if(!strcmp(s,"Solenoid")) return SOLENOID;
   if(!strcmp(s,"LBend")) return LBEND;
   if(!strcmp(s,"Impedance")) return IMPEDANCE;
   if(!strcmp(s,"Chromaticity")) return CHROM;
   if(!strcmp(s,"Exit")) return EXIT;
   cout << "The element " << s << " is not identified\n";
   return 0;
}

void OperatingParameterSet
(const char* bfile,Beam& BEAM,EMIT& SAD,Line_list& ACC)
{
  char b_res[1024],res_type[16],name[16];
  int el_type;
  double driftspace=0.;
  Line_element* Last_element=0;

  N_cv=6;
  N_var=6;
  N_ord=4;

  ifstream op_resource(bfile,ios::in);
  if(!op_resource) {
    cerr << "cannot open input file : " << bfile << '\n';
    cout << " Default values are used as operating parameters.\n";
    return;
  }

  char *bp,*p;
  int n,ERR=0,i,COMMENT=0;
  char dlim;
  while(!op_resource.eof()) {
    COMMENT=0;
    op_resource.get(b_res,Bres_max,';');
    if(!op_resource) ERR=1;
    op_resource >> dlim;
    if(dlim!=';' || !op_resource) ERR=1;
    n=strlen(b_res);
    if(n==0) ERR=2;
    if(n>=Bres_max) {
      cout << " operator resource content is too long\n";
      exit(1);
    }
    for(i=0;i<n;i++) {      /* Set blank for control sequence */
      if(iscntrl(b_res[i])) b_res[i]=' ';
    }
    for(i=0,bp=b_res;i<n;i++) {      /* Search head position of resource */
      if(b_res[i]!=' ') {
	bp=b_res+i;
	break;
      }
    }
/*   Resource check  */
    if(*bp=='/' && *(bp+1)=='/') COMMENT=1;
    if(*bp=='#') COMMENT=1;
    if(strlen(bp)==0 || strchr(bp,'(')==NULL || strchr(bp,')')==NULL) {
      cout << "Element could not be interpleted " << bp << '\n';
      ERR=3;
    }

    if(!ERR && !COMMENT) {

      for(i=0;i<n;i++) {      /* Search head position of resource */
	if(bp[i]==' ') {
	  break;
	}
      }

      strncpy(res_type,bp,i);
      res_type[i]='\0';
      bp+=i;
      get_name(bp,name);

//      str
      if(!strcmp(res_type,"Operator")) N_ord=(int)get_parm(bp,"Nord");
      else if(!strcmp(res_type,"Beam")) ReadBeam(bp,BEAM);

      else if(!strcmp(res_type,"Linmap"))
	ReadMathList(bp,*SAD.Linear_map);
      else if(!strcmp(res_type,"LinmapD")) 
	ReadMathList(bp,*SAD.Damping_matrix);
      else if(!strcmp(res_type,"Bmatrix")) 
	ReadMathList(bp,*SAD.Diffusion_matrix);
      else if(!strcmp(res_type,"Benv"))
	ReadMathList(bp,*SAD.Beam_envelope);
      else if(!strcmp(res_type,"PtoN"))
	ReadMathList(bp,*SAD.PtoN);
      else if(!strcmp(res_type,"COD"))
	ReadMathList(bp,SAD.cod);
      else if((el_type=el_type_number(res_type))==DRIFT) {
	p=strchr(bp,'{')+1;
	driftspace+=atof(p);
      }
      else if(el_type==EXIT) {
	if(driftspace!=0.) {
	  if(Last_element!=0) {
	    Last_element->M->set_exit_edge();}
	  ACC.AppendDrift(driftspace);
	}
	break;
      }
      else if(el_type) {
	Last_element=ReadLineList(name,el_type,bp,ACC,driftspace);
	driftspace=0;
      }
      else {
	cout << "******** Undefined element *********\n";
      }
    }
  }
  if(driftspace!=0.) {
    if(Last_element!=0) {
      Last_element->M->set_exit_edge();}
    ACC.AppendDrift(driftspace);
  }
  op_resource.close();

  ACC.beam_set(BEAM);
//  SAD.set_EMIT();
  SAD.set_EMIT_N();
  SetRandum(17);
//  ACC.egde_set();
}

int ReadLine(const char* name,int eltype,char* resource,Line_list& ACC,
	     double driftspace)
{
// Read Machine resource

  char *p3,*setting_resource,*elem_resource;
  static Line_element* ptrstore=0;
  static int eltstore=0;

  Line_element* ptr=ACC.last;

  p3=strchr(resource,')');
  setting_resource=p3+1;
  elem_resource=resource;

  ptr->dx=get_parm(setting_resource,"dx");
  ptr->dy=get_parm(setting_resource,"dy");
  ptr->dtheta=get_parm(setting_resource,"dtheta");
  ptr->DR_space=get_parm(setting_resource,"drift");
  ptr->cost=cos(ptr->dtheta);
  ptr->sint=sin(ptr->dtheta);
  if(fabs(ptr->dtheta-pi/2.)<1.e-13) {
    ptr->cost=0.; ptr->sint=1.;
  }
  if(fabs(ptr->dtheta+pi/2.)<1.e-13) {
    ptr->cost=0.; ptr->sint=-1.;
  }

  switch(eltype) {
      case DRIFT:	ptr->M=new Drift(elem_resource);	break;
      case BEND:	ptr->M=new Bend(elem_resource);		break;
      case QUAD:	ptr->M=new Quad(elem_resource);		break;
      case SEXT:	ptr->M=new Sext(elem_resource);		break;
      case OCTU:	ptr->M=new Octu(elem_resource);		break;
      case THIN:	ptr->M=new Thin(elem_resource);		break;
      case CAV:		ptr->M=new Cavity(elem_resource);	break;
      case INTP:	ptr->M=new IP(elem_resource);		break;
      case PHROT:	ptr->M=new Ph_rot(elem_resource);	break;
      case WIGGLER:	ptr->M=new Wiggler(elem_resource);	break;
      case SOLENOID:	ptr->M=new Solenoid(elem_resource);	break;
      case LBEND:	ptr->M=new LBend(elem_resource);	break;
      case IMPEDANCE:	ptr->M=new Impedance(elem_resource);	break;
      case CHROM:	ptr->M=new Chromaticity(elem_resource);	break;
      case EXIT:        ptr->M=0;
  }
  if(ptr->M==0) return 0;
  ptr->M->SetName(name);

  ACC.last=new Line_element;
  ptr->next=ACC.last;
  ptr->entrance=ACC.length;
  ptr->exit=ACC.length+ptr->M->get_length()+ptr->DR_space;
  ACC.length=ptr->exit;

//  Qedge  set  
  if(eltype==QUAD && driftspace!=0.) ptr->M->set_ent_edge();
  if(eltstore==QUAD && driftspace!=0.) ptrstore->M->set_exit_edge();
  
  ptrstore=ptr;
  eltstore=eltype;
  return 1;
}

Line_element* ReadLineList(const char* name,int eltype,char* resource,
			   Line_list& ACC,double driftspace)
{
// Read Machine resource

  double K[Npmax];
  static Line_element* ptrstore=0;
  static int eltstore=0;

  Line_element* ptr=ACC.last;

  char* p  = strchr(resource,'{') + 1;
#if 0
  char* pf = strchr(resource,'}');
#endif
  char* p1;
  int i = 0;
  while((p1 = strchr(p,',')) != NULL) {
    K[i]=atof(p);
    p=p1+1;
    i++;
    if(i>Npmax) {
      cout << " Number of parameter exceed a limit " << Npmax;
      exit(1);
    }
  }
  K[i]=atof(p);

  ptr->dx=0.;
  ptr->dy=0.;
  ptr->dtheta=0.;
  ptr->DR_space=0.;

  switch(eltype) {
      case BEND: ptr->dx=K[7]; ptr->dy=K[8]; ptr->dtheta=K[4];
		ptr->M=new Bend(K);
		break;
      case QUAD:	ptr->M=new Quad(K);
		ptr->dx=K[3]; ptr->dy=K[4]; ptr->dtheta=K[2];
		break;
      case SEXT:	ptr->M=new Sext(K);
		ptr->dx=K[3]; ptr->dy=K[4]; ptr->dtheta=K[2];
		break;
      case OCTU:	ptr->M=new Octu(K);
		ptr->dx=K[3]; ptr->dy=K[4]; ptr->dtheta=K[2];
		break;
      case THIN:	ptr->M=new Thin(K);
		break;
      case CAV:		ptr->M=new Cavity(K);
		ptr->dx=K[9]; ptr->dy=K[10]; ptr->dtheta=K[11];
		break;
      case INTP:	ptr->M=new IP(K);
		ptr->dx=0.;
		ptr->dy=0.;
		ptr->dtheta=0.;
		break;
      case PHROT:	ptr->M=new Ph_rot(K);
		break;
      case WIGGLER:	ptr->M=new Wiggler(K);	break;
      case SOLENOID:	ptr->M=new Solenoid(K);	break;
      case LBEND:	ptr->M=new LBend(K);	break;
      case IMPEDANCE:	ptr->M=new Impedance(K);	break;
      case CHROM:	ptr->M=new Chromaticity(K);	break;
//      case EXIT:        if(driftspace!=0.) {
//	                  ptr->M=new Drift(driftspace);
//			}
  }

  ptr->M->SetName(name);
  ptr->cost=cos(ptr->dtheta);
  ptr->sint=sin(ptr->dtheta);
  if(fabs(ptr->dtheta-pi/2.)<1.e-13) {
    ptr->cost=0.; ptr->sint=1.;
  }
  if(fabs(ptr->dtheta+pi/2.)<1.e-13) {
    ptr->cost=0.; ptr->sint=-1.;
  }
  if(eltype!=DRIFT) {
    ptr->DR_space=driftspace;
  }

  ACC.last=new Line_element;
  ptr->next=ACC.last;
  ptr->entrance=ACC.length+ptr->DR_space;
  ptr->exit=ACC.length+ptr->M->get_length()+ptr->DR_space;
  ACC.length=ptr->exit;

//  Qedge  set  
  if(eltype==QUAD && driftspace!=0.) ptr->M->set_ent_edge();
  if(eltstore==QUAD && driftspace!=0.) ptrstore->M->set_exit_edge();

  ptrstore=ptr;
  eltstore=eltype;

  if(eltype==EXIT) return 0;
  return ptr;
}

void Line_list::AppendDrift(double L)
{

  Line_element* ptr=last;

  ptr->dx=0.;
  ptr->dy=0.;
  ptr->dtheta=0.;
  ptr->DR_space=0.;
  ptr->M=new Drift(L);

  ptr->M->SetName("Exit");
  ptr->cost=cos(ptr->dtheta);
  ptr->sint=sin(ptr->dtheta);

  last=new Line_element;
  ptr->next=last;
  ptr->entrance=length;
  ptr->exit=length+ptr->M->get_length()+ptr->DR_space;
  length=ptr->exit;

}

