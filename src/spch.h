#ifndef SPCH_H
#define SPCH_H
#include <element.h>
#include <matrix.h>

void InitPotentialSolver(int nx,int ny,int nz,
                         double dx,double dy,double dz);
void ClearPotentialSolver(void);
void SetCalcPotential(int it);
void UnSetCalcPotentialSolver(int it);
//void SpaceChargeParamSet(double rnp,double Mass,double Gamma);
//void SpaceChargeParamSet(double*,double**);
//void SpaceChargeParamSet(double**);
double **SpaceChargeParamSet(int nlatend, int isp);
//int SpaceChargeElementSize();
void SpaceChargeBPSet(double rnp,double mass,double Gamma)
//void SpaceChargeBeamSet(struct Beam*);
//void SpaceChargeMapping(double* x,double* px,double* y,double* py,
//                        double* z,double* pz,int np,
//                        double TravelLength);
//void SpaceChargeMapping(double*,double*,int);
void SpaceChargeMapping(double** SP,int isp,double* x,int np,
			double TravelLength)
void PotentialSolverPrint(void);
void epotkick(double* x,double* px,double* y,double* py,
                        double* z,double* pz,int np,
                        double TravelLength);
//void PrintPot(FILE*);
//void PrintPot(FILE*,double*);


#endif
