/*
//Solver.h
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#include"Main.h"
#ifndef Sai_SOLVER_H_
#define Sai_SOLVER_H_

#define MAXPPECOUNT 1500
#define MAXDIV 0.0001

void Solver(int NoOfElements, int NoOfNodes, int NoOfBodies, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, BODY_ELEM *BE, int **BodyElems, int Re, double DelT,int Count);
void SolvePPE(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node);
//void ApplyBC(int NoOfElements, CELLDETAILS *CC, CELLDATAS *C);
void WallBC(CELLDATAS *C, int Elem);
void DirichletBC(CELLDATAS *C, int Elem);
void InletBC(CELLDATAS *C, int Elem);
void OutletBC(CELLDATAS *C, int Elem);
void UpdateDerivatives(int NoOfElements, CELLDETAILS *CC, CELLDATAS* C, NODE *Node);
double CalculateP(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT,int Count);
void CalculateSource(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node);
STENCIL UpdateStencil(CELLDATAS *C, int Elem, CELLDETAILS *CC, NODE *N);
void UpdatePressure(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO);
void SolveCFRNS(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT);
void SolveNS(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int count);
void DirichletExtra(CELLDETAILS *CC, CELLDATAS *C, int Elem);

void ApplyBC(int NoOfElements, int NoOfNodes, int NoOfBodies, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, BODY_ELEM *BE, int **BodyElems, int Count);
void ComputeNodalDatas(int NoOfNodes, CELLDATAS *CDO, NODE *Node);
void UpdateVelocity(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD);

void SolverHO(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT,int Count);
void SolvePPEHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node);
double CalculatePHO(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT,int Count);
void CalculateSourceHO(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node);
void SolveCFRNSHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT);
void SolveNSHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int count);

void UpdateTemperature(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CDO);
void SolveTempHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS * CDO, NODE *Node, int Re, double DelT);
void SolveTemp(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS * CDO, NODE *Node, int Re, double DelT);

void MassBalance(int NoOfElements,CELLDETAILS *CC, CELLDATAS *CD);
#endif
