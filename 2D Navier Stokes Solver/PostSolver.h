/*
//PostSolver.h
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#include"Main.h"
#ifndef Sai_POSTSOLVER_H_
#define Sai_POSTSOLVER_H_

void PostSolver(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int Iter, int Comp);
void WritePltFile(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Iter, int Comp);
void WriteRestartData(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, int Iter, int Comp);
void CopyDatas(int NoOfElements, CELLDATAS *CD, CELLDATAS *CDO);
void WriteCPPlot(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Iter, int Comp);
void WriteClCdCm(CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Re, int Iter, int Comp);
void WriteWallStress(CELLDETAILS *CC, CELLDATAS *CD, int Iter, int Re, int Comp);
void WritePhase(CELLDATAS *CD, int Iter, int Comp);

#endif
