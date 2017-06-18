/*
//PreSolver.h
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/


#include"Main.h"
#ifndef Sai_PRESOLVER_H_
#define Sai_PRESOLVER_H_

#define CODE 5

#define INLET 1
#define OUTLET 2
#define TOP 3
#define BOTTOM 4
#define BODY1 5
#define DOM 0
//#define NNUNAVAILABLE 100

#define INNNUN 11
#define OUTNNUN 12
#define TOPNNUN 13
#define BOTNNUN 14
#define BODYNNUN 15

/*#define TOPRTCOR 302
#define TOPLTCOR 301
#define BOTRTCOR 401
#define BOTLTCOR 402
#define INTOPCOR 102
#define INBOTCOR 101
#define OUTTOPCOR 201 
#define OUTBOTCOR 202*/

#define UNAVAILABLEDATA -100

void PreSolver(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, NODE *Node, int Comp);
void ReadNeutralFiles(char *ElementFileName, char *NodeFileName, int NoOfElements, int NoOfNodes, CELLDETAILS *CD, NODE *Node);
void OrientCounterClockwise(int NoOfElements, CELLDETAILS *CC, NODE *Node);
void FindArea(int NoOfElements, CELLDETAILS *CC, NODE *Node);
void FindNeighbour(int Elem, int NoOfElements, CELLDETAILS *CC);
void ComputeNeighDetails(int NoOfElements, CELLDETAILS *CC, NODE *Node);
int CommonNodes(int CellP, int NodeA, int NodeB, CELLDETAILS *CC);
int FindCell(int CellP, int CellA, int NodeA, CELLDETAILS *CC);
int FindMissingNode(int P, int NodeA, int NodeB, CELLDETAILS *CC);
int FindIndex(int Cell, int Node, CELLDETAILS *CC);
void FindKs(int NoOfElements, CELLDETAILS *CC);
void WritePLTFILE(int NoOfElements, int NoOfNodes, CELLDETAILS *CC,  NODE *Node);
void ComputeNodeDistance(int NoOfNodes, CELLDETAILS *CC, NODE *Node);
#endif
