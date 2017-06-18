/*
//Main.h
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#ifndef Sai_MAIN_H_
#define Sai_MAIN_H_
#include <stdio.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 Every Node is represented P(x,y)
 The Type represents the location of the point like,
 DOM --> Fully immersed in the flow
 TOP --> point lies on the Top Boundary
 BOTTOM --> point lies on the Bottom Boundary
 INLET --> point lies on the Inlet Boundary
 OUTLET --> point lies on the Outlet Boundary
 BODY1 --> Point lies on the Body. The Number represents the Body Number.
 */

typedef struct
{
	double x,y;
	int Type;
	int ShareElements[10];
	int ShareCount;
	/*double CC2NodeDist[10];
	/double u,v;*/
	double CC2NodeDist[10];
	double Sum;
	double u,v;
	double T;
}NODE;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 *Every Triangular Cell is represented by the connection of 3 nodes named as "a", "b", "c"
 *The connectivity information is stored in Connect[3]. Arranged in counterclockwise orientation

 Xab, Xbc, Xca represents the X component distance of the line joining ab, bc and ca respectively.		<-------------------
 Yab, Ybc, Yca represents the Y component distance of the line joining ab, bc and ca respectively.		<---------------	|
																														|	|
 Area Represents the area of the cell.																					|	|
																														|	|
 Neigh[3] holds the three neighbouring cells. if any neighour is available it will be assigned to NEIGHUNAVAILABLE		|	|
 CellD, CellE, CellG, CellF, CellI, CellJ, CellL represents the cell informations in the stencil						|	|
																														|	|
 Ybd, Ydc, Ycj, Yja, Yag, Ygb; Represents similar to --------------------------------------------------------------------	|
 Xbd, Xdc, Xcj, Xja, Xag, Xgb; Represents similar to ------------------------------------------------------------------------
 Apa, Apb, Apc represents Area of parent cell + Area of CELLA,Area of parent cell + Area of CELLB, Area of parent cell + Area of CELLCrespectively

 Index4, Index5, Index11, Index12, Index18, Index19 represents the location of the parematers in the stencil

 Ka, Kb, Kc, Kd, Kf, Kg, Ki, Kj, Kl, Kp represents the Constants used in the final form of the Pressure poisson Equation.

 The Type represents the location of the Cell like,
 DOM --> Fully immersed in the flow
 TOP --> one cell face (FACE1) lies on the Top Boundary
 BOTTOM --> one cell face (FACE1) lies on the Bottom Boundary
 INLET --> one cell face (FACE1) lies on the Inlet Boundary
 OUTLET --> one cell face (FACE1) lies on the Outlet Boundary
 BODY1 --> one cell face (FACE1) lies on the Body. The Number represents the Body Number.


 Note:
	 If Two of the cell faces are lying on the boundary the code crashes.
 */


typedef struct
{
	int Connect[3];
	int Neighbour[3];
	int CellD, CellL, CellJ, CellI, CellG, CellF;
	int Index4, Index5, Index11, Index12, Index18, Index19;
	int Type;
		
/*	double Xab, Xbc, Xca;
	double Yab, Ybc, Yca;
	double Area;
	double Xp, Yp;
	double Xag, Xbd, Xcj, Xdc, Xgb, Xja;
	double Yag, Ybd, Ycj, Ydc, Ygb, Yja;
	double Apa, Apb, Apc;
	double XAP, XBP, XCP;
	double YAP, YBP, YCP;
	double K1, K2, K3, K4, K5, K6, K7, K8, K9, K10;
*/
	double Xab, Xbc, Xca;
	double Yab, Ybc, Yca;
	double Area;
	double Xp, Yp;
	double Xag, Xbd, Xcj, Xdc, Xgb, Xja;
	double Yag, Ybd, Ycj, Ydc, Ygb, Yja;
	double Apa, Apb, Apc;
	double XAP, XBP, XCP;
	double YAP, YBP, YCP;
	double Ka, Kb, Kc, Kd, Kf, Kg, Ki, Kj, Kl, Kp;

}CELLDETAILS;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * in General CC represents the Cell Centered values and FC represents the Face Centered
 *
 * CCVel represents The Cell Centered Values of the Velocities u,v
 * CCP represents the Cell Centered Pressure
 * du_dx = ∂u/∂x, dv_dx = ∂v/∂x, du_dy = ∂u/∂y, dv_dy = ∂v/∂y
 */

typedef struct
{
	double CCVel[2];
	double FCVel[3][2];
	double CCP;
	double FCP[3];
	double CCT;
	double FCT[3];
	double FCdu_dx[3], FCdu_dy[3], FCdv_dx[3], FCdv_dy[3];
	double FCdT_dx[3], FCdT_dy[3];
	double Source;
	double XCFLUX1, YCFLUX1, XDFLUX1, YDFLUX1;
	double XCFLUX2, YCFLUX2, XDFLUX2, YDFLUX2;
	double XCFLUX3, YCFLUX3, XDFLUX3, YDFLUX3;
	double XPFLUX1, XPFLUX2, XPFLUX3;
	double YPFLUX1, YPFLUX2, YPFLUX3;
	double CCdudx, CCdudy, CCdvdx, CCdvdy;
	double Vorticity, Divergence;
}CELLDATAS;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Represents the corresponding parameters in the stencil
 */

typedef struct
{
	double u1, u2, u3, u4, u5, u11, u12, u18, u19;
	double v1, v2, v3, v4, v5, v11, v12, v18, v19;
	double p1, p2, p3, p4, p5, p11, p12, p18, p19;
	double dudx1, dudx2, dudx3, dudx4, dudx5, dudx11, dudx12, dudx18, dudx19;
	double dvdx1, dvdx2, dvdx3, dvdx4, dvdx5, dvdx11, dvdx12, dvdx18, dvdx19;
	double dudy1, dudy2, dudy3, dudy4, dudy5, dudy11, dudy12, dudy18, dudy19;
	double dvdy1, dvdy2, dvdy3, dvdy4, dvdy5, dvdy11, dvdy12, dvdy18, dvdy19;
}STENCIL;

typedef struct
{
	int LeadingEdge, TrailingEdge;
	double Amplitude, Frequency;
	int MaxElems;
}BODY_ELEM;


void GetNoOfElemAndNodes(int *, int *);
void InitialConditions(int NoOfElements, int NoOfNodes, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node);
int RestartCheck(int NoOfElements, int NoOfNodes, int Re);
#define MAXITER 1250000
#define DELT 0.0005
#define RE 100
#define Pr 0.71
//#define ALPHA 1 // NonDimensional Tangential Velocity
//#define FREQ 1 //Frequency of Oscilation
#define Uinf 1.0
#define Vinf 0.0
#define Pinf 0.0
#define Tinf 0.0

#define Tinit 0.0
#define Uinit 1.0
#define Vinit 0.0
#define Pinit 0.0
#define TSource 1.0
#define U 0
#define V 1
#define One 0
#define Two 1
#define Three 2 
#define PI 3.141592654
#define MAXBODIES 2
#define MAXELEMSONBODY 200
#endif

