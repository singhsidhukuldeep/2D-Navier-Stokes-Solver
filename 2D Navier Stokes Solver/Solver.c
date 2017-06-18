/*
//Solver.c
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#include"Solver.h"
#include"PreSolver.h"
#include"Main.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<conio.h>

			//SO = UpdateStencil(C, Elem, CC);Gaussian Approximation for derivatives

			/*C[Elem].FCdu_dx[0] = (SO.u2 * CC[Elem].Ybc + SO.u3 * CC[Elem].Yca + SO.u12 * CC[Elem].Yag + SO.u11 * CC[Elem].Ygb) / CC[Elem].Apc ;
			C[Elem].FCdu_dx[1] = (SO.u1 * CC[Elem].Yab + SO.u5 * CC[Elem].Ybd + SO.u4 * CC[Elem].Ydc + SO.u3 * CC[Elem].Yca) / CC[Elem].Apa ;
			C[Elem].FCdu_dx[2] = (SO.u1 * CC[Elem].Yab + SO.u2 * CC[Elem].Ybc + SO.u19 * CC[Elem].Ycj + SO.u18 * CC[Elem].Yja) / CC[Elem].Apb ;

			C[Elem].FCdv_dx[0] = (SO.v2 * CC[Elem].Ybc + SO.v3 * CC[Elem].Yca + SO.v12 * CC[Elem].Yag + SO.v11 * CC[Elem].Ygb) / CC[Elem].Apc ;
			C[Elem].FCdv_dx[1] = (SO.v1 * CC[Elem].Yab + SO.v5 * CC[Elem].Ybd + SO.v4 * CC[Elem].Ydc + SO.v3 * CC[Elem].Yca) / CC[Elem].Apa ;
			C[Elem].FCdv_dx[2] = (SO.v1 * CC[Elem].Yab + SO.v2 * CC[Elem].Ybc + SO.v19 * CC[Elem].Ycj + SO.v18 * CC[Elem].Yja) / CC[Elem].Apb ;

			C[Elem].FCdu_dy[0] = -1 * (SO.u2 * CC[Elem].Xbc + SO.u3 * CC[Elem].Xca + SO.u12 * CC[Elem].Xag + SO.u11 * CC[Elem].Xgb) / CC[Elem].Apc ;
			C[Elem].FCdu_dy[1] = -1 * (SO.u1 * CC[Elem].Xab + SO.u5 * CC[Elem].Xbd + SO.u4 * CC[Elem].Xdc + SO.u3 * CC[Elem].Xca) / CC[Elem].Apa ;
			C[Elem].FCdu_dy[2] = -1 * (SO.u1 * CC[Elem].Xab + SO.u2 * CC[Elem].Xbc + SO.u19 * CC[Elem].Xcj + SO.u18 * CC[Elem].Xja) / CC[Elem].Apb ;

			C[Elem].FCdv_dy[0] = -1 * (SO.v2 * CC[Elem].Xbc + SO.v3 * CC[Elem].Xca + SO.v12 * CC[Elem].Xag + SO.v11 * CC[Elem].Xgb) / CC[Elem].Apc ;
			C[Elem].FCdv_dy[1] = -1 * (SO.v1 * CC[Elem].Xab + SO.v5 * CC[Elem].Xbd + SO.v4 * CC[Elem].Xdc + SO.v3 * CC[Elem].Xca) / CC[Elem].Apa ;
			C[Elem].FCdv_dy[2] = -1 * (SO.v1 * CC[Elem].Xab + SO.v2 * CC[Elem].Xbc + SO.v19 * CC[Elem].Xcj + SO.v18 * CC[Elem].Xja) / CC[Elem].Apb ;
				*/

//void Solver(int NoOfElements, int NoOfNodes, int NoOfBodies, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, BODY_ELEM *BE, int **BodyElems, int Re, double DelT,int Count)
void Solver(int NoOfElements, int NoOfNodes, int NoOfBodies, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, BODY_ELEM *BE, int **BodyElems, int Re, double DelT,int Count)
{
	int Elem;
	Re = RE;
	DelT = DELT;
	ComputeNodalDatas(NoOfNodes, CDO, Node);
	ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
	UpdateVelocity(NoOfElements, CC, CDO);
	ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
	ComputeNodalDatas(NoOfNodes, CDO, Node);
	ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
	UpdateDerivatives(NoOfElements, CC, CDO, Node);
	
	switch(Count)
	{
	case 1:
		for(Elem=0;Elem<NoOfElements;Elem++)
		{
			if(CC[Elem].Type == OUTLET)
				CD[Elem].FCVel[0][0] = Uinf;
				CD[Elem].FCVel[0][1] = Vinf;

				CDO[Elem].FCVel[0][0] = Uinf;
				CDO[Elem].FCVel[0][1] = Vinf;
		}
		SolvePPE(NoOfElements, CC, CD, CDO, Re, DelT, Node);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		UpdatePressure(NoOfElements, CC, CD, CDO);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		SolveCFRNS(NoOfElements, CC, CD, CDO, Node, Re, DelT);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		SolveNS(NoOfElements, CC, CD, CDO, Node, Re, DelT, Count);
		SolveTemp(NoOfElements, CC, CD, CDO, Node, Re, DelT);
		break;
	default:
		SolvePPEHO(NoOfElements, CC, CD, CDO, Re, DelT, Node);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		UpdatePressure(NoOfElements, CC, CD, CDO);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		SolveCFRNSHO(NoOfElements, CC, CD, CDO, Node, Re, DelT);
		ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Count);
		SolveNSHO(NoOfElements, CC, CD, CDO, Node, Re, DelT, Count);
		SolveTempHO(NoOfElements, CC, CD, CDO, Node, Re, DelT);
		break;
	}
	UpdateTemperature(NoOfElements, CC, CD);
	MassBalance(NoOfElements,CC, CD);
}
void ComputeNodalDatas(int NoOfNodes, CELLDATAS *CDO, NODE *Node)
{
	int Nod;
	int Count;
	double SumXiByDi, SumYiByDi, SumTiByDi, Sum1ByDi;
	for(Nod=0;Nod<NoOfNodes;Nod++)
	{
		SumXiByDi = 0.0;
		SumYiByDi = 0.0;
		SumTiByDi = 0.0;
		Sum1ByDi = 0.0;
		
		for(Count = 0;Count <Node[Nod].ShareCount; Count++)
		{
			SumXiByDi = SumXiByDi + (CDO[Node[Nod].ShareElements[Count]].CCVel[0] / Node[Nod].CC2NodeDist[Count]);
			SumYiByDi = SumYiByDi + (CDO[Node[Nod].ShareElements[Count]].CCVel[1] / Node[Nod].CC2NodeDist[Count]);
			SumTiByDi = SumTiByDi + (CDO[Node[Nod].ShareElements[Count]].CCT / Node[Nod].CC2NodeDist[Count]);
			Sum1ByDi = Sum1ByDi + (1/Node[Nod].CC2NodeDist[Count]);
		}
		Node[Nod].Sum = Sum1ByDi;
		Node[Nod].u = SumXiByDi / Sum1ByDi;
		Node[Nod].v = SumYiByDi / Sum1ByDi;
		Node[Nod].T = SumTiByDi / Sum1ByDi;
	}
}
void ApplyBC(int NoOfElements, int NoOfNodes, int NoOfBodies, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, BODY_ELEM *BE, int **BodyElems, int Count)
{
	int Nod, Elem;
	int Body;
	int Body_Elems;
	double Xb, Yb;
	double XB, YB, XA, YA;
	double u,v;
	double Velocity;
	for(Elem=0;Elem<NoOfElements; Elem++)
	{	
		
		switch(CC[Elem].Type)
		{
		case BODY1:
			/*No Slip Boundary Conditions
			CD[Elem].FCVel[0][0] = 0.0;
			CD[Elem].FCVel[0][1] = 0.0;
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = CD[Elem].CCT;
			
			CDO[Elem].FCVel[0][0] = 0.0;
			CDO[Elem].FCVel[0][1] = 0.0;
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = CDO[Elem].CCT;
			*/
			//Rotating Cylinder BC Assuming 0,0 as the centre of the cylinder Not a general case
			
		/*	Velocity = ALPHA * cos (2.0 * PI * FREQ * (DELT * Count));
			Xb = (Node[CC[Elem].Connect[0]-1].x + Node[CC[Elem].Connect[1]-1].x)/2.0;
			Yb = (Node[CC[Elem].Connect[0]-1].y + Node[CC[Elem].Connect[1]-1].y)/2.0;
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			CD[Elem].FCVel[0][0] = u;
			CD[Elem].FCVel[0][1] = v;
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = CD[Elem].CCT;
			
			CDO[Elem].FCVel[0][0] = u;
			CDO[Elem].FCVel[0][1] = v;
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = CDO[Elem].CCT;*/
			break;
		case INLET:
			CD[Elem].FCVel[0][0] = Uinf;
			CD[Elem].FCVel[0][1] = Vinf;
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = Tinit;

			CDO[Elem].FCVel[0][0] = Uinf;
			CDO[Elem].FCVel[0][1] = Vinf;
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = Tinit;
			break;
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			break;
		case OUTLET:
			//CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];
			//CD[Elem].FCVel[0][1] = CD[Elem].CCVel[1];
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = CD[Elem].CCT;
			
			//CDO[Elem].FCVel[0][0] = CDO[Elem].CCVel[0];
			//CDO[Elem].FCVel[0][1] = CDO[Elem].CCVel[1];
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = CDO[Elem].CCT;
			break;
		case TOP:
		case BOTTOM:
			CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];
			CD[Elem].FCVel[0][1] = Vinf;
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = CD[Elem].CCT;
			//CD[Elem].FCT[0] = Tinit;
			
			CDO[Elem].FCVel[0][0] = CDO[Elem].CCVel[0];
			CDO[Elem].FCVel[0][1] = Vinf;
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = CDO[Elem].CCT;
			//CDO[Elem].FCT[0] = Tinit;
			break;
		}
	}
	for(Body=0;Body<NoOfBodies;Body++)
	{
		for(Body_Elems=0; Body_Elems<BE[Body].MaxElems;Body_Elems++)
		{
			Elem = BodyElems[Body][Body_Elems];
			XB = (Node[CC[Elem].Connect[0]-1].x + Node[CC[Elem].Connect[1]-1].x)/2.0;
			YB = (Node[CC[Elem].Connect[0]-1].y + Node[CC[Elem].Connect[1]-1].y)/2.0;
			XA = (Node[BE[Body].LeadingEdge -1].x + Node[BE[Body].TrailingEdge -1].x)/2.0;
			YA = (Node[BE[Body].LeadingEdge -1].y + Node[BE[Body].TrailingEdge -1].y)/2.0;
			Xb = XB - XA;
			Yb = YB - YA;
			Velocity = BE[Body].Amplitude * cos (2.0 * PI * BE[Body].Frequency * (DELT * Count));
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			CD[Elem].FCVel[0][0] = u;
			CD[Elem].FCVel[0][1] = v;
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCT[0] = CD[Elem].CCT;
		
			Velocity = BE[Body].Amplitude * cos (2.0 * PI * BE[Body].Frequency * (DELT * (Count-1)));	
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			CDO[Elem].FCVel[0][0] = u;
			CDO[Elem].FCVel[0][1] = v;
			CDO[Elem].FCP[0] = CDO[Elem].CCP;
			CDO[Elem].FCT[0] = CDO[Elem].CCT;
		
			Nod = CC[Elem].Connect[0] -1;
			XB = Node[Nod].x;
			YB = Node[Nod].y;
			Xb = XB - XA;
			Yb = YB - YA;
			Velocity = BE[Body].Amplitude * cos (2.0 * PI * BE[Body].Frequency * (DELT * Count));
			
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
				
			Node[Nod].u = u;
			Node[Nod].v = v;
			Node[Nod].T = TSource;
			
			Nod = CC[Elem].Connect[1] -1;
			XB = Node[Nod].x;
			YB = Node[Nod].y;
			Xb = XB - XA;
			Yb = YB - YA;
			Velocity = BE[Body].Amplitude * cos (2.0 * PI * BE[Body].Frequency * (DELT * Count));
			
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
				
			Node[Nod].u = u;
			Node[Nod].v = v;
			Node[Nod].T = TSource;
		}
		
	}
	
	
	for(Nod=0;Nod<NoOfNodes;Nod++)
	{
		int MinElem;
		int j;
		switch(Node[Nod].Type)
		{
		case BODY1:
		/* No slip BC
			Node[Nod].u = 0.0;
			Node[Nod].v = 0.0;
			Node[Nod].T = TSource;
		*/
		//Rotating Cylinder
		/*	Xb = Node[Nod].x;
			Yb = Node[Nod].y;
			Velocity = ALPHA * cos (2.0 * PI * FREQ * (DELT * Count));
			
			u = Yb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			v = -1.0 * Xb * Velocity / sqrt((Yb * Yb) + (Xb * Xb));
			
			Node[Nod].u = u;
			Node[Nod].v = v;
			Node[Nod].T = TSource;*/
			break;
		case TOP:
		
		
		case BOTTOM:
			MinElem = Node[Nod].ShareElements[0];
			for(j=1;j<Node[Nod].ShareCount;j++)
			{
				if(fabs(CC[MinElem].Xp - Node[Nod].x) > fabs(CC[Node[Nod].ShareElements[j]].Xp - Node[Nod].x))
				{
					MinElem = Node[Nod].ShareElements[j];
				}
			}
			Node[Nod].u = CD[MinElem].CCVel[0];
			Node[Nod].v = Vinf;
			Node[Nod].T = CD[MinElem].CCT;
			break;
		case INLET:
			Node[Nod].u = Uinf;
			Node[Nod].v = Vinf;
			Node[Nod].T = Tinf;
			break;
		case OUTLET:
			MinElem = Node[Nod].ShareElements[0];
			for(j=1;j<Node[Nod].ShareCount;j++)
			{
				if(fabs(CC[MinElem].Yp - Node[Nod].y) > fabs(CC[Node[Nod].ShareElements[j]].Yp - Node[Nod].y))
				{
					MinElem = Node[Nod].ShareElements[j];
				}
			}
			//Node[Nod].u = CD[MinElem].CCVel[0];
			//Node[Nod].v = CD[MinElem].CCVel[1];
			Node[Nod].T = CD[MinElem].CCT;
			break;
		}
	}
	
}
///////////////////////
void CalculateSource(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node)
{
	STENCIL S;
	int Elem;
	int j, MinElem;
	int Nod;
	for(Elem=0; Elem<NoOfElements;Elem++)
	{
		CD[Elem].Source = 0.0;
		
		switch(CC[Elem].Type)
		{
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			S = UpdateStencil(CDO, Elem, CC, Node);
			
			CD[Elem].XDFLUX1 = S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca + S.dudx12 * CC[Elem].Yag - S.dudy12 * CC[Elem].Xag;
			CD[Elem].XDFLUX1 += S.dudx11 * CC[Elem].Ygb - S.dudy11 * CC[Elem].Xgb + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc;
		
			CD[Elem].YDFLUX1 = S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca + S.dvdx12 * CC[Elem].Yag - S.dvdy12 * CC[Elem].Xag;
			CD[Elem].YDFLUX1 += S.dvdx11 * CC[Elem].Ygb - S.dvdy11 * CC[Elem].Xgb + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc;
			
			CD[Elem].XCFLUX1 = pow(S.u3,2) * CC[Elem].Yca - S.u3 * S.v3 *CC[Elem].Xca + pow(S.u12,2) * CC[Elem].Yag - S.u12 * S.v12 * CC[Elem].Xag;
			CD[Elem].XCFLUX1 += pow(S.u11,2) * CC[Elem].Ygb - S.u11 * S.v11 * CC[Elem].Xgb + pow(S.u2,2) * CC[Elem].Ybc - S.u2 * S.v2 * CC[Elem].Xbc;
		
			CD[Elem].YCFLUX1 = S.u3 * S.v3 * CC[Elem].Yca - pow(S.v3,2)* CC[Elem].Xca + S.u12 * S.v12 * CC[Elem].Yag - pow(S.v12,2) * CC[Elem].Xag;
			CD[Elem].YCFLUX1 += S.u11 * S.v11 * CC[Elem].Ygb - pow(S.v11,2) * CC[Elem].Xgb + S.u2 * S.v2 * CC[Elem].Ybc - pow(S.v2,2) * CC[Elem].Xbc;
				
			break;
		case OUTLET:
			Nod = CC[Elem].Connect[2];
			Nod = Nod -1;

			
			if(Elem == Node[Nod].ShareElements[0])
			{
				MinElem = Node[Nod].ShareElements[1];
				j = 2;
			}
			else
			{
				MinElem = Node[Nod].ShareElements[0];
				j = 1;
			}
			for(;j<Node[Nod].ShareCount;j++)
			{
				if(fabs(CC[MinElem].Yp - CC[Elem].Yp) > fabs(CC[Node[Nod].ShareElements[j]].Yp - CC[Elem].Yp))
				{
					if(Elem != Node[Nod].ShareElements[j])
						MinElem = Node[Nod].ShareElements[j];
				}
			}
			
			CD[Elem].XDFLUX1 = (Node[Nod].u - CDO[MinElem].CCVel[0]) / (Node[Nod].x - CC[MinElem].Xp);
			CD[Elem].YDFLUX1 = (Node[Nod].v - CDO[MinElem].CCVel[1]) / (Node[Nod].x - CC[MinElem].Xp);
			CD[Elem].XCFLUX1 = 0.0;
			CD[Elem].YCFLUX1 = 0.0;
			break;

		case BODY1:
		case INLET:

		case TOP:
		case BOTTOM:
			S = UpdateStencil(CDO, Elem, CC, Node);
			
			CD[Elem].XCFLUX1 = 0.0;
			CD[Elem].YCFLUX1 = 0.0;
			CD[Elem].XDFLUX1 = 0.0;
			CD[Elem].YDFLUX1 = 0.0;
			break;
		}
			
		CD[Elem].XDFLUX2 = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx5 * CC[Elem].Ybd - S.dudy5 * CC[Elem].Xbd;
		CD[Elem].XDFLUX2+= S.dudx4 * CC[Elem].Ydc - S.dudy4 * CC[Elem].Xdc + S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca;
	
		CD[Elem].YDFLUX2 = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx5 * CC[Elem].Ybd - S.dvdy5 * CC[Elem].Xbd;
		CD[Elem].YDFLUX2+= S.dvdx4 * CC[Elem].Ydc - S.dvdy4 * CC[Elem].Xdc + S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca;
		
		CD[Elem].XDFLUX3 = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc;
		CD[Elem].XDFLUX3+= S.dudx19 * CC[Elem].Ycj - S.dudy19 * CC[Elem].Xcj + S.dudx18 * CC[Elem].Yja - S.dudy18 * CC[Elem].Xja;

		CD[Elem].YDFLUX3 = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc;
		CD[Elem].YDFLUX3+= S.dvdx19 * CC[Elem].Ycj - S.dvdy19 * CC[Elem].Xcj + S.dvdx18 * CC[Elem].Yja - S.dvdy18 * CC[Elem].Xja;
	
		CD[Elem].XCFLUX2 = pow(S.u1,2) * CC[Elem].Yab - S.u1 * S.v1 * CC[Elem].Xab + pow(S.u5,2) * CC[Elem].Ybd - S.u5 * S.v5 * CC[Elem].Xbd;
		CD[Elem].XCFLUX2+= pow(S.u4,2) * CC[Elem].Ydc - S.u4 * S.v4 * CC[Elem].Xdc + pow(S.u3,2) * CC[Elem].Yca - S.u3 * S.v3 * CC[Elem].Xca;
		
		CD[Elem].YCFLUX2 = S.u1 * S.v1 * CC[Elem].Yab - pow(S.v1,2) * CC[Elem].Xab + S.u5 * S.v5 * CC[Elem].Ybd - pow(S.v5,2) * CC[Elem].Xbd;
		CD[Elem].YCFLUX2+= S.u4 * S.v4 * CC[Elem].Ydc - pow(S.v4,2) * CC[Elem].Xdc + S.u3 * S.v3 * CC[Elem].Yca - pow(S.v3,2) * CC[Elem].Xca;
	
		CD[Elem].XCFLUX3 = pow(S.u1,2) * CC[Elem].Yab - S.u1 * S.v1 * CC[Elem].Xab + pow(S.u2,2) * CC[Elem].Ybc - S.u2 * S.v2 * CC[Elem].Xbc;
		CD[Elem].XCFLUX3+= pow(S.u19,2) * CC[Elem].Ycj - S.u19 * S.v19 * CC[Elem].Xcj + pow(S.u18,2) * CC[Elem].Yja - S.u18 * S.v18 * CC[Elem].Xja;
		
		CD[Elem].YCFLUX3 = S.u1 * S.v1 * CC[Elem].Yab - pow(S.v1,2) * CC[Elem].Xab + S.u2 * S.v2 * CC[Elem].Ybc - pow(S.v2,2) * CC[Elem].Xbc;
		CD[Elem].YCFLUX3+= S.u19 * S.v19 * CC[Elem].Ycj - pow(S.v19,2) * CC[Elem].Xcj + S.u18 * S.v18 * CC[Elem].Yja - pow(S.v18,2) * CC[Elem].Xja;

		CD[Elem].Source = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca - S.v1 * CC[Elem].Xab - S.v2 * CC[Elem].Xbc - S.v3 * CC[Elem].Xca)/ DelT;
				
		if(CC[Elem].Type == OUTLET)
		{
			CD[Elem].Source += CD[Elem].YDFLUX1 * CC[Elem].Xab / CC[Elem].Area;
			CD[Elem].Source -= CD[Elem].XDFLUX1 * CC[Elem].Yab / CC[Elem].Area;
		}
		else
		{
			CD[Elem].Source += ( (CD[Elem].XDFLUX1/Re) - (CD[Elem].XCFLUX1) ) * (CC[Elem].Yab / CC[Elem].Apc);
			CD[Elem].Source -= ( (CD[Elem].YDFLUX1/Re)- (CD[Elem].YCFLUX1) ) * (CC[Elem].Xab / CC[Elem].Apc);
		}

		CD[Elem].Source += ( (CD[Elem].XDFLUX2/Re) - (CD[Elem].XCFLUX2) ) * (CC[Elem].Ybc / CC[Elem].Apa);
		CD[Elem].Source += ( (CD[Elem].XDFLUX3/Re) - (CD[Elem].XCFLUX3) ) * (CC[Elem].Yca / CC[Elem].Apb);
		CD[Elem].Source -= ( (CD[Elem].YDFLUX2/Re)- (CD[Elem].YCFLUX2) ) * (CC[Elem].Xbc / CC[Elem].Apa);
		CD[Elem].Source -= ( (CD[Elem].YDFLUX3/Re)- (CD[Elem].YCFLUX3) ) * (CC[Elem].Xca / CC[Elem].Apb);
	}
	
}

void CalculateSourceHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node)
{
	STENCIL S;
	int Elem;
	int j, MinElem;
	int Nod;
	double NX1, N_1X1;
	double NX2, N_1X2;
	double NX3, N_1X3;

	double NY1, N_1Y1;
	double NY2, N_1Y2;
	double NY3, N_1Y3;
	//FILE *TMP;
	//TMP = fopen("Source.Csv","w+");
	for(Elem=0; Elem<NoOfElements;Elem++)
	{
		CD[Elem].Source = 0.0;
		
		switch(CC[Elem].Type)
		{
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:

			S = UpdateStencil(CDO, Elem, CC, Node);
			
			CD[Elem].XDFLUX1 = S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca + S.dudx12 * CC[Elem].Yag - S.dudy12 * CC[Elem].Xag;
			CD[Elem].XDFLUX1 += S.dudx11 * CC[Elem].Ygb - S.dudy11 * CC[Elem].Xgb + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc;
		
			CD[Elem].YDFLUX1 = S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca + S.dvdx12 * CC[Elem].Yag - S.dvdy12 * CC[Elem].Xag;
			CD[Elem].YDFLUX1 += S.dvdx11 * CC[Elem].Ygb - S.dvdy11 * CC[Elem].Xgb + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc;
			
			CD[Elem].XCFLUX1 = pow(S.u3,2) * CC[Elem].Yca - S.u3 * S.v3 *CC[Elem].Xca + pow(S.u12,2) * CC[Elem].Yag - S.u12 * S.v12 * CC[Elem].Xag;
			CD[Elem].XCFLUX1 += pow(S.u11,2) * CC[Elem].Ygb - S.u11 * S.v11 * CC[Elem].Xgb + pow(S.u2,2) * CC[Elem].Ybc - S.u2 * S.v2 * CC[Elem].Xbc;
		
			CD[Elem].YCFLUX1 = S.u3 * S.v3 * CC[Elem].Yca - pow(S.v3,2)* CC[Elem].Xca + S.u12 * S.v12 * CC[Elem].Yag - pow(S.v12,2) * CC[Elem].Xag;
			CD[Elem].YCFLUX1 += S.u11 * S.v11 * CC[Elem].Ygb - pow(S.v11,2) * CC[Elem].Xgb + S.u2 * S.v2 * CC[Elem].Ybc - pow(S.v2,2) * CC[Elem].Xbc;
			
			N_1X1 = (CDO[Elem].XDFLUX1/Re) - CDO[Elem].XCFLUX1 - CDO[Elem].XPFLUX1;
			N_1Y1 = (CDO[Elem].YDFLUX1/Re) - CDO[Elem].YCFLUX1 - CDO[Elem].YPFLUX1;
			
			break;
		case OUTLET:
			Nod = CC[Elem].Connect[2];
			Nod = Nod -1;

			if(Elem == Node[Nod].ShareElements[0])
			{
				MinElem = Node[Nod].ShareElements[1];
				j = 2;
			}
			else
			{
				MinElem = Node[Nod].ShareElements[0];
				j = 1;
			}
			for(;j<Node[Nod].ShareCount;j++)
			{
				if(fabs(CC[MinElem].Yp - CC[Elem].Yp) > fabs(CC[Node[Nod].ShareElements[j]].Yp - CC[Elem].Yp))
				{
					if(Elem != Node[Nod].ShareElements[j])
						MinElem = Node[Nod].ShareElements[j];
				}
			}
			
			CD[Elem].XDFLUX1 = (CDO[Elem].CCVel[0] - CDO[MinElem].CCVel[0]) / (CC[Elem].Xp - CC[MinElem].Xp);
			CD[Elem].YDFLUX1 = (CDO[Elem].CCVel[1] - CDO[MinElem].CCVel[1]) / (CC[Elem].Xp - CC[MinElem].Xp);
			CD[Elem].XCFLUX1 = 0.0;
			CD[Elem].YCFLUX1 = 0.0;
			N_1X1 = CDO[Elem].XDFLUX1;
			N_1Y1 = CDO[Elem].YDFLUX1;
			break;

		case BODY1:
		case INLET:
		case TOP:
		case BOTTOM:

			S = UpdateStencil(CDO, Elem, CC, Node);
			
			CD[Elem].XCFLUX1 = 0.0;
			CD[Elem].YCFLUX1 = 0.0;
			CD[Elem].XDFLUX1 = 0.0;
			CD[Elem].YDFLUX1 = 0.0;

			N_1X1 = 0.0;
			N_1Y1 = 0.0;
			break;
		
		}
			
		CD[Elem].XDFLUX2 = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx5 * CC[Elem].Ybd - S.dudy5 * CC[Elem].Xbd;
		CD[Elem].XDFLUX2+= S.dudx4 * CC[Elem].Ydc - S.dudy4 * CC[Elem].Xdc + S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca;
	
		CD[Elem].YDFLUX2 = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx5 * CC[Elem].Ybd - S.dvdy5 * CC[Elem].Xbd;
		CD[Elem].YDFLUX2+= S.dvdx4 * CC[Elem].Ydc - S.dvdy4 * CC[Elem].Xdc + S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca;
		
		CD[Elem].XDFLUX3 = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc;
		CD[Elem].XDFLUX3+= S.dudx19 * CC[Elem].Ycj - S.dudy19 * CC[Elem].Xcj + S.dudx18 * CC[Elem].Yja - S.dudy18 * CC[Elem].Xja;

		CD[Elem].YDFLUX3 = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc;
		CD[Elem].YDFLUX3+= S.dvdx19 * CC[Elem].Ycj - S.dvdy19 * CC[Elem].Xcj + S.dvdx18 * CC[Elem].Yja - S.dvdy18 * CC[Elem].Xja;
	
		CD[Elem].XCFLUX2 = pow(S.u1,2) * CC[Elem].Yab - S.u1 * S.v1 * CC[Elem].Xab + pow(S.u5,2) * CC[Elem].Ybd - S.u5 * S.v5 * CC[Elem].Xbd;
		CD[Elem].XCFLUX2+= pow(S.u4,2) * CC[Elem].Ydc - S.u4 * S.v4 * CC[Elem].Xdc + pow(S.u3,2) * CC[Elem].Yca - S.u3 * S.v3 * CC[Elem].Xca;
		
		CD[Elem].YCFLUX2 = S.u1 * S.v1 * CC[Elem].Yab - pow(S.v1,2) * CC[Elem].Xab + S.u5 * S.v5 * CC[Elem].Ybd - pow(S.v5,2) * CC[Elem].Xbd;
		CD[Elem].YCFLUX2+= S.u4 * S.v4 * CC[Elem].Ydc - pow(S.v4,2) * CC[Elem].Xdc + S.u3 * S.v3 * CC[Elem].Yca - pow(S.v3,2) * CC[Elem].Xca;
	
		CD[Elem].XCFLUX3 = pow(S.u1,2) * CC[Elem].Yab - S.u1 * S.v1 * CC[Elem].Xab + pow(S.u2,2) * CC[Elem].Ybc - S.u2 * S.v2 * CC[Elem].Xbc;
		CD[Elem].XCFLUX3+= pow(S.u19,2) * CC[Elem].Ycj - S.u19 * S.v19 * CC[Elem].Xcj + pow(S.u18,2) * CC[Elem].Yja - S.u18 * S.v18 * CC[Elem].Xja;
		
		CD[Elem].YCFLUX3 = S.u1 * S.v1 * CC[Elem].Yab - pow(S.v1,2) * CC[Elem].Xab + S.u2 * S.v2 * CC[Elem].Ybc - pow(S.v2,2) * CC[Elem].Xbc;
		CD[Elem].YCFLUX3+= S.u19 * S.v19 * CC[Elem].Ycj - pow(S.v19,2) * CC[Elem].Xcj + S.u18 * S.v18 * CC[Elem].Yja - pow(S.v18,2) * CC[Elem].Xja;

		CD[Elem].Source = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca - S.v1 * CC[Elem].Xab - S.v2 * CC[Elem].Xbc - S.v3 * CC[Elem].Xca)/ DelT;
		if(CC[Elem].Type == OUTLET)
		{
			NX1   = CD[Elem].XDFLUX1;
			NY1   = CD[Elem].YDFLUX1;
			CD[Elem].Source += ( 3.0 * NX1 - N_1X1 ) * (CC[Elem].Yab / (2.0 * CC[Elem].Area));
			CD[Elem].Source -= ( 3.0 * NY1 - N_1Y1 ) * (CC[Elem].Xab / (2.0 * CC[Elem].Area));
		}
		else
		{
			NX1   = (CD[Elem].XDFLUX1/Re) - CD[Elem].XCFLUX1;
			NY1   = (CD[Elem].YDFLUX1/Re) - CD[Elem].YCFLUX1;
			CD[Elem].Source += ( 3.0 * NX1 - N_1X1 ) * (CC[Elem].Yab / (2.0 * CC[Elem].Apc));
			CD[Elem].Source -= ( 3.0 * NY1 - N_1Y1 ) * (CC[Elem].Xab / (2.0 * CC[Elem].Apc));
		}
		
		NX2   = (CD[Elem].XDFLUX2/Re) - CD[Elem].XCFLUX2;
		N_1X2 = (CDO[Elem].XDFLUX2/Re) - CDO[Elem].XCFLUX2 - CDO[Elem].XPFLUX2;
		NX3   = (CD[Elem].XDFLUX3/Re) - CD[Elem].XCFLUX3;
		N_1X3 = (CDO[Elem].XDFLUX3/Re) - CDO[Elem].XCFLUX3 - CDO[Elem].XPFLUX3;
		
		NY2   = (CD[Elem].YDFLUX2/Re) - CD[Elem].YCFLUX2;
		N_1Y2 = (CDO[Elem].YDFLUX2/Re) - CDO[Elem].YCFLUX2 - CDO[Elem].YPFLUX2;
		NY3   = (CD[Elem].YDFLUX3/Re) - CD[Elem].YCFLUX3;
		N_1Y3 = (CDO[Elem].YDFLUX3/Re) - CDO[Elem].YCFLUX3 - CDO[Elem].YPFLUX3;

		CD[Elem].Source += ( 3.0 * NX2 - N_1X2 ) * (CC[Elem].Ybc / (2.0 * CC[Elem].Apa));
		CD[Elem].Source += ( 3.0 * NX3 - N_1X3 ) * (CC[Elem].Yca / (2.0 *CC[Elem].Apb));

		CD[Elem].Source -= ( 3.0 * NY2 - N_1Y2 ) * (CC[Elem].Xbc / (2.0 * CC[Elem].Apa));
		CD[Elem].Source -= ( 3.0 * NY3 - N_1Y3 ) * (CC[Elem].Xca / (2.0 * CC[Elem].Apb));
		//fprintf(TMP,"%d,%e\n",Elem,CD[Elem].Source);
	}
	//fclose(TMP);
}
void SolvePPE(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node)
{
	int Elem, Count, MaxDivLocation=0;
	double MaxDiv, Div= 0.0;
	//int tmp;
	Count = 0;
	MaxDiv = 0.0;
	
	CalculateSource(NoOfElements, CC, CD, CDO, Re, DelT, Node);
	do
	{
		MaxDiv = 0.0;
		for(Elem=0; Elem<NoOfElements; Elem++)
		{
			switch(CC[Elem].Type)
			{
			case DOM:
			case BODY1:
			case BODYNNUN:
			case INLET:
			case INNNUN:
			//case OUTLET:
			case OUTNNUN:
			case TOP:
			case BOTTOM:
			case TOPNNUN:
			case BOTNNUN:
				Div = CalculateP(Elem, CC, CD, CDO, Re, DelT,Count);

				if(Div < 0.0)
				Div *= -1.0;

				if(MaxDiv < Div)
				{
					MaxDiv = Div;
					MaxDivLocation = Elem;
				}
				/*printf("\nElem %d Type %d Div %e",Elem+1, CC[Elem].Type, Div);
				if(CC[Elem].Type ==  INLET || CC[Elem].Type == INNNUN)
					getch();*/
				break;
			}
		}
					
		for(Elem=0;Elem<NoOfElements;Elem++)
			CDO[Elem].CCP = CD[Elem].CCP;
		//printf("\nCount = %d, MaxDiv = %e, MaxDiv Location = %d, Location Type = %d", Count, MaxDiv, MaxDivLocation, CC[MaxDivLocation].Type);
		
		Count++;
		
		
	}while( (Count <= MAXPPECOUNT) && (MaxDiv >= MAXDIV) );
	printf(", Count = %d, MaxDiv = %.3e", Count, MaxDiv);
	//printf(", Count = %d, MaxDiv = %.3e, MaxDiv Location = %d, Location Type = %d", Count, MaxDiv, MaxDivLocation, CC[MaxDivLocation].Type);
	//printf(", Count = %d, ", Count);
}	

void SolvePPEHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, NODE *Node)
{
	int Elem, Count, MaxDivLocation=0;
	double MaxDiv, Div= 0.0;
	//int tmp;
	Count = 0;
	MaxDiv = 0.0;
	
	CalculateSourceHO(NoOfElements, CC, CD, CDO, Re, DelT, Node);
	do
	{
		MaxDiv = 0.0;
		for(Elem=0; Elem<NoOfElements; Elem++)
		{
			switch(CC[Elem].Type)
			{
			case DOM:
			case BODY1:
			case BODYNNUN:
			case INLET:
			case INNNUN:
			//case OUTLET:
			case OUTNNUN:
			case TOP:
			case BOTTOM:
			case TOPNNUN:
			case BOTNNUN:
				Div = CalculatePHO(Elem, CC, CD, CDO, Re, DelT,Count);

				if(Div < 0.0)
				Div *= -1.0;

				if(MaxDiv < Div)
				{
					MaxDiv = Div;
					MaxDivLocation = Elem;
				}
				/*printf("\nElem %d Type %d Div %e",Elem+1, CC[Elem].Type, Div);
				if(CC[Elem].Type ==  INLET || CC[Elem].Type == INNNUN)
					getch();*/
				break;
			}
		}
					
		for(Elem=0;Elem<NoOfElements;Elem++)
			CDO[Elem].CCP = CD[Elem].CCP;
		//printf("\nCount = %d, MaxDiv = %e, MaxDiv Location = %d, Location Type = %d", Count, MaxDiv, MaxDivLocation, CC[MaxDivLocation].Type);
		
		Count++;
		
		
	}while( (Count <= MAXPPECOUNT) && (MaxDiv >= MAXDIV) );
	printf(", Count = %d, MaxDiv = %.3e", Count, MaxDiv);
	//printf(", Count = %d, MaxDiv = %.3e, MaxDiv Location = %d, Location Type = %d", Count, MaxDiv, MaxDivLocation, CC[MaxDivLocation].Type);
	//printf(", Count = %d, ", Count);
}

void UpdateDerivatives(int NoOfElements, CELLDETAILS *CC, CELLDATAS* C, NODE *Node)
{
///////////////////////////////////////
	int Elem;
	double uPrime, vPrime;
	int j, MinElem, Nod;
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		
		
			switch(CC[Elem].Type)
			{
			case DOM:
			case BODYNNUN:
			case INNNUN:
			case OUTNNUN:
			case TOPNNUN:
			case BOTNNUN:
				C[Elem].FCdu_dx[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].YCP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[2]].CCVel[0]) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdu_dx[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].YAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdu_dx[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].YBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdv_dx[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].YCP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[2]].CCVel[1]) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdv_dx[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].YAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdv_dx[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].YBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdu_dy[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].XCP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[2]].CCVel[0]) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdu_dy[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].XAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdu_dy[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].XBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				C[Elem].FCdv_dy[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].XCP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[2]].CCVel[1]) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdv_dy[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].XAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdv_dy[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].XBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				C[Elem].FCdT_dx[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].YCP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[2]].CCT) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdT_dx[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].YAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdT_dx[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].YBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdT_dy[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].XCP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[2]].CCT) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdT_dy[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].XAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdT_dy[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].XBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				break;
			case INLET:
				C[Elem].FCdu_dx[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].YCP) - ((C[Elem].CCVel[0] - Uinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdu_dx[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].YAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdu_dx[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].YBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdv_dx[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].YCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdv_dx[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].YAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdv_dx[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].YBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdu_dy[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].XCP) - ((C[Elem].CCVel[0] - Uinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdu_dy[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].XAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdu_dy[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].XBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				C[Elem].FCdv_dy[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].XCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdv_dy[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].XAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdv_dy[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].XBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				
				C[Elem].FCdT_dx[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].YCP) - ((C[Elem].CCT - Tinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdT_dx[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].YAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdT_dx[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].YBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdT_dy[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].XCP) - ((C[Elem].CCT - Tinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdT_dy[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].XAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdT_dy[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].XBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				break;
			case TOP:
			case BOTTOM:
				C[Elem].FCdu_dx[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].YCP) - ((C[Elem].CCVel[0] - C[Elem].CCVel[0]) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdu_dx[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].YAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdu_dx[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].YBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdv_dx[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].YCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdv_dx[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].YAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdv_dx[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].YBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				//C[Elem].FCdu_dy[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].XCP) - ((C[Elem].CCVel[0] - Uinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdu_dy[0] = 0.0;
				C[Elem].FCdu_dy[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].XAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdu_dy[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].XBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				//C[Elem].FCdv_dy[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].XCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdv_dy[0] = 0.0;
				C[Elem].FCdv_dy[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].XAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdv_dy[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].XBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				
				C[Elem].FCdT_dx[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].YCP) - ((C[Elem].CCT - C[Elem].CCT) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdT_dx[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].YAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdT_dx[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].YBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdT_dy[0] = 0.0;
				C[Elem].FCdT_dy[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].XAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdT_dy[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].XBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				break;
			case OUTLET:
				Nod = CC[Elem].Connect[2];
				Nod = Nod -1;

			
				if(Elem == Node[Nod].ShareElements[0])
				{
					MinElem = Node[Nod].ShareElements[1];
					j = 2;
				}
				else
				{
					MinElem = Node[Nod].ShareElements[0];
					j = 1;
				}
				for(;j<Node[Nod].ShareCount;j++)
				{
					if(fabs(CC[MinElem].Yp - CC[Elem].Yp) > fabs(CC[Node[Nod].ShareElements[j]].Yp - CC[Elem].Yp))
					{
						if(Elem != Node[Nod].ShareElements[j])
							MinElem = Node[Nod].ShareElements[j];
					}
				}
			
				//CD[Elem].XDFlux1 = 
				//CD[Elem].YDFlux1 = 
				
				//C[Elem].FCdu_dx[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].YCP) - ((C[Elem].CCVel[0] - Uinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				//C[Elem].FCdu_dx[0] = 0.0;
				C[Elem].FCdu_dx[0] = (C[Elem].CCVel[0] - C[MinElem].CCVel[0]) / (CC[Elem].Xp - CC[MinElem].Xp);
				C[Elem].FCdu_dx[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].YAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdu_dx[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].YBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				//C[Elem].FCdv_dx[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].YCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				//C[Elem].FCdv_dx[0] = 0.0;
				C[Elem].FCdv_dx[0] = (C[Elem].CCVel[1] - C[MinElem].CCVel[1]) / (CC[Elem].Xp - CC[MinElem].Xp);
				C[Elem].FCdv_dx[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].YAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdv_dx[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].YBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				//C[Elem].FCdu_dy[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].XCP) - ((C[Elem].CCVel[0] - Uinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdu_dy[0] = (Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u ) / CC[Elem].Yab;
				C[Elem].FCdu_dy[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].XAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdu_dy[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].XBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				//C[Elem].FCdv_dy[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].XCP) - ((C[Elem].CCVel[1] - Vinf) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdv_dy[0] = (Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v ) / CC[Elem].Yab;
				C[Elem].FCdv_dy[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].XAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdv_dy[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].XBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				
				C[Elem].FCdT_dx[0] = 0.0;
				C[Elem].FCdT_dx[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].YAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdT_dx[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].YBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);
				
				C[Elem].FCdT_dy[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].XCP) - ((C[Elem].CCT - C[Elem].CCT) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdT_dy[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].XAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdT_dy[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].XBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				break;
			
			case BODY1:
				uPrime = 2.0 * C[Elem].FCVel[0][0] - C[Elem].CCVel[0];
				vPrime = 2.0 * C[Elem].FCVel[0][1] - C[Elem].CCVel[1];
				C[Elem].FCdu_dx[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].YCP) - ((C[Elem].CCVel[0] - uPrime) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdu_dx[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].YAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdu_dx[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].YBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdv_dx[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].YCP) - ((C[Elem].CCVel[1] - vPrime) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdv_dx[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].YAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdv_dx[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].YBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdu_dy[0] = (((Node[CC[Elem].Connect[1]-1].u - Node[CC[Elem].Connect[0]-1].u) * CC[Elem].XCP) - ((C[Elem].CCVel[0] - uPrime) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdu_dy[1] = (((Node[CC[Elem].Connect[2]-1].u - Node[CC[Elem].Connect[1]-1].u) * CC[Elem].XAP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdu_dy[2] = (((Node[CC[Elem].Connect[0]-1].u - Node[CC[Elem].Connect[2]-1].u) * CC[Elem].XBP) - ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);

				C[Elem].FCdv_dy[0] = (((Node[CC[Elem].Connect[1]-1].v - Node[CC[Elem].Connect[0]-1].v) * CC[Elem].XCP) - ((C[Elem].CCVel[1] - vPrime) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdv_dy[1] = (((Node[CC[Elem].Connect[2]-1].v - Node[CC[Elem].Connect[1]-1].v) * CC[Elem].XAP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1]) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdv_dy[2] = (((Node[CC[Elem].Connect[0]-1].v - Node[CC[Elem].Connect[2]-1].v) * CC[Elem].XBP) - ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1]) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				
				C[Elem].FCdT_dx[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].YCP) - ((C[Elem].CCT - TSource) * CC[Elem].Yab))/(CC[Elem].Xab * CC[Elem].YCP - CC[Elem].XCP * CC[Elem].Yab);
				C[Elem].FCdT_dx[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].YAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Ybc))/(CC[Elem].Xbc * CC[Elem].YAP - CC[Elem].XAP * CC[Elem].Ybc);
				C[Elem].FCdT_dx[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].YBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Yca))/(CC[Elem].Xca * CC[Elem].YBP - CC[Elem].XBP * CC[Elem].Yca);

				C[Elem].FCdT_dy[0] = (((Node[CC[Elem].Connect[1]-1].T - Node[CC[Elem].Connect[0]-1].T) * CC[Elem].XCP) - ((C[Elem].CCT - TSource) * CC[Elem].Xab))/(CC[Elem].XCP * CC[Elem].Yab - CC[Elem].Xab * CC[Elem].YCP);
				C[Elem].FCdT_dy[1] = (((Node[CC[Elem].Connect[2]-1].T - Node[CC[Elem].Connect[1]-1].T) * CC[Elem].XAP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[0]].CCT) * CC[Elem].Xbc))/(CC[Elem].XAP * CC[Elem].Ybc - CC[Elem].Xbc * CC[Elem].YAP);
				C[Elem].FCdT_dy[2] = (((Node[CC[Elem].Connect[0]-1].T - Node[CC[Elem].Connect[2]-1].T) * CC[Elem].XBP) - ((C[Elem].CCT - C[CC[Elem].Neighbour[1]].CCT) * CC[Elem].Xca))/(CC[Elem].XBP * CC[Elem].Yca - CC[Elem].Xca * CC[Elem].YBP);
				break;
			default:
				printf("\n\nCase not written in update derivatives. \npress any key to exit...");
				exit(9);
			}
	}
}


STENCIL UpdateStencil(CELLDATAS *C, int Elem, CELLDETAILS *CC, NODE *N)
{
	STENCIL S;
	
	S.u1 = C[Elem].FCVel[0][0];
	S.u2 = C[Elem].FCVel[1][0];
	S.u3 = C[Elem].FCVel[2][0];

	S.v1 = C[Elem].FCVel[0][1];
	S.v2 = C[Elem].FCVel[1][1];
	S.v3 = C[Elem].FCVel[2][1];

	S.p1 = C[Elem].FCP [0];
	S.p2 = C[Elem].FCP [1];
	S.p3 = C[Elem].FCP [2];

	S.dudx1 = C[Elem].FCdu_dx [0];
	S.dudx2 = C[Elem].FCdu_dx [1];
	S.dudx3 = C[Elem].FCdu_dx [2];

	S.dudy1 = C[Elem].FCdu_dy [0];
	S.dudy2 = C[Elem].FCdu_dy [1];
	S.dudy3 = C[Elem].FCdu_dy [2];

	S.dvdx1 = C[Elem].FCdv_dx [0];
	S.dvdx2 = C[Elem].FCdv_dx [1];
	S.dvdx3 = C[Elem].FCdv_dx [2];

	S.dvdy1 = C[Elem].FCdv_dy [0];
	S.dvdy2 = C[Elem].FCdv_dy [1];
	S.dvdy3 = C[Elem].FCdv_dy [2];
	
	S.u4 = C[CC[Elem].Neighbour[0]].FCVel[CC[Elem].Index4][0];
	S.u5 = C[CC[Elem].Neighbour[0]].FCVel[CC[Elem].Index5][0];
	S.u18 = C[CC[Elem].Neighbour[1]].FCVel[CC[Elem].Index18][0];
	S.u19 = C[CC[Elem].Neighbour[1]].FCVel[CC[Elem].Index19][0];
	
	S.v4 = C[CC[Elem].Neighbour[0]].FCVel[CC[Elem].Index4][1];
	S.v5 = C[CC[Elem].Neighbour[0]].FCVel[CC[Elem].Index5][1];
	S.v18 = C[CC[Elem].Neighbour[1]].FCVel[CC[Elem].Index18][1];
	S.v19 = C[CC[Elem].Neighbour[1]].FCVel[CC[Elem].Index19][1];


	S.p4 = C[CC[Elem].Neighbour[0]].FCP[CC[Elem].Index4];
	S.p5 = C[CC[Elem].Neighbour[0]].FCP[CC[Elem].Index5];
	S.p18 = C[CC[Elem].Neighbour[1]].FCP[CC[Elem].Index18];
	S.p19 = C[CC[Elem].Neighbour[1]].FCP[CC[Elem].Index19];

	S.dudx4 = C[CC[Elem].Neighbour[0]].FCdu_dx[CC[Elem].Index4];
	S.dudx5 = C[CC[Elem].Neighbour[0]].FCdu_dx[CC[Elem].Index5];
	S.dudx18 = C[CC[Elem].Neighbour[1]].FCdu_dx[CC[Elem].Index18];
	S.dudx19 = C[CC[Elem].Neighbour[1]].FCdu_dx[CC[Elem].Index19];

	S.dvdx4 = C[CC[Elem].Neighbour[0]].FCdv_dx[CC[Elem].Index4];
	S.dvdx5 = C[CC[Elem].Neighbour[0]].FCdv_dx[CC[Elem].Index5];
	S.dvdx18 = C[CC[Elem].Neighbour[1]].FCdv_dx[CC[Elem].Index18];
	S.dvdx19 = C[CC[Elem].Neighbour[1]].FCdv_dx[CC[Elem].Index19];

	S.dudy4 = C[CC[Elem].Neighbour[0]].FCdu_dy[CC[Elem].Index4];
	S.dudy5 = C[CC[Elem].Neighbour[0]].FCdu_dy[CC[Elem].Index5];
	S.dudy18 = C[CC[Elem].Neighbour[1]].FCdu_dy[CC[Elem].Index18];
	S.dudy19 = C[CC[Elem].Neighbour[1]].FCdu_dy[CC[Elem].Index19];

	S.dvdy4 = C[CC[Elem].Neighbour[0]].FCdv_dy[CC[Elem].Index4];
	S.dvdy5 = C[CC[Elem].Neighbour[0]].FCdv_dy[CC[Elem].Index5];
	S.dvdy18 = C[CC[Elem].Neighbour[1]].FCdv_dy[CC[Elem].Index18];
	S.dvdy19 = C[CC[Elem].Neighbour[1]].FCdv_dy[CC[Elem].Index19];


	switch(CC[Elem].Type)
	{
	case DOM:
	case BODYNNUN:
	case TOPNNUN:
	case BOTNNUN:
	case INNNUN:
	case OUTNNUN:
		S.u11 = C[CC[Elem].Neighbour[2]].FCVel[CC[Elem].Index11][0];
		S.u12 = C[CC[Elem].Neighbour[2]].FCVel[CC[Elem].Index12][0];
		S.v11 = C[CC[Elem].Neighbour[2]].FCVel[CC[Elem].Index11][1];
		S.v12 = C[CC[Elem].Neighbour[2]].FCVel[CC[Elem].Index12][1];
		S.p11 = C[CC[Elem].Neighbour[2]].FCP[CC[Elem].Index11];
		S.p12 = C[CC[Elem].Neighbour[2]].FCP[CC[Elem].Index12];
		S.dudx11 = C[CC[Elem].Neighbour[2]].FCdu_dx[CC[Elem].Index11];
		S.dudx12 = C[CC[Elem].Neighbour[2]].FCdu_dx[CC[Elem].Index12];
		S.dvdx11 = C[CC[Elem].Neighbour[2]].FCdv_dx[CC[Elem].Index11];
		S.dvdx12 = C[CC[Elem].Neighbour[2]].FCdv_dx[CC[Elem].Index12];
		S.dudy11 = C[CC[Elem].Neighbour[2]].FCdu_dy[CC[Elem].Index11];
		S.dudy12 = C[CC[Elem].Neighbour[2]].FCdu_dy[CC[Elem].Index12];
		S.dvdy11 = C[CC[Elem].Neighbour[2]].FCdv_dy[CC[Elem].Index11];
		S.dvdy12 = C[CC[Elem].Neighbour[2]].FCdv_dy[CC[Elem].Index12];
			break;
	case BODY1:
		S.u11 = -1.0 * S.u2;
		S.u12 = -1.0 * S.u3;
		S.v11 = -1.0 * S.v2;
		S.v12 = -1.0 * S.v3;
			break;
	case OUTLET:
		/*S.u11 = S.u2;
		S.u12 = S.u3;
		S.v11 = S.v2;
		S.v12 = S.v3;
		S.p11 = S.p2;
		S.p12 = S.p3;
		S.dudx11 = 0.0;
		S.dudx12 = 0.0;
		S.dvdx11 = 0.0;
		S.dvdx12 = 0.0;
		S.dudy11 = ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0])*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(C[Elem].CCVel[0] - N[CC[Elem].Connect[1]-1].u)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[0]].Yp)*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(N[CC[Elem].Connect[2]-1].y - N[CC[Elem].Connect[1]-1].y)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp));
		S.dudy12 = ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0])*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].u - C[Elem].CCVel[0])*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[1]].Yp)*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].y - N[CC[Elem].Connect[2]-1].y)*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp));
		S.dvdy11 = ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1])*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(C[Elem].CCVel[1] - N[CC[Elem].Connect[1]-1].v)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[0]].Yp)*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(N[CC[Elem].Connect[2]-1].y - N[CC[Elem].Connect[1]-1].y)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp));
		S.dvdy12 = ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1])*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].v - C[Elem].CCVel[1])*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[1]].Yp)*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].y - N[CC[Elem].Connect[2]-1].y)*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp));
		*/
		break;
	case TOP:
	case BOTTOM:
	case INLET:
	
		S.u11 = S.u2;
		S.u12 = S.u3;
		S.v11 = 0.0;
		S.v12 = 0.0;
		S.p11 = S.p2;
		S.p12 = S.p3;
		/*S.dudx11 = 0.0;
		S.dudx12 = 0.0;
		S.dvdx11 = 0.0;
		S.dvdx12 = 0.0;
		S.dudy11 = ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[0]].CCVel[0])*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(C[Elem].CCVel[0] - N[CC[Elem].Connect[1]-1].u)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[0]].Yp)*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(N[CC[Elem].Connect[2]-1].y - N[CC[Elem].Connect[1]-1].y)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp));
		S.dudy12 = ((C[Elem].CCVel[0] - C[CC[Elem].Neighbour[1]].CCVel[0])*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].u - C[Elem].CCVel[0])*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[1]].Yp)*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].y - N[CC[Elem].Connect[2]-1].y)*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp));
		S.dvdy11 = ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[0]].CCVel[1])*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(C[Elem].CCVel[1] - N[CC[Elem].Connect[1]-1].v)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[0]].Yp)*(N[CC[Elem].Connect[0]-1].x - N[CC[Elem].Connect[2]-1].x)-(N[CC[Elem].Connect[2]-1].y - N[CC[Elem].Connect[1]-1].y)*(CC[CC[Elem].Neighbour[0]].Xp - CC[Elem].Xp));
		S.dvdy12 = ((C[Elem].CCVel[1] - C[CC[Elem].Neighbour[1]].CCVel[1])*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].v - C[Elem].CCVel[1])*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp))/((CC[Elem].Yp - CC[CC[Elem].Neighbour[1]].Yp)*(N[CC[Elem].Connect[2]-1].x - N[CC[Elem].Connect[0]-1].x)-(N[CC[Elem].Connect[0]-1].y - N[CC[Elem].Connect[2]-1].y)*(CC[CC[Elem].Neighbour[1]].Xp - CC[Elem].Xp));
		*/
			break; 
	}
	return S;
}


double CalculateP(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT,int Count)
{
	double PrevP;
	double Div;
	int i;
	i=Elem+1;
	PrevP = CD[Elem].CCP;
	
	switch(CC[Elem].Type)
	{
	case DOM:
		CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP - CC[Elem].Kg * CDO[CC[Elem].CellG].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		break;
	case BODYNNUN:
	case INNNUN:
	case OUTNNUN:
	case TOPNNUN:
	case BOTNNUN:
		if(CC[Elem].CellG == UNAVAILABLEDATA && CC[Elem].CellI == UNAVAILABLEDATA)
		{
			CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		if(CC[Elem].CellG != UNAVAILABLEDATA)
		{
			CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP - CC[Elem].Kg * CDO[CC[Elem].CellG].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		if(CC[Elem].CellI != UNAVAILABLEDATA)
		{
			CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		break;
	case BODY1:
		CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;
		break;
	case INLET:
	case TOP:
	case BOTTOM:
	//case OUTLET:
		CD[Elem].CCP = CD[Elem].Source - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;
		break;
	}
	
	CD[Elem].CCP = CD[Elem].CCP / CC[Elem].Kp;
	Div = PrevP - CD[Elem].CCP;
	return Div;
}

double CalculatePHO(int Elem, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT,int Count)
{
	double PrevP;
	double Div;
	int i;
	i=Elem+1;
	PrevP = CD[Elem].CCP;
	
	switch(CC[Elem].Type)
	{
	case DOM:
		CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP - CC[Elem].Kg * CDO[CC[Elem].CellG].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		break;
	case BODYNNUN:
	case INNNUN:
	case OUTNNUN:
	case TOPNNUN:
	case BOTNNUN:
		if(CC[Elem].CellG == UNAVAILABLEDATA && CC[Elem].CellI == UNAVAILABLEDATA)
		{
			CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		if(CC[Elem].CellG != UNAVAILABLEDATA)
		{
			CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP - CC[Elem].Kg * CDO[CC[Elem].CellG].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		if(CC[Elem].CellI != UNAVAILABLEDATA)
		{
			CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kc * CDO[CC[Elem].Neighbour [2]].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP - CC[Elem].Kf * CDO[CC[Elem].CellF].CCP;
			CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;	
		}
		break;
	case BODY1:
		CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;
		break;
	case INLET:
	case TOP:
	case BOTTOM:
	//case OUTLET:
		CD[Elem].CCP = (2.0 * CD[Elem].Source / 3.0) - CC[Elem].Ka * CDO[CC[Elem].Neighbour [0]].CCP - CC[Elem].Kb * CDO[CC[Elem].Neighbour [1]].CCP - CC[Elem].Kd * CDO[CC[Elem].CellD].CCP;
		CD[Elem].CCP = CD[Elem].CCP - CC[Elem].Ki * CDO[CC[Elem].CellI].CCP - CC[Elem].Kj * CDO[CC[Elem].CellJ].CCP - CC[Elem].Kl * CDO[CC[Elem].CellL].CCP;
		break;
	}
	
	CD[Elem].CCP = CD[Elem].CCP / CC[Elem].Kp;
	Div = PrevP - CD[Elem].CCP;
	return Div;
}
void UpdatePressure(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO)
{
	int Elem;
	//FILE *TMP;
	//TMP = fopen("p.csv","w+");
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		//fprintf(TMP,"%d,%e\n",Elem,CD[Elem].CCP);	
		switch(CC[Elem].Type)
		{
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			CD[Elem].FCP[0] = (CD[Elem].CCP * CC[CC[Elem].Neighbour[2]].Area + CD[CC[Elem].Neighbour[2]].CCP * CC[Elem].Area )/ CC[Elem].Apc;	
			CD[Elem].FCP[1] = (CD[Elem].CCP * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCP * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCP[2] = (CD[Elem].CCP * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCP * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case BODY1:
		case INLET:
		case OUTLET:
		case TOP:
		case BOTTOM:
			CD[Elem].FCP[0] = CD[Elem].CCP;
			CD[Elem].FCP[1] = (CD[Elem].CCP * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCP * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCP[2] = (CD[Elem].CCP * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCP * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		}	
	}
	//fclose(TMP);
}

void UpdateTemperature(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD)
{
	int Elem;
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		switch(CC[Elem].Type)
		{
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			CD[Elem].FCT[0] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[2]].Area + CD[CC[Elem].Neighbour[2]].CCT * CC[Elem].Area )/ CC[Elem].Apc;	
			CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case BODY1:
			CD[Elem].FCT[0] = TSource;
			CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case INLET:
			CD[Elem].FCT[0] = Tinf;
			CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case OUTLET:
		case TOP:
		case BOTTOM:
			CD[Elem].FCT[0] = CD[Elem].CCT;
			CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		}	
		Elem=Elem;
	}
}

void SolveCFRNS(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT)
{
	int Elem;
	//STENCIL S, SN;
	STENCIL SN;
	double MaxDiv = 0.0, Div = 0.0;
	
	for(Elem = 0; Elem <NoOfElements; Elem++)
	{			
		switch(CC[Elem].Type)
		{
		case OUTNNUN:
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case TOPNNUN: 
		case BOTNNUN:
			//S = UpdateStencil(CDO, Elem, CC, Node);
			SN = UpdateStencil(CD, Elem, CC, Node);
			
			CD[Elem].XPFLUX1 = SN.p3 * CC[Elem].Yca + SN.p12 * CC[Elem].Yag + SN.p11 * CC[Elem].Ygb + SN.p2 * CC[Elem].Ybc;
			CD[Elem].YPFLUX1 = -1.0 * (SN.p3 * CC[Elem].Xca + SN.p12 * CC[Elem].Xag + SN.p11 * CC[Elem].Xgb + SN.p2 * CC[Elem].Xbc);
			
			CD[Elem].XPFLUX2 = SN.p1 * CC[Elem].Yab + SN.p5 * CC[Elem].Ybd + SN.p4 * CC[Elem].Ydc + SN.p3 * CC[Elem].Yca;
			CD[Elem].YPFLUX2 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p5 * CC[Elem].Xbd + SN.p4 * CC[Elem].Xdc + SN.p3 * CC[Elem].Xca);
			
			CD[Elem].XPFLUX3 = SN.p1 * CC[Elem].Yab + SN.p2 * CC[Elem].Ybc + SN.p19 * CC[Elem].Ycj + SN.p18 * CC[Elem].Yja;
			CD[Elem].YPFLUX3 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p2 * CC[Elem].Xbc + SN.p19 * CC[Elem].Xcj + SN.p18 * CC[Elem].Xja);

			
			CD[Elem].FCVel[1][0] = CDO[Elem].FCVel[1][0] + DelT * (CD[Elem].XDFLUX2/Re - CD[Elem].XCFLUX2 - CD[Elem].XPFLUX2)/CC[Elem].Apa;
			Div = CDO[Elem].FCVel[1][0] - CD[Elem].FCVel[1][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[1][1] = CDO[Elem].FCVel[1][1] + DelT * (CD[Elem].YDFLUX2/Re - CD[Elem].YCFLUX2 - CD[Elem].YPFLUX2)/CC[Elem].Apa;
			Div = CDO[Elem].FCVel[1][1] - CD[Elem].FCVel[1][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][0] = CDO[Elem].FCVel[2][0] + DelT * (CD[Elem].XDFLUX3/Re - CD[Elem].XCFLUX3 - CD[Elem].XPFLUX3)/CC[Elem].Apb;
			Div = CDO[Elem].FCVel[2][0] - CD[Elem].FCVel[2][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][1] = CDO[Elem].FCVel[2][1] + DelT * (CD[Elem].YDFLUX3/Re - CD[Elem].YCFLUX3 - CD[Elem].YPFLUX3)/CC[Elem].Apb;
			Div = CDO[Elem].FCVel[2][1] - CD[Elem].FCVel[2][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[0][0] = CDO[Elem].FCVel[0][0] + DelT * (( (CD[Elem].XDFLUX1/Re) - CD[Elem].XCFLUX1 - CD[Elem].XPFLUX1 )/CC[Elem].Apc);
			Div = CDO[Elem].FCVel[0][0] - CD[Elem].FCVel[0][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}
				
			CD[Elem].FCVel[0][1] = CDO[Elem].FCVel[0][1] + DelT * (((CD[Elem].YDFLUX1/Re) - CD[Elem].YCFLUX1 - CD[Elem].YPFLUX1)/CC[Elem].Apc);
			Div = CDO[Elem].FCVel[0][1] - CD[Elem].FCVel[0][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			break;
			
		case OUTLET:
		case BODY1:
		case INLET:
		case TOP:
		case BOTTOM:
		//case OUTLET:
			//S = UpdateStencil(CDO, Elem, CC, Node);
			SN = UpdateStencil(CD, Elem, CC, Node);
			
			if(CC[Elem].Type == INLET)
			{
				CD[Elem].FCVel[0][0] = Uinf;
				CD[Elem].FCVel[0][1] = Vinf;
			}
			else if(CC[Elem].Type == TOP || CC[Elem].Type == BOTTOM)
			{
				CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];
				CD[Elem].FCVel[0][1] = Vinf;
			}
			else if(CC[Elem].Type == OUTLET)
			{
				CD[Elem].FCVel[0][0] = CDO[Elem].CCVel[0] - (CD[Elem].XDFLUX1 * DelT);
				CD[Elem].FCVel[0][1] = CDO[Elem].CCVel[1] - (CD[Elem].YDFLUX1 * DelT);
			}
			else 
			{
				CD[Elem].FCVel[0][0] = 0.0;
				CD[Elem].FCVel[0][1] = 0.0;
			}
			
			CD[Elem].XPFLUX2 = SN.p1 * CC[Elem].Yab + SN.p5 * CC[Elem].Ybd + SN.p4 * CC[Elem].Ydc + SN.p3 * CC[Elem].Yca;
			CD[Elem].YPFLUX2 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p5 * CC[Elem].Xbd + SN.p4 * CC[Elem].Xdc + SN.p3 * CC[Elem].Xca);
			
			CD[Elem].XPFLUX3 = SN.p1 * CC[Elem].Yab + SN.p2 * CC[Elem].Ybc + SN.p19 * CC[Elem].Ycj + SN.p18 * CC[Elem].Yja;
			CD[Elem].YPFLUX3 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p2 * CC[Elem].Xbc + SN.p19 * CC[Elem].Xcj + SN.p18 * CC[Elem].Xja);

			
			CD[Elem].FCVel[1][0] = CDO[Elem].FCVel[1][0] + DelT * (CD[Elem].XDFLUX2/Re - CD[Elem].XCFLUX2 - CD[Elem].XPFLUX2)/CC[Elem].Apa;
			
			Div = CDO[Elem].FCVel[1][0] - CD[Elem].FCVel[1][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[1][1] = CDO[Elem].FCVel[1][1] + DelT * (CD[Elem].YDFLUX2/Re - CD[Elem].YCFLUX2 - CD[Elem].YPFLUX2)/CC[Elem].Apa;

			Div = CDO[Elem].FCVel[1][1] - CD[Elem].FCVel[1][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][0] = CDO[Elem].FCVel[2][0] + DelT * (CD[Elem].XDFLUX3/Re - CD[Elem].XCFLUX3 - CD[Elem].XPFLUX3)/CC[Elem].Apb;

			Div = CDO[Elem].FCVel[2][0] - CD[Elem].FCVel[2][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][1] = CDO[Elem].FCVel[2][1] + DelT * (CD[Elem].YDFLUX3/Re - CD[Elem].YCFLUX3 - CD[Elem].YPFLUX3)/CC[Elem].Apb;

			Div = CDO[Elem].FCVel[2][1] - CD[Elem].FCVel[2][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}
			break;
		}	
		
	}
	
}

void SolveNS(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int count)
{
	int Elem;
	double XDIFF, YDIFF, XCONV, YCONV, XPRESS, YPRESS;
	STENCIL S;

	UpdateDerivatives(NoOfElements, CC, CD, Node);
	Re = RE;
	DelT = DELT;
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		S = UpdateStencil(CD, Elem, CC, Node);
		
		CD[Elem].CCdudx = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca)/CC[Elem].Area;
		CD[Elem].CCdvdx = (S.v1 * CC[Elem].Yab + S.v2 * CC[Elem].Ybc + S.v3 * CC[Elem].Yca)/CC[Elem].Area;
		CD[Elem].CCdudy = -1 * (S.u1 * CC[Elem].Xab + S.u2 * CC[Elem].Xbc + S.u3 * CC[Elem].Xca)/CC[Elem].Area;
		CD[Elem].CCdvdy = -1 * (S.v1 * CC[Elem].Xab + S.v2 * CC[Elem].Xbc + S.v3 * CC[Elem].Xca)/CC[Elem].Area;
		CD[Elem].Vorticity = CD[Elem].CCdvdx - CD[Elem].CCdudy;
		CD[Elem].Divergence = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca) - (S.v1 * CC[Elem].Xab + S.v2 * CC[Elem].Xbc + S.v3 * CC[Elem].Xca);
		switch(CC[Elem].Type)
		{
		
		case INLET:
		case OUTLET:
		case BODY1:
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
		case TOP:
		case BOTTOM:
			XDIFF = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc + S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca;
			YDIFF = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc + S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca;

			XCONV = pow(S.u1,2)  * CC[Elem].Yab + pow(S.u2,2) * CC[Elem].Ybc + pow(S.u3,2) * CC[Elem].Yca - S.u1 * S.v1 * CC[Elem]. Xab - S.u2 * S.v2 * CC[Elem].Xbc - S.u3 * S.v3 * CC[Elem].Xca;
			YCONV = S.u1 * S.v1 * CC[Elem]. Yab + S.u2 * S.v2 * CC[Elem].Ybc + S.u3 * S.v3 * CC[Elem].Yca - pow(S.v1,2)  * CC[Elem].Xab - pow(S.v2,2) * CC[Elem].Xbc - pow(S.v3,2) * CC[Elem].Xca;

			XPRESS = S.p1 * CC[Elem].Yab + S.p2 * CC[Elem]. Ybc + S.p3 * CC[Elem].Yca;
			YPRESS = -1.0 * ( S.p1 * CC[Elem].Xab + S.p2 * CC[Elem]. Xbc + S.p3 * CC[Elem].Xca );

			CD[Elem].CCVel[0] = CDO[Elem].CCVel[0] + DelT * (XDIFF/Re - XCONV - XPRESS)/CC[Elem].Area;
			CD[Elem].CCVel[1] = CDO[Elem].CCVel[1] + DelT * (YDIFF/Re - YCONV - YPRESS)/CC[Elem].Area;
			break;
		}
		//fprintf(FP,"%d, %d, %e, %e, %e, %e\n",Elem+1,CC[Elem].Type, CD[Elem].CCVel[0],CD[Elem].CCVel[1], CD[Elem].CCP, CDO[Elem].Source);
	}	
	
}


void SolveCFRNSHO (int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT)
{
	int Elem;
	STENCIL S, SN;
//	FILE *FPu, *FPv;
	double NX1, N_1X1;
	double NX2, N_1X2;
	double NX3, N_1X3;

	double NY1, N_1Y1;
	double NY2, N_1Y2;
	double NY3, N_1Y3;
	double MaxDiv = 0.0, Div = 0.0;
	
	for(Elem = 0; Elem <NoOfElements; Elem++)
	{			
		NX1   = (CD[Elem].XDFLUX1/Re) - CD[Elem].XCFLUX1;
		N_1X1 = (CDO[Elem].XDFLUX1/Re) - CDO[Elem].XCFLUX1 - CDO[Elem].XPFLUX1;
		NX2   = (CD[Elem].XDFLUX2/Re) - CD[Elem].XCFLUX2;
		N_1X2 = (CDO[Elem].XDFLUX2/Re) - CDO[Elem].XCFLUX2 - CDO[Elem].XPFLUX2;
		NX3   = (CD[Elem].XDFLUX3/Re) - CD[Elem].XCFLUX3;
		N_1X3 = (CDO[Elem].XDFLUX3/Re) - CDO[Elem].XCFLUX3 - CDO[Elem].XPFLUX3;

		NY1   = (CD[Elem].YDFLUX1/Re) - CD[Elem].YCFLUX1;
		N_1Y1 = (CDO[Elem].YDFLUX1/Re) - CDO[Elem].YCFLUX1 - CDO[Elem].YPFLUX1;
		NY2   = (CD[Elem].YDFLUX2/Re) - CD[Elem].YCFLUX2;
		N_1Y2 = (CDO[Elem].YDFLUX2/Re) - CDO[Elem].YCFLUX2 - CDO[Elem].YPFLUX2;
		NY3   = (CD[Elem].YDFLUX3/Re) - CD[Elem].YCFLUX3;
		N_1Y3 = (CDO[Elem].YDFLUX3/Re) - CDO[Elem].YCFLUX3 - CDO[Elem].YPFLUX3;

		switch(CC[Elem].Type)
		{
		
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			S = UpdateStencil(CDO, Elem, CC, Node);
			SN = UpdateStencil(CD, Elem, CC, Node);

			
			CD[Elem].XPFLUX1 = SN.p3 * CC[Elem].Yca + SN.p12 * CC[Elem].Yag + SN.p11 * CC[Elem].Ygb + SN.p2 * CC[Elem].Ybc;
			CD[Elem].YPFLUX1 = -1.0 * (SN.p3 * CC[Elem].Xca + SN.p12 * CC[Elem].Xag + SN.p11 * CC[Elem].Xgb + SN.p2 * CC[Elem].Xbc);
			
			CD[Elem].XPFLUX2 = SN.p1 * CC[Elem].Yab + SN.p5 * CC[Elem].Ybd + SN.p4 * CC[Elem].Ydc + SN.p3 * CC[Elem].Yca;
			CD[Elem].YPFLUX2 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p5 * CC[Elem].Xbd + SN.p4 * CC[Elem].Xdc + SN.p3 * CC[Elem].Xca);
			
			CD[Elem].XPFLUX3 = SN.p1 * CC[Elem].Yab + SN.p2 * CC[Elem].Ybc + SN.p19 * CC[Elem].Ycj + SN.p18 * CC[Elem].Yja;
			CD[Elem].YPFLUX3 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p2 * CC[Elem].Xbc + SN.p19 * CC[Elem].Xcj + SN.p18 * CC[Elem].Xja);

			
			CD[Elem].FCVel[1][0] = CDO[Elem].FCVel[1][0] + DelT * (3.0*(NX2 - CD[Elem].XPFLUX1)-N_1X2 )/(2.0 *CC[Elem].Apa);
			Div = CDO[Elem].FCVel[1][0] - CD[Elem].FCVel[1][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[1][1] = CDO[Elem].FCVel[1][1] + DelT * (3.0*(NY2 - CD[Elem].YPFLUX1)-N_1Y2 )/(2.0*CC[Elem].Apa);
			Div = CDO[Elem].FCVel[1][1] - CD[Elem].FCVel[1][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][0] = CDO[Elem].FCVel[2][0] + DelT * (3.0*(NX3 - CD[Elem].XPFLUX3)-N_1X3 )/(2.0*CC[Elem].Apb);
			Div = CDO[Elem].FCVel[2][0] - CD[Elem].FCVel[2][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][1] = CDO[Elem].FCVel[2][1] + DelT * (3.0*(NY3 - CD[Elem].YPFLUX3)-N_1Y3 )/(2.0*CC[Elem].Apb);
			Div = CDO[Elem].FCVel[2][1] - CD[Elem].FCVel[2][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[0][0] = CDO[Elem].FCVel[0][0] + DelT * (3.0*(NX1 - CD[Elem].XPFLUX1)-N_1X1 )/(2.0*CC[Elem].Apc);
			Div = CDO[Elem].FCVel[0][0] - CD[Elem].FCVel[0][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}
				
			CD[Elem].FCVel[0][1] = CDO[Elem].FCVel[0][1] + DelT * (3.0*(NY1 - CD[Elem].YPFLUX1)-N_1Y1 )/(2.0*CC[Elem].Apc);
			Div = CDO[Elem].FCVel[0][1] - CD[Elem].FCVel[0][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			break;
			
		case BODY1:
		case INLET:
		case OUTLET:
		case TOP:
		case BOTTOM:
			S = UpdateStencil(CDO, Elem, CC, Node);
			SN = UpdateStencil(CD, Elem, CC, Node);

			if(CC[Elem].Type == INLET)
			{
				CD[Elem].FCVel[0][0] = Uinf;
				CD[Elem].FCVel[0][1] = Vinf;
			}
			else if(CC[Elem].Type == OUTLET)
			{
				CD[Elem].FCVel[0][0] = CDO[Elem].CCVel[0] - (( 3.0 * CD[Elem].XDFLUX1 - CDO[Elem].XDFLUX1) * DelT / 2.0);
				CD[Elem].FCVel[0][1] = CDO[Elem].CCVel[1] - (( 3.0 * CD[Elem].YDFLUX1 - CDO[Elem].YDFLUX1) * DelT / 2.0);
				//CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];
				//CD[Elem].FCVel[0][1] = CD[Elem].CCVel[1];
			}
			else if(CC[Elem].Type == TOP || CC[Elem].Type == BOTTOM)
			{
				CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];
				CD[Elem].FCVel[0][1] = Vinf;
			}
			else
			{
				CD[Elem].FCVel[0][0] = 0.0;
				CD[Elem].FCVel[0][1] = 0.0;
			}
			
			CD[Elem].XPFLUX2 = SN.p1 * CC[Elem].Yab + SN.p5 * CC[Elem].Ybd + SN.p4 * CC[Elem].Ydc + SN.p3 * CC[Elem].Yca;
			CD[Elem].YPFLUX2 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p5 * CC[Elem].Xbd + SN.p4 * CC[Elem].Xdc + SN.p3 * CC[Elem].Xca);
			
			CD[Elem].XPFLUX3 = SN.p1 * CC[Elem].Yab + SN.p2 * CC[Elem].Ybc + SN.p19 * CC[Elem].Ycj + SN.p18 * CC[Elem].Yja;
			CD[Elem].YPFLUX3 = -1.0 * (SN.p1 * CC[Elem].Xab + SN.p2 * CC[Elem].Xbc + SN.p19 * CC[Elem].Xcj + SN.p18 * CC[Elem].Xja);

			
			CD[Elem].FCVel[1][0] = CDO[Elem].FCVel[1][0] + DelT * (3.0*(NX2 - CD[Elem].XPFLUX2)-N_1X2 )/(2.0*CC[Elem].Apa);
			
			Div = CDO[Elem].FCVel[1][0] - CD[Elem].FCVel[1][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[1][1] = CDO[Elem].FCVel[1][1] + DelT * (3.0*(NY2 - CD[Elem].YPFLUX2)-N_1Y2 )/(2.0*CC[Elem].Apa);

			Div = CDO[Elem].FCVel[1][1] - CD[Elem].FCVel[1][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][0] = CDO[Elem].FCVel[2][0] + DelT * (3.0*(NX3 - CD[Elem].XPFLUX3)-N_1X3 )/(2.0*CC[Elem].Apb);

			Div = CDO[Elem].FCVel[2][0] - CD[Elem].FCVel[2][0];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}

			CD[Elem].FCVel[2][1] = CDO[Elem].FCVel[2][1] + DelT * (3.0*(NY3 - CD[Elem].YPFLUX3)-N_1Y3 )/(2.0*CC[Elem].Apb);

			Div = CDO[Elem].FCVel[2][1] - CD[Elem].FCVel[2][1];
			if(Div < 0.0)
			Div *= -1.0;
			if(MaxDiv < Div)
			{
				MaxDiv = Div;
			}
			break;
		}	
		
	}
	
	//exit(0);
}

void SolveNSHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int count)
{
	int Elem;
	double XDIFF, YDIFF, XCONV, YCONV, XPRESS, YPRESS;
	double XDIFFO, YDIFFO, XCONVO, YCONVO, XPRESSO, YPRESSO;
	//STENCIL S;
	STENCIL S, SO;//int tmp;
//	FILE *Deepan;
	UpdateDerivatives(NoOfElements, CC, CD, Node);
	Re = RE;
	DelT = DELT;
	//FILE *FP;
	//FP=fopen("DeepanCellCent.csv","w+");
	//ApplyBC(NoOfElements, CC, CD);
	//ApplyBC(NoOfElements, CC, CDO);
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		S = UpdateStencil(CD, Elem, CC, Node);
		SO = UpdateStencil(CDO, Elem, CC, Node);
		
		CD[Elem].CCdudx = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca)/CC[Elem].Area;
		CD[Elem].CCdvdx = (S.v1 * CC[Elem].Yab + S.v2 * CC[Elem].Ybc + S.v3 * CC[Elem].Yca)/CC[Elem].Area;
		CD[Elem].CCdudy = -1 * (S.u1 * CC[Elem].Xab + S.u2 * CC[Elem].Xbc + S.u3 * CC[Elem].Xca)/CC[Elem].Area;
		CD[Elem].CCdvdy = -1 * (S.v1 * CC[Elem].Xab + S.v2 * CC[Elem].Xbc + S.v3 * CC[Elem].Xca)/CC[Elem].Area;
		CD[Elem].Vorticity = CD[Elem].CCdvdx - CD[Elem].CCdudy;
		CD[Elem].Divergence = (S.u1 * CC[Elem].Yab + S.u2 * CC[Elem].Ybc + S.u3 * CC[Elem].Yca) - (S.v1 * CC[Elem].Xab + S.v2 * CC[Elem].Xbc + S.v3 * CC[Elem].Xca);
		switch(CC[Elem].Type)
		{
		case TOPNNUN:
		case BOTNNUN:
		case TOP:
		case BOTTOM:
		case INLET:
		case BODY1:
		case OUTLET:
		case DOM:
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:
			XDIFF = S.dudx1 * CC[Elem].Yab - S.dudy1 * CC[Elem].Xab + S.dudx2 * CC[Elem].Ybc - S.dudy2 * CC[Elem].Xbc + S.dudx3 * CC[Elem].Yca - S.dudy3 * CC[Elem].Xca;
			YDIFF = S.dvdx1 * CC[Elem].Yab - S.dvdy1 * CC[Elem].Xab + S.dvdx2 * CC[Elem].Ybc - S.dvdy2 * CC[Elem].Xbc + S.dvdx3 * CC[Elem].Yca - S.dvdy3 * CC[Elem].Xca;

			XCONV = pow(S.u1,2)  * CC[Elem].Yab + pow(S.u2,2) * CC[Elem].Ybc + pow(S.u3,2) * CC[Elem].Yca - S.u1 * S.v1 * CC[Elem]. Xab - S.u2 * S.v2 * CC[Elem].Xbc - S.u3 * S.v3 * CC[Elem].Xca;
			YCONV = S.u1 * S.v1 * CC[Elem]. Yab + S.u2 * S.v2 * CC[Elem].Ybc + S.u3 * S.v3 * CC[Elem].Yca - pow(S.v1,2)  * CC[Elem].Xab - pow(S.v2,2) * CC[Elem].Xbc - pow(S.v3,2) * CC[Elem].Xca;

			XPRESS = S.p1 * CC[Elem].Yab + S.p2 * CC[Elem]. Ybc + S.p3 * CC[Elem].Yca;
			YPRESS = -1.0 * ( S.p1 * CC[Elem].Xab + S.p2 * CC[Elem]. Xbc + S.p3 * CC[Elem].Xca );

			XDIFFO = SO.dudx1 * CC[Elem].Yab - SO.dudy1 * CC[Elem].Xab + SO.dudx2 * CC[Elem].Ybc - SO.dudy2 * CC[Elem].Xbc + SO.dudx3 * CC[Elem].Yca - SO.dudy3 * CC[Elem].Xca;
			YDIFFO = SO.dvdx1 * CC[Elem].Yab - SO.dvdy1 * CC[Elem].Xab + SO.dvdx2 * CC[Elem].Ybc - SO.dvdy2 * CC[Elem].Xbc + SO.dvdx3 * CC[Elem].Yca - SO.dvdy3 * CC[Elem].Xca;

			XCONVO = pow(SO.u1,2)  * CC[Elem].Yab + pow(SO.u2,2) * CC[Elem].Ybc + pow(SO.u3,2) * CC[Elem].Yca - SO.u1 * SO.v1 * CC[Elem]. Xab - SO.u2 * SO.v2 * CC[Elem].Xbc - SO.u3 * SO.v3 * CC[Elem].Xca;
			YCONVO = SO.u1 * SO.v1 * CC[Elem]. Yab + SO.u2 * SO.v2 * CC[Elem].Ybc + SO.u3 * SO.v3 * CC[Elem].Yca - pow(SO.v1,2)  * CC[Elem].Xab - pow(SO.v2,2) * CC[Elem].Xbc - pow(SO.v3,2) * CC[Elem].Xca;

			XPRESSO = SO.p1 * CC[Elem].Yab + SO.p2 * CC[Elem]. Ybc + SO.p3 * CC[Elem].Yca;
			YPRESSO = -1.0 * ( SO.p1 * CC[Elem].Xab + SO.p2 * CC[Elem]. Xbc + SO.p3 * CC[Elem].Xca );

			CD[Elem].CCVel[0] = CDO[Elem].CCVel[0] + DelT * (3.0 * (XDIFF/Re - XCONV - XPRESS) - (XDIFFO/Re - XCONVO - XPRESSO))/(2.0 * CC[Elem].Area);
			CD[Elem].CCVel[1] = CDO[Elem].CCVel[1] + DelT * (3.0 * (YDIFF/Re - YCONV - YPRESS) - (YDIFFO/Re - YCONVO - YPRESSO))/(2.0 * CC[Elem].Area);
			break;
		}
		//fprintf(FP,"%d, %d, %e, %e, %e, %e\n",Elem+1,CC[Elem].Type, CD[Elem].CCVel[0],CD[Elem].CCVel[1], CD[Elem].CCP, CDO[Elem].Source);
	}	
	
}


void UpdateVelocity(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD)
{
	int Elem;
	
	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		
		switch(CC[Elem].Type)
		{
		case DOM:
		case BODYNNUN:	
		case INNNUN:
		case OUTNNUN:
		case TOPNNUN:
		case BOTNNUN:
			CD[Elem].FCVel[0][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[2]].Area + CD[CC[Elem].Neighbour[2]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apc;	
			CD[Elem].FCVel[0][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[2]].Area + CD[CC[Elem].Neighbour[2]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apc;	
			CD[Elem].FCVel[1][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[1][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[2][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apb;
			CD[Elem].FCVel[2][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case INLET:
			CD[Elem].FCVel[0][0] = Uinf;	
			CD[Elem].FCVel[0][1] = Vinf;	
			CD[Elem].FCVel[1][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[1][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[2][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apb;
			CD[Elem].FCVel[2][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case OUTLET:
			//CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];	
			//CD[Elem].FCVel[0][1] = CD[Elem].CCVel[1];	
			CD[Elem].FCVel[1][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[1][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[2][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apb;
			CD[Elem].FCVel[2][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case BODY1:
			CD[Elem].FCVel[0][0] = 0.0;	
			CD[Elem].FCVel[0][1] = 0.0;	
			CD[Elem].FCVel[1][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[1][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[2][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apb;
			CD[Elem].FCVel[2][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		case TOP:
		case BOTTOM:
			CD[Elem].FCVel[0][0] = CD[Elem].CCVel[0];	
			CD[Elem].FCVel[0][1] = Vinf;	
			CD[Elem].FCVel[1][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[1][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apa;
			CD[Elem].FCVel[2][0] = (CD[Elem].CCVel[0] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[0] * CC[Elem].Area )/ CC[Elem].Apb;
			CD[Elem].FCVel[2][1] = (CD[Elem].CCVel[1] * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCVel[1] * CC[Elem].Area )/ CC[Elem].Apb;
			break;
		}	
	}
	
}

void MassBalance(int NoOfElements,CELLDETAILS *CC, CELLDATAS *CD)
{
	int Elem;
	double MassIn, MassOut;
	double MassDefect;
	double OutletLength;
	MassIn = 0.0;
	MassOut = 0.0;

	OutletLength = 0.0;
	for(Elem=0;Elem<NoOfElements;Elem++)
	{
		switch(CC[Elem].Type)
		{
		case INLET:
			MassIn = MassIn + CD[Elem].FCVel[0][0] * fabs(CC[Elem].Yab);
			break;
		case OUTLET:
			MassOut = MassOut + CD[Elem].CCVel[0] * fabs(CC[Elem].Yab);
			//nOutlet++;
			OutletLength += fabs(CC[Elem].Yab);
			break;
		}
	}
	MassDefect = MassOut - MassIn;
	//printf(" \nMass In %.10lf, Mass Out %.10lf, Mass Defect %e\n",MassIn, MassOut, MassDefect);
	//printf(" Outlength %e",OutletLength);
	for(Elem=0;Elem<NoOfElements;Elem++)
	{
		switch(CC[Elem].Type)
		{
		case OUTLET:
			CD[Elem].CCVel[0] = CD[Elem].CCVel[0] - (MassDefect / OutletLength);
			break;
		}
	}

	MassIn = 0.0;
	MassOut = 0.0;

	for(Elem=0;Elem<NoOfElements;Elem++)
	{
		switch(CC[Elem].Type)
		{
		case INLET:
			MassIn = MassIn + CD[Elem].FCVel[0][0] * fabs(CC[Elem].Yab);
			break;
		case OUTLET:
			MassOut = MassOut + CD[Elem].CCVel[0] * fabs(CC[Elem].Yab);
			//nOutlet++;
			break;
		}
	}
	MassDefect = MassOut - MassIn;
	printf(", MassDefect %.2e",MassDefect);
}

void SolveTempHO(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS * CDO, NODE *Node, int Re, double DelT)
{
	int Elem;
	double DIFF, CONV;
	double TO, T;
	double u1, v1, t1, u2, v2, t2, u3, v3, t3;
	double dTdx1, dTdx2, dTdx3, dTdy1, dTdy2, dTdy3;
	double PrevT;
	double MaxDiffT = 0.0;

	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		PrevT = CD[Elem].CCT;
		u1 = CD[Elem].FCVel[0][0];
		u2 = CD[Elem].FCVel[1][0];
		u3 = CD[Elem].FCVel[2][0];

		v1 = CD[Elem].FCVel[0][1];
		v2 = CD[Elem].FCVel[1][1];
		v3 = CD[Elem].FCVel[2][1];

		t1 = CD[Elem].FCT[0];
		t2 = CD[Elem].FCT[1];
		t3 = CD[Elem].FCT[2];

		dTdx1 = CD[Elem].FCdT_dx[0];
		dTdx2 = CD[Elem].FCdT_dx[1];
		dTdx3 = CD[Elem].FCdT_dx[2];

		dTdy1 = CD[Elem].FCdT_dy[0];
		dTdy2 = CD[Elem].FCdT_dy[1];
		dTdy3 = CD[Elem].FCdT_dy[2];

		
		CONV = (u1 * t1 * CC[Elem].Yab + u2 * t2 * CC[Elem].Ybc + u3 * t3 * CC[Elem].Yca) - (v1 * t1 * CC[Elem].Xab + v2 * t2 * CC[Elem].Xbc + v3 * t3 * CC[Elem].Xca);
		DIFF = (dTdx1 * CC[Elem].Yab + dTdx2 * CC[Elem].Ybc + dTdx3 * CC[Elem].Yca) - (dTdy1 * CC[Elem].Xab + dTdy2 * CC[Elem].Xbc + dTdy3 * CC[Elem].Xca);
		T = ( DIFF /(Re * Pr) ) - CONV;
		
		u1 = CDO[Elem].FCVel[0][0];
		u2 = CDO[Elem].FCVel[1][0];
		u3 = CDO[Elem].FCVel[2][0];

		v1 = CDO[Elem].FCVel[0][1];
		v2 = CDO[Elem].FCVel[1][1];
		v3 = CDO[Elem].FCVel[2][1];

		t1 = CDO[Elem].FCT[0];
		t2 = CDO[Elem].FCT[1];
		t3 = CDO[Elem].FCT[2];

		dTdx1 = CDO[Elem].FCdT_dx[0];
		dTdx2 = CDO[Elem].FCdT_dx[1];
		dTdx3 = CDO[Elem].FCdT_dx[2];

		dTdy1 = CDO[Elem].FCdT_dy[0];
		dTdy2 = CDO[Elem].FCdT_dy[1];
		dTdy3 = CDO[Elem].FCdT_dy[2];
		
		CONV = (u1 * t1 * CC[Elem].Yab + u2 * t2 * CC[Elem].Ybc + u3 * t3 * CC[Elem].Yca) - (v1 * t1 * CC[Elem].Xab + v2 * t2 * CC[Elem].Xbc + v3 * t3 * CC[Elem].Xca);
		DIFF = (dTdx1 * CC[Elem].Yab + dTdx2 * CC[Elem].Ybc + dTdx3 * CC[Elem].Yca) - (dTdy1 * CC[Elem].Xab + dTdy2 * CC[Elem].Xbc + dTdy3 * CC[Elem].Xca);

		TO = ( DIFF /(Re * Pr) ) - CONV;

		CD[Elem].CCT = 	CDO[Elem].CCT + DelT * (3.0 * T - TO)/(2.0 * CC[Elem].Area);	
		
		if( MaxDiffT < fabs(PrevT - CD[Elem].CCT))
		{
			MaxDiffT = fabs(PrevT - CD[Elem].CCT);
		}
		
	}
	printf(", TempConv = %.3e",MaxDiffT);

}

void SolveTemp(int NoOfElements, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS * CDO, NODE *Node, int Re, double DelT)
{
	int Elem;
	double DIFF, CONV;
	double u1, v1, t1, u2, v2, t2, u3, v3, t3;
	double dTdx1, dTdx2, dTdx3, dTdy1, dTdy2, dTdy3;
	double PrevT;
	double MaxDiffT = 0.0;

	for(Elem = 0; Elem<NoOfElements; Elem++)
	{
		PrevT = CD[Elem].CCT;
		u1 = CD[Elem].FCVel[0][0];
		u2 = CD[Elem].FCVel[1][0];
		u3 = CD[Elem].FCVel[2][0];

		v1 = CD[Elem].FCVel[0][1];
		v2 = CD[Elem].FCVel[1][1];
		v3 = CD[Elem].FCVel[2][1];

		t1 = CD[Elem].FCT[0];
		t2 = CD[Elem].FCT[1];
		t3 = CD[Elem].FCT[2];
		/*printf("\n%d, %lf, %lf, %lf",Elem, t1, t2, t3);
		if(Elem%50==0)
		scanf("%lf",&dTdx1);*/

		dTdx1 = CD[Elem].FCdT_dx[0];
		dTdx2 = CD[Elem].FCdT_dx[1];
		dTdx3 = CD[Elem].FCdT_dx[2];

		dTdy1 = CD[Elem].FCdT_dy[0];
		dTdy2 = CD[Elem].FCdT_dy[1];
		dTdy3 = CD[Elem].FCdT_dy[2];

		
		CONV = (u1 * t1 * CC[Elem].Yab + u2 * t2 * CC[Elem].Ybc + u3 * t3 * CC[Elem].Yca) - (v1 * t1 * CC[Elem].Xab + v2 * t2 * CC[Elem].Xbc + v3 * t3 * CC[Elem].Xca);
		DIFF = (dTdx1 * CC[Elem].Yab + dTdx2 * CC[Elem].Ybc + dTdx3 * CC[Elem].Yca) - (dTdy1 * CC[Elem].Xab + dTdy2 * CC[Elem].Xbc + dTdy3 * CC[Elem].Xca);
		
		CD[Elem].CCT = 	CDO[Elem].CCT + DelT * ((DIFF/(Re * Pr)) - CONV)/CC[Elem].Area;
		
		if( MaxDiffT < fabs(PrevT - CD[Elem].CCT))
		{
			MaxDiffT = fabs(PrevT - CD[Elem].CCT);
		}
			
		/*printf("\n%d\t%d\t%lf",Elem,CC[Elem].Type,CD[Elem].CCT);
		if(Elem%400==0)
			getch();*/
	}
	printf(", TempConv = %.3e",MaxDiffT);	
}

