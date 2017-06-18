/*
//Main.c
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#ifdef _MSC_VER
	#define COMP 1
#else
	#ifdef __GNUC__
		#define COMP 2
	#else
		#ifdef __INTEL_COMPILER
			#define COMP 3
		#else
			#define COMP 0
		#endif		
	#endif
#endif


#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#include<time.h>
#include"Main.h"
#include"PreSolver.h"
#include"Solver.h"
#include"PostSolver.h"

#include"PreSolver.c"
#include"Solver.c"
#include"PostSolver.c"


int main()
{
	//FILE *FP;
	time_t start,end;
	double dif;
	FILE *FPressPoints;
	char CPressPoints[100];
	NODE *Node;
	CELLDETAILS *CC;
	CELLDATAS *CD,*CDO; //CD is current, CDO in old datas

	BODY_ELEM *BE;
	int NoOfBodies;
	int i,j;
	int NoOfElements,NoOfNodes;
	double DelT;
	int Re;
	int Iter;
	//int **BodyElems;
        int **BodyElems;
	
	time(&start);
	Re = RE;
	DelT = DELT;
	
	
////////////////////////////
//
	if((BodyElems=(int **)malloc(MAXBODIES*sizeof(int *))) == NULL)
	{
		printf("\nError");
		exit(0);
	}
        for (i=0;i<=MAXBODIES;i++)
        {   
        	if((BodyElems[i]=(int *)malloc(MAXELEMSONBODY*sizeof(int)))==NULL)
        	{
			printf("\nError");
			exit(0);
		}
        }
        
         
	GetNoOfElemAndNodes(&NoOfElements, &NoOfNodes);

	if((Node =(NODE*)malloc((NoOfNodes)*sizeof(NODE)))==NULL)
	{
		printf("Not Enough Memory for Allocating Nodes. \nProgram Terminated");
		exit(2);
	}
	for(Iter=0;Iter <NoOfNodes;Iter++)
	{
		Node[Iter].ShareElements[0] = Node[Iter].ShareElements[1] = Node[Iter].ShareElements[2] = Node[Iter].ShareElements[3] = Node[Iter].ShareElements[4] = -100;
		Node[Iter].ShareElements[5] = Node[Iter].ShareElements[6] = Node[Iter].ShareElements[7] = Node[Iter].ShareElements[8] = Node[Iter].ShareElements[9] = -100;
		Node[Iter].ShareCount = 0;
	}
	Iter = 0;
	if((CC =(CELLDETAILS*)malloc((NoOfElements)*sizeof(CELLDETAILS)))==NULL)
	{
		printf("Not Enough Memory for Allocating Nodes. \nProgram Terminated");
		exit(2);
	}
	if((CD =(CELLDATAS*)malloc((NoOfElements)*sizeof(CELLDATAS)))==NULL)
	{
		printf("Not Enough Memory for Allocating Nodes. \nProgram Terminated");
		exit(2);
	}
	if((CDO =(CELLDATAS*)malloc((NoOfElements)*sizeof(CELLDATAS)))==NULL)
	{
		printf("Not Enough Memory for Allocating Nodes. \nProgram Terminated");
		exit(2);
	}
////////////////////////////////

	
	PreSolver(NoOfElements, NoOfNodes, CC, Node, COMP); 
	
	switch (COMP)
	{
	case 1:
		sprintf(CPressPoints,"Mesh\\PressPoints.dat");
		break;
	case 2:
		sprintf(CPressPoints,".//Mesh//PressPoints.dat");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}
	if( (FPressPoints = fopen(CPressPoints,"r+")) == NULL)
	{
		printf("\nUnable to open Restart File\nProgram Terminated\n");
		exit(0);
	}
	fscanf(FPressPoints,"%d",&NoOfBodies);
	if((BE =(BODY_ELEM*)malloc((NoOfBodies)*sizeof(BODY_ELEM)))==NULL)
	{
		printf("Not Enough Memory for Allocating Nodes. \nProgram Terminated");
		exit(2);
	}	
	
	for(i=0;i<NoOfBodies;i++)
	{
		int Upper, Lower;
		int Elem;
		double Xc, Yc;
		fscanf(FPressPoints,"%d",&Upper);
		fscanf(FPressPoints,"%d",&Lower);
		
		BE[i].MaxElems = Upper + Lower;
		fscanf(FPressPoints,"%d",&Elem);
		
		Elem = Elem-1;
		BE[i].LeadingEdge = CC[Elem].Connect[0];
		printf("\n1 LE :%d, %lf",BE[i].LeadingEdge,Node[BE[i].LeadingEdge-1].x);		
		BodyElems[i][0]=Elem;
		
		for(j=1;j<Upper;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem -= 1;
			BodyElems[i][j] = Elem;
		}

		BE[i].TrailingEdge = CC[Elem].Connect[1];
		for(j=Upper;j<Upper+Lower;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			BodyElems[i][j] = Elem-1;
		}
		Xc = (Node[BE[i].LeadingEdge-1].x+Node[BE[i].TrailingEdge-1].x)/2.0;
		Yc = (Node[BE[i].LeadingEdge-1].y+Node[BE[i].TrailingEdge-1].y)/2.0;
		printf("\n\nCentre of %d of %d Bodies is %lf,%lf",i+1,NoOfBodies,Xc, Yc);
		
		printf("\nEnter the Amplitude of Rotation of %d Body:",i+1);
		scanf("%lf",&BE[i].Amplitude);
		printf("Enter the Frequency of rotational Oscilation of %d Body:",i+1);
		scanf("%lf",&BE[i].Frequency);
	}
	printf("\n\n");
	fclose(FPressPoints);	
	/*for(i=0;i<MAXELEMSONBODY;i++)
	{
		printf("\n%d",i+1);
		for(j=0;j< MAXBODIES;j++)
			printf("\t%d",BodyElems[j][i]);
	}
	exit(0);*/
	Iter = RestartCheck(NoOfElements, NoOfNodes, Re);
	if(Iter)
	{
		FILE *Restart;
		char FileRestart[100];
		int Elem;
		printf("\n%d Iter Restart File Found...",Iter);
		switch (COMP)
		{
		case 1:
			sprintf(FileRestart,"Restart\\CD.dat");
			break;
		case 2:
			sprintf(FileRestart,".//Restart//CD.dat");
			break;
		default:
			printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
			exit(0);
		}

		if( (Restart = fopen(FileRestart,"r+")) == NULL)
		{
			printf("\nUnable to open Restart File\nProgram Terminated\n");
			exit(0);
		}
		
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		for(Elem=0; Elem<NoOfElements; Elem++)
		{
			
			fscanf(Restart,"%lf\t",&CDO[Elem].CCP);
			fscanf(Restart,"%lf\t",&CDO[Elem].CCT);
			fscanf(Restart,"%lf\t",&CDO[Elem].CCVel[0]);
			fscanf(Restart,"%lf\t",&CDO[Elem].CCVel[1]);

			fscanf(Restart,"%lf\t",&CDO[Elem].XCFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].XCFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].XCFLUX3);
			fscanf(Restart,"%lf\t",&CDO[Elem].YCFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].YCFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].YCFLUX3);

			fscanf(Restart,"%lf\t",&CDO[Elem].XDFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].XDFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].XDFLUX3);
			fscanf(Restart,"%lf\t",&CDO[Elem].YDFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].YDFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].YDFLUX3);

			fscanf(Restart,"%lf\t",&CDO[Elem].XPFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].XPFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].XPFLUX3);
			fscanf(Restart,"%lf\t",&CDO[Elem].YPFLUX1);
			fscanf(Restart,"%lf\t",&CDO[Elem].YPFLUX2);
			fscanf(Restart,"%lf\t",&CDO[Elem].YPFLUX3);
		}
		for(Elem=0; Elem<NoOfElements; Elem++)
		{
			
			fscanf(Restart,"%lf\t",&CD[Elem].CCP);
			fscanf(Restart,"%lf\t",&CD[Elem].CCT);
			fscanf(Restart,"%lf\t",&CD[Elem].CCVel[0]);
			fscanf(Restart,"%lf\t",&CD[Elem].CCVel[1]);
			
		}

		for(Elem=0; Elem<NoOfElements; Elem++)
		{
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

					CDO[Elem].FCT[0] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[2]].Area + CDO[CC[Elem].Neighbour[2]].CCT * CC[Elem].Area )/ CC[Elem].Apc;	
					CDO[Elem].FCT[1] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CDO[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CDO[Elem].FCT[2] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CDO[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
					break;
				case BODY1:
					CD[Elem].FCT[0] = TSource;
					CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;

					CDO[Elem].FCT[0] = TSource;
					CDO[Elem].FCT[1] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CDO[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CDO[Elem].FCT[2] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CDO[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
					break;
				case INLET:
					CD[Elem].FCT[0] = Tinf;
					CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;

					CDO[Elem].FCT[0] = Tinf;
					CDO[Elem].FCT[1] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CDO[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CDO[Elem].FCT[2] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CDO[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
					break;
				case OUTLET:
				case TOP:
				case BOTTOM:
					CD[Elem].FCT[0] = CD[Elem].CCT;
					CD[Elem].FCT[1] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CD[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CD[Elem].FCT[2] = (CD[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CD[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;

					CDO[Elem].FCT[0] = CDO[Elem].CCT;
					CDO[Elem].FCT[1] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[0]].Area + CDO[CC[Elem].Neighbour[0]].CCT * CC[Elem].Area )/ CC[Elem].Apa;
					CDO[Elem].FCT[2] = (CDO[Elem].CCT * CC[CC[Elem].Neighbour[1]].Area + CDO[CC[Elem].Neighbour[1]].CCT * CC[Elem].Area )/ CC[Elem].Apb;
				break;
				}
			}	
		}


		for(Elem=0;Elem<NoOfElements;Elem++)
		{
			if(CC[Elem].Type == OUTLET)
			{
				fscanf(Restart,"%lf\t",&CDO[Elem].FCVel[0][0]);
				fscanf(Restart,"%lf\t",&CD[Elem].FCVel[0][0]);
				fscanf(Restart,"%lf\t",&CDO[Elem].FCVel[0][1]);
				fprintf(Restart,"%.9lf\n",CD[Elem].FCVel[0][1]);
			}
		}
		UpdatePressure(NoOfElements, CC, CDO, CD);
		UpdateVelocity(NoOfElements,CC,CDO);
		UpdateVelocity(NoOfElements,CC,CD);
		fclose(Restart);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		printf("\nDatas Succesfully Taken and Restarted\n");
	}
	else
	{
		Iter = 0;
		/*switch(COMP)
		{
		case 1:
			system("del ScreenOut.txt");
			break;
		case 2:
			system("rm -f ScreenOut.txt");
			break;
		default:
			printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
			exit(0);
		}*/
		
		InitialConditions(NoOfElements, NoOfNodes, CD, CDO, Node);
	}
	
	do
	{
		Iter++;
		printf("\nIter: %d ", Iter);
		Solver(NoOfElements, NoOfNodes, NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Re, DelT,Iter);
		//ApplyBC(NoOfElements, NoOfNodes,  NoOfBodies, CC, CD, CDO, Node, BE, BodyElems, Iter);
		PostSolver(NoOfElements, NoOfNodes, CC, CD, CDO, Node, Re, DelT, Iter, COMP);
	}while(Iter < MAXITER);
	
	printf("\nJai Sai Krishna\n");

	free(Node);
	free(CC);
	free(CD);
	free(CDO);
	free(BE);
	free(BodyElems);
	time (&end);
	dif = difftime (end,start);
	printf ("It took you %e seconds to run your code.\n", dif );
	return 0;
}


void GetNoOfElemAndNodes(int *NElem, int *NNode)
{
	FILE *FPElem;
	char FileElemName[100];

	int E_NUMNP, E_NELEM, E_NGRPS, E_NBSETS, E_NDFCD, E_NDFVL;
	char E_HEADNUT[100];

	FILE *FPNode;
	char FileNodeName[100];
	int N_NUMNP, N_NELEM, N_NGRPS, N_NBSETS, N_NDFCD, N_NDFVL;
	char N_HEADNUT[100];

	int t;

	char *t1;

	switch (COMP)
	{
	case 1:
		sprintf(FileElemName,"Mesh\\Elements.neu");
		sprintf(FileNodeName,"Mesh\\Nodes.neu");
		break;
	case 2:
		sprintf(FileElemName,".//Mesh//Elements.neu");
		sprintf(FileNodeName,".//Mesh//Nodes.neu");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}
	
	if( (FPElem = fopen(FileElemName,"rt") )==NULL)
	{
		printf("Unable to Open Mesh//Elements.neu\n");
		printf("Program Terminated.");
		//getch();
		exit(3);
	}

	if( (FPNode = fopen(FileNodeName,"rt") )==NULL)
	{
		printf("Unable to Open Mesh//Nodes.neu");
		printf("Program Terminated.");
		//getch();
		exit(3);
	}

	//printf("Files %s and %s Opened Succesfully\n", ElementFileName, NodeFileName);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("Control Info Not Match between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("File Header Not Matching between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}

	/*if(strcmp(E_HEADNUT,"** GAMBIT NEUTRAL FILE")!=0)
	{
		printf("The Input File %s is not a GAMBIT NEUTRAL FILE.\nProgram Terminated.", ElementFileName);
		exit(0);
	}*/

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("File ID Not Matching between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("File Version Not Matching between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("File Creation Date Not Matching between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	if(strcmp(E_HEADNUT,N_HEADNUT)!=0)
	{
		printf("File Parameters Not Matching between 2 files\n");
		printf("Program Terminated.");
		exit(1);
	}

	t = fscanf(FPElem, "%d %d %d %d %d %d", &E_NUMNP, &E_NELEM, &E_NGRPS, &E_NBSETS, &E_NDFCD, &E_NDFVL);
	t = fscanf(FPNode, "%d %d %d %d %d %d", &N_NUMNP, &N_NELEM, &N_NGRPS, &N_NBSETS, &N_NDFCD, &N_NDFVL);
	if( (E_NUMNP != N_NUMNP)|| (E_NELEM != N_NELEM) || (E_NGRPS != N_NGRPS) || (E_NBSETS != N_NBSETS) || (E_NDFCD != N_NDFCD) || (E_NDFVL != N_NDFVL) )
	{
		printf("Essential Parameters Missmatch\n");
		printf("Program Terminated");
	}
	*NElem = E_NELEM;
	*NNode = E_NUMNP;
}



void InitialConditions(int NoOfElements, int NoOfNodes, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node)
{
	int i;
	printf("\nDefinig Initial Conditions...");
	for(i=0;i<NoOfElements;i++)
	{
		CDO[i].CCP = Pinit;
		CDO[i].CCT = Tinit;

		CDO[i].FCP[0] = Pinit;
		CDO[i].FCP[1] = Pinit;
		CDO[i].FCP[2] = Pinit;

		CDO[i].CCVel[0] = Uinit;
		CDO[i].CCVel[1] = Vinit;//change to 0.0

		CD[i].FCT[0] = CDO[i].FCT[0] = Tinit;
		CD[i].FCT[0] = CDO[i].FCT[1] = Tinit;
		CD[i].FCT[0] = CDO[i].FCT[2] = Tinit;


		CD[i].CCP = Pinit;
		CD[i].CCT = Tinit;

		CD[i].FCP[0] = Pinit;
		CD[i].FCP[1] = Pinit;
		CD[i].FCP[2] = Pinit;

		CD[i].CCVel[0] = Uinit;
		CD[i].CCVel[1] = Vinit;
	}



	printf("\nInitial Conditions Defined.");
}


int RestartCheck(int NoOfElements, int NoOfNodes, int Re)
{
	FILE *FP;
	int Iter;
	int ConElements, ConNodes, ConRe;
	char FileRestart[100];
	
	switch (COMP)
	{
	case 1:
		sprintf(FileRestart,"Restart\\Config.cnf");
		break;
	case 2:
		sprintf(FileRestart,".//Restart//Config.cnf");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	printf("\n\n\nLooking For Restart Data...\n\n\n");

	if( (FP = fopen(FileRestart,"r+")) == NULL)
	{
		printf("\nUnable to open Restart config File\n");
		return 0;
	}
	fscanf(FP, "%d%d%d%d", &ConElements, &ConNodes, &ConRe, &Iter);
	//printf("%d\n",Iter);
	//getch();
	fclose(FP);
	if(ConElements == NoOfElements && ConNodes == NoOfNodes && ConRe == Re)
		return Iter;
	else
		return 0;
}
