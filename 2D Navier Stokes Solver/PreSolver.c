/*
//PreSolver.c
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/

#include"PreSolver.h"
#include<math.h>
//#include<process.h>
//#include<conio.h>
#include<ctype.h>
#include<string.h>

void PreSolver(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, NODE *Node, int Comp)
{
	int i;
	FILE *Fwrite;
	int nElem, nNode;

	char FileElemName[100], FileNodeName[100], FileNeighName[100];
	switch (Comp)
	{
	case 1:
		sprintf(FileElemName,"Mesh\\Elements.neu");
		sprintf(FileNodeName,"Mesh\\Nodes.neu");
		sprintf(FileNeighName, "Mesh\\Neigh.dat");
		break;

	case 2:
		sprintf(FileElemName,".//Mesh//Elements.neu");
		sprintf(FileNodeName,".//Mesh//Nodes.neu");
		sprintf(FileNeighName, ".//Mesh//Neigh.dat");
		break;

	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	ReadNeutralFiles(FileElemName, FileNodeName, NoOfElements, NoOfNodes, CC, Node);
	//WritePLTFILE(NoOfElements, NoOfNodes, CC, Node);
	//exit(0);
	OrientCounterClockwise(NoOfElements, CC, Node);
	
	
	printf("\nFinding Neighbours...");
	
	
	if( (Fwrite = fopen(FileNeighName,"r")) == NULL)
	{
		printf("\nNeigh.dat File not Found...\n");
		for(i=0;i<NoOfElements;i++)
		{
			FindNeighbour(i, NoOfElements, CC);
		}
	
		if( (Fwrite = fopen("Neigh.dat","w+")) == NULL)
		{
			printf("Unable to open Out File");
			exit(0) ;
		}
		fprintf(Fwrite,"%d\t%d\n",NoOfElements, NoOfNodes);
		for(i=0;i<NoOfElements;i++)
		{
			fprintf(Fwrite,"%d\t%d\t%d\n",CC[i].Neighbour[0], CC[i].Neighbour[1], CC[i].Neighbour[2]);
		}
	}
	else
	{
		fscanf(Fwrite,"%d\t%d\n",&nElem, &nNode);
		if(nElem == NoOfElements && nNode == NoOfNodes)
		{
			for(i=0;i<NoOfElements;i++)
			{
				fscanf(Fwrite,"%d\t%d\t%d\n",&CC[i].Neighbour[0],&CC[i].Neighbour[1],&CC[i].Neighbour[2]);
			}
		}
		else
		{
			fclose(Fwrite);
			if(Comp == 1)
				system ("del Mesh\\Neigh.dat");
			if(Comp == 2)
				system ("rm -f .//Mesh//Neigh.dat");
			for(i=0;i<NoOfElements;i++)
			{
				FindNeighbour(i, NoOfElements, CC);
			}
		
			if( (Fwrite = fopen(FileNeighName,"w+")) == NULL)
			{
				printf("Unable to open Out File");
				exit(0) ;
			}
			fprintf(Fwrite,"%d\t%d\n",NoOfElements, NoOfNodes);
			for(i=0;i<NoOfElements;i++)
			{
				fprintf(Fwrite,"%d\t%d\t%d\n",CC[i].Neighbour[0], CC[i].Neighbour[1], CC[i].Neighbour[2]);
			}
		}
	}
	fclose(Fwrite);
	
	printf("\nNeighbours for all Elements Found."); 
	FindArea(NoOfElements, CC, Node);
	WritePLTFILE(NoOfElements, NoOfNodes, CC, Node);
	ComputeNeighDetails(NoOfElements, CC, Node);
	//getch();
	FindKs(NoOfElements, CC);
	printf("Pre solver over");
	
	
	ComputeNodeDistance(NoOfNodes, CC, Node);
	//getch();
	//exit(0);
}

void ReadNeutralFiles(char *ElementFileName, char *NodeFileName, int NoOfElements, int NoOfNodes, CELLDETAILS *CD, NODE *Node)
{
	FILE *FPElem;
	unsigned int E_NUMNP, E_NELEM, E_NGRPS, E_NBSETS, E_NDFCD, E_NDFVL;
	char E_HEADNUT[100];

	FILE *FPNode;
	unsigned int N_NUMNP, N_NELEM, N_NGRPS, N_NBSETS, N_NDFCD, N_NDFVL;
	char N_HEADNUT[100];

	unsigned int Elem, Id, t;
	int a, b;
	char *t1;

	unsigned int i;
	unsigned j;
	unsigned int ELEM, ELEMTYPE, FACE;
	char Code[CODE][7]={"INLET", "OUTLET", "TOP", "BOTTOM", "BODY1"};

	int NODE;

	if( (FPElem = fopen(ElementFileName,"rt") )==NULL)
	{
		printf("Unable to Open %s",ElementFileName);
		printf("Program Terminated.");
		exit(3);
	}

	if( (FPNode = fopen(NodeFileName,"rt") )==NULL)
	{
		printf("Unable to Open %s",NodeFileName);
		printf("Program Terminated.");
		exit(3);
	}

	printf("Files %s and \n%s Opened Succesfully\n", ElementFileName, NodeFileName);

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

	printf("Preliminary tests on the files were succesful.\nExtracting Grid From the files...\n");
	printf("Total number of Elements is %d\nTotal number of Nodes is %d\n",E_NELEM, E_NUMNP);
	//fgets(FP,)
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);


	for(Elem=0;Elem<E_NUMNP;Elem++)
	{
		t = fscanf(FPElem, "%d %lf %lf", &Id, &Node[Elem].x, &Node[Elem].y );
		t = fscanf(FPNode, "%d %lf %lf", &Id, &Node[Elem].x, &Node[Elem].y );
		//printf("%d\t%lf\t%lf\n", Id, Node[Elem].x, Node[Elem].y );
		Node[Elem].Type = 0;
	}

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);

	for(Elem=0;Elem<E_NELEM;Elem++)
	{
		t = fscanf(FPElem, "%d %d %d %d %d %d", &Id, &a, &b, &CD[Elem].Connect[0], &CD[Elem].Connect[1], &CD[Elem].Connect[2] );
		t = fscanf(FPNode, "%d %d %d %d %d %d", &Id, &a, &b, &CD[Elem].Connect[0], &CD[Elem].Connect[1], &CD[Elem].Connect[2]);
		
		Node[CD[Elem].Connect[0]-1].ShareElements[Node[CD[Elem].Connect[0]-1].ShareCount] = Elem;
		Node[CD[Elem].Connect[0]-1].ShareCount ++;
		Node[CD[Elem].Connect[1]-1].ShareElements[Node[CD[Elem].Connect[1]-1].ShareCount] = Elem;
		Node[CD[Elem].Connect[1]-1].ShareCount ++;
		Node[CD[Elem].Connect[2]-1].ShareElements[Node[CD[Elem].Connect[2]-1].ShareCount] = Elem;
		Node[CD[Elem].Connect[2]-1].ShareCount ++;

		CD[Elem].Neighbour[0] = CD[Elem].Neighbour[2] = CD[Elem].Neighbour[1] = -2;
		CD[Elem].Type = 0;
		//printf("%d\t %d\t %d\t %d\t %d\t %d\n", Id, a, b, CD[Elem].Connect[0], CD[Elem].Connect[1], CD[Elem].Connect[2]);
	}
//printf("Over");
	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	t1 = fgets(E_HEADNUT,100, FPElem);
	t1 = fgets(N_HEADNUT,100, FPNode);
	//puts(N_HEADNUT);

	for(Elem=0;Elem<E_NELEM;Elem++)
	{
		t = fscanf(FPElem, "%d", &Id);
		t = fscanf(FPNode, "%d", &Id);
		//printf("%d\n", Id);
	}

///printf("E_NBSETS:%d",E_NBSETS);
	for(i=0;i<E_NBSETS;i++)
	{
		unsigned int ITYPE, NENTRY, NVALUES, IBCODE1;
		char NAME[20];

		t1 = fgets(E_HEADNUT,100, FPElem);
		t1 = fgets(E_HEADNUT,100, FPElem);
		t1 = fgets(E_HEADNUT,100, FPElem);
		t1 = fgets(E_HEADNUT,100, FPElem);

		//puts(E_HEADNUT);
		sscanf(E_HEADNUT,"%s %d %d %d %d",NAME, &ITYPE, &NENTRY, &NVALUES, &IBCODE1);
		//printf("\n%s %d %d %d %d ",NAME, ITYPE, NENTRY, NVALUES, IBCODE1);
		for(j=0;j<strlen(NAME);j++)
		{
			NAME[j]= toupper(NAME[j]);
		}
		for(j=0;j<CODE;j++)
		{
			if(strcmp(NAME,Code[j])==0)
			{
				IBCODE1 = j+1;
				break;
			}
		}
		//printf("IBCODE :%d\n", IBCODE1);
		
		for(j=0;j<NENTRY;j++)
		{
			t = fscanf(FPElem,"%d %d %d", &ELEM, &ELEMTYPE, &FACE);
			CD[ELEM-1].Type = IBCODE1;
			//printf("%d\n",ELEM);
		}
		//getch();
		//scanf("%d",&t);
	}

	for(i=0;i<N_NBSETS;i++)
	{
		unsigned int ITYPE, NENTRY, NVALUES, IBCODE1;
		char NAME[20];

		t1 = fgets(N_HEADNUT,100, FPNode);
		//puts(N_HEADNUT);
		t1 = fgets(N_HEADNUT,100, FPNode);
		//puts(N_HEADNUT);
		t1 = fgets(N_HEADNUT,100, FPNode);
		//puts(N_HEADNUT);
		t1 = fgets(N_HEADNUT,100, FPNode);
		//puts(N_HEADNUT);
		//scanf("%d",&t);
		sscanf(N_HEADNUT,"%s %d %d %d %d",NAME, &ITYPE, &NENTRY, &NVALUES, &IBCODE1);
		//printf("%s %d %d %d %d ",NAME, ITYPE, NENTRY, NVALUES, IBCODE1);
		//scanf("%d",&t);
		for(j=0;j<strlen(NAME);j++)
		{
			NAME[j]=toupper(NAME[j]);
		}
		for(j=0;j<CODE;j++)
		{
			if(strcmp(NAME,Code[j])==0)
			{
				IBCODE1 = j+1;
				break;
			}
		}
		//printf("IBCODE :%d\n", IBCODE1);
		//scanf("%d",&t);
		
		for(j=0;j<NENTRY;j++)
		{
			t = fscanf(FPNode,"%d", &NODE);
			Node[NODE-1].Type = IBCODE1;
			//printf("%d\n",NODE);
			//getch();
		}
		//scanf("%d",&t);
		//printf("Hello");
	}
	fclose(FPElem);
	fclose(FPNode);
	if((NoOfNodes != E_NUMNP) ||(NoOfElements != E_NELEM))
	{
		printf("Miss match in the end of presolver File reading.\nProgram Terminated");
		exit(0);
	}
}

void OrientCounterClockwise(int NoOfElements, CELLDETAILS *CC, NODE *Node)
{
	int i, j;
	int temp;
	printf("\nOrienting To Clockwise Direction...\n");
	for(i=0;i<NoOfElements;i++)
	{
		double ABx, ABy, ACx, ACy, BCx, BCy;
		
		if(i==0 || i==318||i==837)
			i=i;
		CC[i].Xp = (Node[CC[i].Connect[0]-1].x + Node[CC[i].Connect[1]-1].x + Node[CC[i].Connect[2]-1].x)/3.0;
		CC[i].Yp = (Node[CC[i].Connect[0]-1].y + Node[CC[i].Connect[1]-1].y + Node[CC[i].Connect[2]-1].y)/3.0;

		ABx = Node[CC[i].Connect[1]-1].x - Node[CC[i].Connect[0]-1].x;
		ABy = Node[CC[i].Connect[1]-1].y - Node[CC[i].Connect[0]-1].y;
		ACx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[0]-1].x;
		ACy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[0]-1].y;
		BCx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[1]-1].x;
		BCy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[1]-1].y;
		
		if( (ABx*BCy - ABy*BCx) < 0)
		{
			if((ABx*ACy - ABy*ACx)<0)
			{
				temp = CC[i].Connect[0];
				CC[i].Connect[0] = CC[i].Connect[1];
				CC[i].Connect[1] = temp;
				i--;
				continue;
			}
			else
			{
				temp = CC[i].Connect[2];
				CC[i].Connect[2] = CC[i].Connect[1];
				CC[i].Connect[1] = temp;
				i--;
				continue;
			}
		}
		else
		{
			if((ABx*ACy - ABy*ACx)<0)
			{
				temp = CC[i].Connect[0];
				CC[i].Connect[0] = CC[i].Connect[2];
				CC[i].Connect[2] = temp;
				i--;
				continue;
			}
		}
		CC[i].Xab = ABx = Node[CC[i].Connect[1]-1].x - Node[CC[i].Connect[0]-1].x;
		CC[i].Yab = ABy = Node[CC[i].Connect[1]-1].y - Node[CC[i].Connect[0]-1].y;
		CC[i].Xca = ACx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[0]-1].x;
		CC[i].Xca = -CC[i].Xca;
		CC[i].Yca = ACy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[0]-1].y;
		CC[i].Yca = -CC[i].Yca;
		CC[i].Xbc = BCx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[1]-1].x;
		CC[i].Ybc = BCy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[1]-1].y;
		if( (ABx*BCy - ABy*BCx) < 0)
		{
			printf("Checking Falied at 1 for Element %d",i);
			exit(1);
		}
		if((ABx*ACy - ABy*ACx)<0)
		{
			printf("Checking Falied at 2 for Element %d",i);
			exit(1);
		}
			
		if((ACx*BCy - ACy*BCx)<0)
		{
			printf("Checking Falied at 3 for Element %d",i);
			exit(1);
		}

		if(CC[i].Type == 1 || CC[i].Type == 2 || CC[i].Type == 3 || CC[i].Type == 4 || CC[i].Type == 5)
		{
			for(j=0;j<2;j++)
				if (Node[CC[i].Connect[2]-1].Type != 0) //Third node not on any of the boundary
				{
					temp = CC[i].Connect[0];
					CC[i].Connect[0] = CC[i].Connect[1];
					CC[i].Connect[1] = CC[i].Connect[2];
					CC[i].Connect[2] = temp;
				}

			CC[i].Xab = ABx = Node[CC[i].Connect[1]-1].x - Node[CC[i].Connect[0]-1].x;
			CC[i].Yab = ABy = Node[CC[i].Connect[1]-1].y - Node[CC[i].Connect[0]-1].y;
			CC[i].Xca = ACx = Node[CC[i].Connect[0]-1].x - Node[CC[i].Connect[2]-1].x; //CAx
			CC[i].Yca = ACy = Node[CC[i].Connect[0]-1].y - Node[CC[i].Connect[2]-1].y;//CAy
			CC[i].Xbc = BCx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[1]-1].x;
			CC[i].Ybc = BCy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[1]-1].y;
		}
		if(CC[i].Type==0)
		{
			
			/*if( Node[CC[i].Connect[0]-1].Type != 0)
			{
				CC[i].Type = 10 + Node[CC[i].Connect[0]-1].Type;
				
			}
			else if( Node[CC[i].Connect[1]-1].Type != 0)
				{
					CC[i].Type = 10 + Node[CC[i].Connect[1]-1].Type;
					
				}
			else if( Node[CC[i].Connect[2]-1].Type != 0)
				{
					CC[i].Type = 10 + Node[CC[i].Connect[2]-1].Type;
					
				}
			*/
			int j;
			for(j=0; j<3; j++)
			{
				if(Node[CC[i].Connect[j]-1].Type != 0)
					CC[i].Type = 10 + Node[CC[i].Connect[j]-1].Type;
			}
		}
		if(CC[i].Type == 11 || CC[i].Type == 12 || CC[i].Type == 13 || CC[i].Type == 14 || CC[i].Type == 15)
		{
			int j, temp;
			double ABx, ABy, ACx, ACy, BCx, BCy;
			for(j=0;j<2;j++)
				if (Node[CC[i].Connect[0]-1].Type == 0) //First Node on any of the boundary
				{
					temp = CC[i].Connect[0];
					CC[i].Connect[0] = CC[i].Connect[1];
					CC[i].Connect[1] = CC[i].Connect[2];
					CC[i].Connect[2] = temp;
				}

			CC[i].Xab = ABx = Node[CC[i].Connect[1]-1].x - Node[CC[i].Connect[0]-1].x;
			CC[i].Yab = ABy = Node[CC[i].Connect[1]-1].y - Node[CC[i].Connect[0]-1].y;
			CC[i].Xca = ACx = Node[CC[i].Connect[0]-1].x - Node[CC[i].Connect[2]-1].x; //CAx
			CC[i].Yca = ACy = Node[CC[i].Connect[0]-1].y - Node[CC[i].Connect[2]-1].y;//CAy
			CC[i].Xbc = BCx = Node[CC[i].Connect[2]-1].x - Node[CC[i].Connect[1]-1].x;
			CC[i].Ybc = BCy = Node[CC[i].Connect[2]-1].y - Node[CC[i].Connect[1]-1].y;
		}
	}
	printf("Orientation Completed.");
}

void FindArea(int NoOfElements, CELLDETAILS *CC, NODE *Node)
{
	int i;
	double S, AB, BC, CA;
	double Total = 0.0;
	printf("\nCalculating Area...");
	for(i=0;i<NoOfElements;i++)
	{
		AB = pow( (pow(CC[i].Xab,2)+pow(CC[i].Yab,2)), 0.5);
		BC = pow( (pow(CC[i].Xbc,2)+pow(CC[i].Ybc,2)), 0.5);
		CA = pow( (pow(CC[i].Xca,2)+pow(CC[i].Yca,2)), 0.5);
		S = (AB + BC + CA)/2;
		CC[i].Area = pow( (S*(S - AB)*(S - BC)*(S - CA)),0.5);
		Total += CC[i].Area;	
	}
	printf("\nArea Computed."); 
}


void FindNeighbour(int Elem, int NoOfElements, CELLDETAILS *CC)
{
	int AssignCount=0;
	int Match=0;
	int i,j,k;
	for(i=0; i<NoOfElements; i++)
	{
		if(i!= Elem)
		{
			Match=0;
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					if(CC[i].Connect[j] == CC[Elem].Connect[k])
					{
						Match++;
						if(Match == 2)
						{
							CC[Elem].Neighbour[AssignCount] = i;
							AssignCount++;
							j=k=10;
								
						}
					}
				}
			}
		}
		if(AssignCount == 3)
			return;
	}	
}

void ComputeNeighDetails(int NoOfElements, CELLDETAILS *CC, NODE *Node)
{
	int i;
	int CellA, CellB, CellC, CellD, CellL, CellJ, CellI, CellF, CellG;
	int NodeA, NodeB, NodeC, NodeD, NodeG, NodeJ;
	double Mab, Mcg; // Slope of AB & CG Where G is the reflection of C on AB
	double Cab, Ccg; //y=mx+c <- of ab, cg
	double Xp, Yp; // Co-ordinate of intersection of AB & CG (P)
	double Xg, Yg; // Co-ordinate of G
	//FILE *FPN, *FPC;
	int count;
	printf("\nComputing Neighbour Details...");
	i=0;

	
	//FPN = fopen("ghostN.csv","w+");
	//FPC = fopen("ghostC.csv","w+");
	
	count=1;

	for(i=0; i<NoOfElements; i++)
	{
		
		switch(CC[i].Type)
		{
		
		case 11:
		case 12:
		case 13:
		case 14:
		//	break;
		case 15:
		case 0:
		
			NodeA = CC[i].Connect[0];
			NodeB = CC[i].Connect[1];
			NodeC = CC[i].Connect[2];
			
			CellA = CommonNodes(i, NodeB, NodeC, CC);
			CellB = CommonNodes(i, NodeC, NodeA, CC);
			CellC = CommonNodes(i, NodeA, NodeB, CC);

			CC[i].Neighbour[0] = CellA;
			CC[i].Neighbour[1] = CellB;
			CC[i].Neighbour[2] = CellC;

			CC[i].CellD = CellD = FindCell(i, CellA, NodeB, CC);
			
			CC[i].CellL = CellL = FindCell(i, CellA, NodeC, CC);
			
			CC[i].CellJ = CellJ = FindCell(i, CellB, NodeC, CC);
			CC[i].CellI = CellI = FindCell(i, CellB, NodeA, CC);
			
			CC[i].CellG = CellG = FindCell(i, CellC, NodeA, CC);
			CC[i].CellF = CellF = FindCell(i, CellC, NodeB, CC);

			NodeD = FindMissingNode(CellA, NodeB, NodeC, CC);
			NodeJ = FindMissingNode(CellB, NodeA, NodeC, CC);
			NodeG = FindMissingNode(CellC, NodeA, NodeB, CC);


			CC[i].Index4 = FindIndex(CellA, NodeD, CC);
			CC[i].Index5 = FindIndex(CellA, NodeB, CC);
			CC[i].Index11 = FindIndex(CellC, NodeG, CC);
			CC[i].Index12 = FindIndex(CellC, NodeA, CC);
			CC[i].Index18 = FindIndex(CellB, NodeJ, CC);
			CC[i].Index19 = FindIndex(CellB, NodeC, CC);

			CC[i].Apa = CC[i].Area + CC[CellA].Area;
			CC[i].Apb = CC[i].Area + CC[CellB].Area;
			CC[i].Apc = CC[i].Area + CC[CellC].Area;

			CC[i].Xab = Node[NodeB-1].x - Node[NodeA-1].x;
			CC[i].Xbc = Node[NodeC-1].x - Node[NodeB-1].x;
			CC[i].Xca = Node[NodeA-1].x - Node[NodeC-1].x;

			CC[i].Xag = Node[NodeG-1].x - Node[NodeA-1].x;
			CC[i].Xbd = Node[NodeD-1].x - Node[NodeB-1].x;
			CC[i].Xcj = Node[NodeJ-1].x - Node[NodeC-1].x;
			CC[i].Xdc = Node[NodeC-1].x - Node[NodeD-1].x;
			CC[i].Xgb = Node[NodeB-1].x - Node[NodeG-1].x;
			CC[i].Xja = Node[NodeA-1].x - Node[NodeJ-1].x;

			CC[i].Yab = Node[NodeB-1].y - Node[NodeA-1].y;
			CC[i].Ybc = Node[NodeC-1].y - Node[NodeB-1].y;
			CC[i].Yca = Node[NodeA-1].y - Node[NodeC-1].y;

			CC[i].Yag = Node[NodeG-1].y - Node[NodeA-1].y;
			CC[i].Ybd = Node[NodeD-1].y - Node[NodeB-1].y;
			CC[i].Ycj = Node[NodeJ-1].y - Node[NodeC-1].y;
			CC[i].Ydc = Node[NodeC-1].y - Node[NodeD-1].y;
			CC[i].Ygb = Node[NodeB-1].y - Node[NodeG-1].y;
			CC[i].Yja = Node[NodeA-1].y - Node[NodeJ-1].y;
			
			CC[i].XAP = CC[i].Xp - CC[CellA].Xp;
			CC[i].YAP = CC[i].Yp - CC[CellA].Yp;

			CC[i].XBP = CC[i].Xp - CC[CellB].Xp;
			CC[i].YBP = CC[i].Yp - CC[CellB].Yp;

			CC[i].XCP = CC[i].Xp - CC[CellC].Xp;
			CC[i].YCP = CC[i].Yp - CC[CellC].Yp;
			break;

		case 1://Inlet
		case 2://Outlet
		case 3://Top
		case 4://Bottom
			//break;
		case 5://body
			//printf("Body  ");
			
			NodeA = CC[i].Connect[0];
			NodeB = CC[i].Connect[1];
			NodeC = CC[i].Connect[2];
			
			if(fabs(Node[NodeB-1].y - Node[NodeA-1].y)<1e-8)
			{
				Xg = Node[NodeC-1].x;
				Yg = 2 * Node[NodeA-1].y - Node[NodeC-1].y;

				CC[i].Xag = Xg - Node[NodeA-1].x;
				CC[i].Xgb = Node[NodeB-1].x - Xg;
				CC[i].Yag = Yg - Node[NodeA-1].y;
				CC[i].Ygb = Node[NodeB-1].y - Yg;

			}
			else if(fabs(Node[NodeB-1].x - Node[NodeA-1].x)<1e-8)
			{
				Xg = 2 * Node[NodeA-1].x - Node[NodeC-1].x;
				Yg = Node[NodeC-1].y;
				

				CC[i].Xag = Xg - Node[NodeA-1].x;
				CC[i].Xgb = Node[NodeB-1].x - Xg;
				CC[i].Yag = Yg - Node[NodeA-1].y;
				CC[i].Ygb = Node[NodeB-1].y - Yg;
			}
			/*else if((Node[NodeB-1].x - Node[NodeA-1].x)<1e-5)
			{
				CC[i].Xag = -CC[i].Xca ;
				CC[i].Xgb = -CC[i].Xbc ;
				CC[i].Yag = -CC[i].Yca ;
				CC[i].Ygb = -CC[i].Ybc ;
			}*/
			else
			{
				Mab = (Node[NodeB-1].y - Node[NodeA-1].y) / (Node[NodeB-1].x - Node[NodeA-1].x);
			
				Mcg = -1/Mab;

				Cab = Node[NodeA-1].y - Mab * Node[NodeA-1].x;
				Ccg = Node[NodeC-1].y - Mcg * Node[NodeC-1].x;

				Xp = (Ccg - Cab) / (Mab - Mcg);
				Yp = Xp*Mab + Cab;
				Xg = 2*Xp-Node[NodeC-1].x;
				Yg = 2*Yp-Node[NodeC-1].y;

				CC[i].Xag = Xg - Node[NodeA-1].x;
				CC[i].Xgb = Node[NodeB-1].x - Xg;
				CC[i].Yag = Yg - Node[NodeA-1].y;
				CC[i].Ygb = Node[NodeB-1].y - Yg;
			}

			
			//fprintf(FPN,"%lf\t%lf\n",Xg, Yg);
			//fprintf(FPC,"%d\t%d\t%d\n",CC[i].Connect[0],CC[i].Connect[1],count);
			//count++;

			CellA = CommonNodes(i, NodeB, NodeC, CC);
			CellB = CommonNodes(i, NodeC, NodeA, CC);
			CellC = UNAVAILABLEDATA;

			CC[i].Neighbour[0] = CellA;
			CC[i].Neighbour[1] = CellB;
			CC[i].Neighbour[2] = CellC;

			CC[i].CellD = CellD = FindCell(i, CellA, NodeB, CC);
			CC[i].CellL = CellL = FindCell(i, CellA, NodeC, CC);
			if(CC[i].CellD == UNAVAILABLEDATA)
			{
				CC[i].CellD = CellA;
			}


			CC[i].CellJ = CellJ = FindCell(i, CellB, NodeC, CC);
			CC[i].CellI = CellI = FindCell(i, CellB, NodeA, CC);
			if(CC[i].CellI == UNAVAILABLEDATA)
			{
				CC[i].CellI = CellB;
			}
			CC[i].CellG = UNAVAILABLEDATA;
			CC[i].CellF = UNAVAILABLEDATA;

			NodeD = FindMissingNode(CellA, NodeB, NodeC, CC);
			NodeG = UNAVAILABLEDATA;
			NodeJ = FindMissingNode(CellB, NodeA, NodeC, CC);
			
			CC[i].Index4 = FindIndex(CellA, NodeD, CC);
			CC[i].Index5 = FindIndex(CellA, NodeB, CC);
			CC[i].Index11 = UNAVAILABLEDATA;
			CC[i].Index12 = UNAVAILABLEDATA;
			CC[i].Index18 = FindIndex(CellB, NodeJ, CC);
			CC[i].Index19 = FindIndex(CellB, NodeC, CC);

			CC[i].Apa = CC[i].Area + CC[CellA].Area;
			CC[i].Apb = CC[i].Area + CC[CellB].Area;
			CC[i].Apc = CC[i].Area * 2.0;


			CC[i].Xab = Node[NodeB-1].x - Node[NodeA-1].x;
			CC[i].Xbc = Node[NodeC-1].x - Node[NodeB-1].x;
			CC[i].Xca = Node[NodeA-1].x - Node[NodeC-1].x;

			CC[i].Xbd = Node[NodeD-1].x - Node[NodeB-1].x;
			CC[i].Xcj = Node[NodeJ-1].x - Node[NodeC-1].x;
			CC[i].Xdc = Node[NodeC-1].x - Node[NodeD-1].x;
			CC[i].Xja = Node[NodeA-1].x - Node[NodeJ-1].x;

			CC[i].Yab = Node[NodeB-1].y - Node[NodeA-1].y;
			CC[i].Ybc = Node[NodeC-1].y - Node[NodeB-1].y;
			CC[i].Yca = Node[NodeA-1].y - Node[NodeC-1].y;

			CC[i].Ybd = Node[NodeD-1].y - Node[NodeB-1].y;
			CC[i].Ycj = Node[NodeJ-1].y - Node[NodeC-1].y;
			CC[i].Ydc = Node[NodeC-1].y - Node[NodeD-1].y;
			CC[i].Yja = Node[NodeA-1].y - Node[NodeJ-1].y;

			CC[i].XAP = CC[i].Xp - CC[CellA].Xp;
			CC[i].YAP = CC[i].Yp - CC[CellA].Yp;

			CC[i].XBP = CC[i].Xp - CC[CellB].Xp;
			CC[i].YBP = CC[i].Yp - CC[CellB].Yp;

			CC[i].XCP = (CC[i].Xp - ((Node[NodeA-1].x + Node[NodeB-1].x + Xg)/3.0));
			CC[i].YCP = (CC[i].Yp - ((Node[NodeA-1].y + Node[NodeB-1].y + Yg)/3.0));

			break;
		default://
			printf("\nExtra Boundary Entering...\nProgram Terminating...CC[%d].Type=%d\n",i,CC[i].Type);



			//exit(10);
		}
		/*count = 28340;
		if(i== count-1)
		{
			printf("Elem %d\n",count);
			printf("NodeA %d\n",NodeA);
			printf("NodeB %d\n",NodeB);
			printf("NodeC %d\n",NodeC);
			printf("NodeD %d\n",NodeD);
			printf("NodeG %d\n",NodeG);
			printf("NodeJ %d\n",NodeJ);
			printf("CellA %d\n",CellA);
			printf("CellB %d\n",CellB);
			printf("CellC %d\n",CellC);
			printf("CellD %d\n",CellD);
			printf("CellF %d\n",CellF);
			printf("CellG %d\n",CellG);
			printf("CellI %d\n",CellI);
			printf("CellJ %d\n",CellJ);
			printf("CellL %d\n",CellL);
		}*/
	}
	printf("\nNeighbour to Neighbour Details Obtained.");
	/*FILE *FP;
	FP = fopen("Cells.csv","w+");
	for(i=0;i<NoOfElements;i++)
	{
		fprintf(FP,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",i+1,CC[i].Neighbour[0]+1,CC[i].Neighbour[1]+1,CC[i].Neighbour[2]+1, CC[i].CellD+1, CC[i].CellL+1, CC[i].CellF+1, CC[i].CellG+1, CC[i].CellI+1, CC[i].CellJ+1);
	}
	fclose(FP);
	printf("Cells Printed");
	exit(0);*/
	//getch();
}

int CommonNodes(int CellP, int NodeA, int NodeB, CELLDETAILS *CC)
{
	int i,j,k;
	for(i=0; i<3; i++)
	{
		for(j=0; j<3; j++)
		{
			if(CC[CellP].Neighbour[i] < 0)
			{
				printf("Error in Common Nodes trying to Access Negative Index");
				printf("CellP:%d\tNeigh[%d]:%d\tType:%d",CellP,i,CC[CellP].Neighbour[i],CC[CellP].Type);
				//getch();
				exit(0);
			}
			if(CC[CC[CellP].Neighbour[i]].Connect[j] == NodeA)
				for(k=0;k<3;k++)
				{
					if(CC[CC[CellP].Neighbour[i]].Connect[k] == NodeB)
						return CC[CellP].Neighbour[i];
				}
		}
	}
	printf("Error in CommonNodes\n");
	exit(9);
}

int FindCell(int CellP, int CellA, int NodeA, CELLDETAILS *CC)
{
	int i,j;
	int Neigh;
	for(i=0; i<3; i++)
	{
		Neigh = CC[CellA].Neighbour[i];
		if(Neigh < 0)
		{
			return UNAVAILABLEDATA;//Previous Statement is "return CellA;"
		}
		for(j=0; j<3; j++)
		{
			if(Neigh != CellP)
				if(CC[Neigh].Connect[j] == NodeA)
					return Neigh;
		}
	}
	printf("Error in FindCell\n");
	exit(9);
}

int FindMissingNode(int P, int NodeA, int NodeB, CELLDETAILS *CC)
{
	int Node, i;
	for(i=0;i<3;i++)
	{
		if(CC[P].Connect[i] != NodeA)
			if(CC[P].Connect[i] != NodeB)
			{
				Node = CC[P].Connect[i];
				i=10;
			}
	}
	return Node;
}

int FindIndex(int Cell, int Node, CELLDETAILS *CC)
{
	int i;
	for(i=0; i<3; i++)
	{
		if(Node == CC[Cell].Connect[i])
		{
			return i;
		}
	}
	printf("Error in Finding Index...\nProgram Terminated.\n");
	exit(9);
}

void FindKs(int NoOfElements, CELLDETAILS *CC)
{
	int Elem;
	double GK1, GK2, GK3, GK4, GK5, GK11, GK12, GK18, GK19;
	double Al, Ai, Aj, Ad, Af, Ag, Aa, Ab, Ac;
	double SumK;
	printf("\nFinding Ks...");

	
	for(Elem=0; Elem<NoOfElements; Elem++)
	{
		switch(CC[Elem].Type)
		{
		case DOM:
			Ac = CC[CC[Elem].Neighbour[2]].Area;
			Ad = CC[CC[Elem].CellD].Area;
			Af = CC[CC[Elem].CellF].Area;
			Ag = CC[CC[Elem].CellG].Area;
			Ai = CC[CC[Elem].CellI].Area;

			Aa = CC[CC[Elem].Neighbour[0]].Area;
			Ab = CC[CC[Elem].Neighbour[1]].Area;
			Al = CC[CC[Elem].CellL].Area;
			Aj = CC[CC[Elem].CellJ].Area;

			GK1 = ((CC[Elem].Yab * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Yab * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xab * CC[Elem].Xbc) / CC[Elem].Apa) + ((CC[Elem].Xab * CC[Elem].Xca) / CC[Elem].Apb);
			GK2 = ((CC[Elem].Ybc * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Ybc * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xbc * CC[Elem].Xab) / CC[Elem].Apc) + ((CC[Elem].Xbc * CC[Elem].Xca) / CC[Elem].Apb);
			GK3 = ((CC[Elem].Yca * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Yca * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xca * CC[Elem].Xab) / CC[Elem].Apc) + ((CC[Elem].Xca * CC[Elem].Xbc) / CC[Elem].Apa);
			GK4 = ((CC[Elem].Ydc * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xdc * CC[Elem].Xbc) / CC[Elem].Apa);
			GK5 = ((CC[Elem].Ybd * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xbd * CC[Elem].Xbc) / CC[Elem].Apa);
			GK11 = ((CC[Elem].Ygb * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Xgb * CC[Elem].Xab) / CC[Elem].Apc);
			GK12 = ((CC[Elem].Yag * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Xag * CC[Elem].Xab) / CC[Elem].Apc);
			GK18 = ((CC[Elem].Yja * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xja * CC[Elem].Xca) / CC[Elem].Apb);
			GK19 = ((CC[Elem].Ycj * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xcj * CC[Elem].Xca) / CC[Elem].Apb);

			break;

		case BODY1:
		case INLET:
		case OUTLET:

		case TOP:
		case BOTTOM:

			Aa = CC[CC[Elem].Neighbour[0]].Area;
			Ab = CC[CC[Elem].Neighbour[1]].Area;
			Al = CC[CC[Elem].CellL].Area;
			Aj = CC[CC[Elem].CellJ].Area;

			Ac = UNAVAILABLEDATA;
			if(CC[Elem].CellD == UNAVAILABLEDATA)
			{
				Ad = Aa;
			}
			else
			{
				Ad = CC[CC[Elem].CellD].Area;
			}
			Af = UNAVAILABLEDATA;
			Ag = UNAVAILABLEDATA;

			if(CC[Elem].CellI == UNAVAILABLEDATA)
			{
				Ai = Ab;
			}
			else
			{
				Ai = CC[CC[Elem].CellI].Area;
			}

			

			GK1 = ((CC[Elem].Yab * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Yab * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xab * CC[Elem].Xbc) / CC[Elem].Apa) + ((CC[Elem].Xab * CC[Elem].Xca) / CC[Elem].Apb);
			GK2 = ((CC[Elem].Ybc * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xbc * CC[Elem].Xca) / CC[Elem].Apb);//
			GK3 = ((CC[Elem].Yca * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xca * CC[Elem].Xbc) / CC[Elem].Apa);//
			GK4 = ((CC[Elem].Ydc * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xdc * CC[Elem].Xbc) / CC[Elem].Apa);
			GK5 = ((CC[Elem].Ybd * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xbd * CC[Elem].Xbc) / CC[Elem].Apa);
			GK11 = 0.0;//
			GK12 = 0.0;//
			GK18 = ((CC[Elem].Yja * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xja * CC[Elem].Xca) / CC[Elem].Apb);
			GK19 = ((CC[Elem].Ycj * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xcj * CC[Elem].Xca) / CC[Elem].Apb);

			break;
		
		
			break;

		case INNNUN:
		case BODYNNUN:
		case OUTNNUN:

		case TOPNNUN:
		case BOTNNUN:
			Aa = CC[CC[Elem].Neighbour[0]].Area;
			Ab = CC[CC[Elem].Neighbour[1]].Area;
			Al = CC[CC[Elem].CellL].Area;
			Aj = CC[CC[Elem].CellJ].Area;
			Ac = CC[CC[Elem].Neighbour[2]].Area;
			Ad = CC[CC[Elem].CellD].Area;
			Af = CC[CC[Elem].CellF].Area;
			
			if(CC[Elem].CellG == UNAVAILABLEDATA)
				Ag = UNAVAILABLEDATA;
			else
				Ag = CC[CC[Elem].CellG].Area;
			if(CC[Elem].CellI == UNAVAILABLEDATA)
				Ai = UNAVAILABLEDATA;
			else
				Ai = CC[CC[Elem].CellI].Area;
			if(CC[Elem].CellG != UNAVAILABLEDATA && CC[Elem].CellI != UNAVAILABLEDATA)
				CC[Elem].Type = DOM;

			GK1 = ((CC[Elem].Yab * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Yab * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xab * CC[Elem].Xbc) / CC[Elem].Apa) + ((CC[Elem].Xab * CC[Elem].Xca) / CC[Elem].Apb);
			GK2 = ((CC[Elem].Ybc * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Ybc * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xbc * CC[Elem].Xab) / CC[Elem].Apc) + ((CC[Elem].Xbc * CC[Elem].Xca) / CC[Elem].Apb);
			GK3 = ((CC[Elem].Yca * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Yca * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xca * CC[Elem].Xab) / CC[Elem].Apc) + ((CC[Elem].Xca * CC[Elem].Xbc) / CC[Elem].Apa);
			GK4 = ((CC[Elem].Ydc * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xdc * CC[Elem].Xbc) / CC[Elem].Apa);
			GK5 = ((CC[Elem].Ybd * CC[Elem].Ybc) / CC[Elem].Apa) + ((CC[Elem].Xbd * CC[Elem].Xbc) / CC[Elem].Apa);
			GK11 = ((CC[Elem].Ygb * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Xgb * CC[Elem].Xab) / CC[Elem].Apc);
			GK12 = ((CC[Elem].Yag * CC[Elem].Yab) / CC[Elem].Apc) + ((CC[Elem].Xag * CC[Elem].Xab) / CC[Elem].Apc);
			GK18 = ((CC[Elem].Yja * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xja * CC[Elem].Xca) / CC[Elem].Apb);
			GK19 = ((CC[Elem].Ycj * CC[Elem].Yca) / CC[Elem].Apb) + ((CC[Elem].Xcj * CC[Elem].Xca) / CC[Elem].Apb);
			break;

		
		default:
			printf("Default Case CC[%d].Type=%d",Elem, CC[Elem].Type);
			//getch();
			exit(0);
		}
		
		
		switch(CC[Elem].Type)
		{
		case DOM:
			CC[Elem].Kp = ((Ac * GK1) / CC[Elem].Apc) + ((Aa * GK2) / CC[Elem].Apa) + ((Ab * GK3) / CC[Elem].Apb);
			CC[Elem].Ka = ((CC[Elem].Area * GK2) / CC[Elem].Apa) + ((Al * GK4) / (Aa + Al)) + ((Ad * GK5) / (Aa+Ad));
			CC[Elem].Kb = (CC[Elem].Area * GK3 / CC[Elem].Apb) + (Ai *  GK18 / (Ab + Ai)) + (Aj * GK19 / (Ab + Aj));
			CC[Elem].Kc = (CC[Elem].Area * GK1 / CC[Elem].Apc) + (Af * GK11 / (Ac + Af)) + (Ag * GK12 / (Ag + Ac));
			CC[Elem].Kd = Aa * GK5 / (Aa + Ad);
			CC[Elem].Kf = Ac * GK11 / (Ac + Af);
			CC[Elem].Kg = Ac * GK12 / (Ac + Ag);
			CC[Elem].Ki = Ab * GK18 / (Ab + Ai);
			CC[Elem].Kj = Ab * GK19 / (Ab + Aj);
			CC[Elem].Kl = Aa * GK4 / (Aa + Al);
			break;
		
		case BODYNNUN:
		case INNNUN:
		case OUTNNUN:

		case TOPNNUN:
		case BOTNNUN:
			if(CC[Elem].CellG == UNAVAILABLEDATA && CC[Elem].CellI == UNAVAILABLEDATA)
			{
				CC[Elem].Kp = ((Ac * GK1) / CC[Elem].Apc) + ((Aa * GK2) / CC[Elem].Apa) + ((Ab * GK3) / CC[Elem].Apb);
				CC[Elem].Ka = ((CC[Elem].Area * GK2) / CC[Elem].Apa) + ((Al * GK4) / (Aa + Al)) + ((Ad * GK5) / (Aa+Ad));
				CC[Elem].Kb = (CC[Elem].Area * GK3 / CC[Elem].Apb) + GK18 + (Aj * GK19 / (Ab + Aj));
				CC[Elem].Kc = (CC[Elem].Area * GK1 / CC[Elem].Apc) + (Af * GK11 / (Ac + Af)) + GK12;
				CC[Elem].Kd = Aa * GK5 / (Aa + Ad);
				CC[Elem].Kf = Ac * GK11 / (Ac + Af);
				CC[Elem].Kg = 0.0;
				CC[Elem].Ki = 0.0;
				CC[Elem].Kj = Ab * GK19 / (Ab + Aj);
				CC[Elem].Kl = Aa * GK4 / (Aa + Al);
			}
			
			if(CC[Elem].CellG == UNAVAILABLEDATA && CC[Elem].CellI != UNAVAILABLEDATA)
			{
				CC[Elem].Kp = ((Ac * GK1) / CC[Elem].Apc) + ((Aa * GK2) / CC[Elem].Apa) + ((Ab * GK3) / CC[Elem].Apb);
				CC[Elem].Ka = ((CC[Elem].Area * GK2) / CC[Elem].Apa) + ((Al * GK4) / (Aa + Al)) + ((Ad * GK5) / (Aa+Ad));
				CC[Elem].Kb = (CC[Elem].Area * GK3 / CC[Elem].Apb) + (Ai *  GK18 / (Ab + Ai)) + (Aj * GK19 / (Ab + Aj));
				CC[Elem].Kc = (CC[Elem].Area * GK1 / CC[Elem].Apc) + (Af * GK11 / (Ac + Af)) + GK12;
				CC[Elem].Kd = Aa * GK5 / (Aa + Ad);
				CC[Elem].Kf = Ac * GK11 / (Ac + Af);
				CC[Elem].Kg = 0.0;
				CC[Elem].Ki = Ab * GK18 / (Ab + Ai);
				CC[Elem].Kj = Ab * GK19 / (Ab + Aj);
				CC[Elem].Kl = Aa * GK4 / (Aa + Al);
			}
			if(CC[Elem].CellG != UNAVAILABLEDATA && CC[Elem].CellI == UNAVAILABLEDATA)
			{
				CC[Elem].Kp = ((Ac * GK1) / CC[Elem].Apc) + ((Aa * GK2) / CC[Elem].Apa) + ((Ab * GK3) / CC[Elem].Apb);
				CC[Elem].Ka = ((CC[Elem].Area * GK2) / CC[Elem].Apa) + ((Al * GK4) / (Aa + Al)) + ((Ad * GK5) / (Aa+Ad));
				CC[Elem].Kb = (CC[Elem].Area * GK3 / CC[Elem].Apb) + GK18 + (Aj * GK19 / (Ab + Aj));
				CC[Elem].Kc = (CC[Elem].Area * GK1 / CC[Elem].Apc) + (Af * GK11 / (Ac + Af)) + (Ag * GK12 / (Ag + Ac));
				CC[Elem].Kd = Aa * GK5 / (Aa + Ad);
				CC[Elem].Kf = Ac * GK11 / (Ac + Af);
				CC[Elem].Kg = Ac * GK12 / (Ac + Ag);
				CC[Elem].Ki = 0.0;
				CC[Elem].Kj = Ab * GK19 / (Ab + Aj);
				CC[Elem].Kl = Aa * GK4 / (Aa + Al);
			}
			break;
		case BODY1:
		case INLET:
		case OUTLET:

		case TOP:
		case BOTTOM:
			CC[Elem].Kp = GK1 + ((Aa * GK2) / CC[Elem].Apa) + ((Ab * GK3) / CC[Elem].Apb);//
			CC[Elem].Ka = ((CC[Elem].Area * GK2) / CC[Elem].Apa) + ((Al * GK4) / (Aa + Al)) + ((Ad * GK5) / (Aa+Ad));
			CC[Elem].Kb = (CC[Elem].Area * GK3 / CC[Elem].Apb) + (Ai *  GK18 / (Ab + Ai)) + (Aj * GK19 / (Ab + Aj));
			CC[Elem].Kc = 0.0;//
			CC[Elem].Kd = Aa * GK5 / (Aa + Ad);
			CC[Elem].Kf = 0.0;//
			CC[Elem].Kg = 0.0;//
			CC[Elem].Ki = Ab * GK18 / (Ab + Ai);
			CC[Elem].Kj = Ab * GK19 / (Ab + Aj);
			CC[Elem].Kl = Aa * GK4 / (Aa + Al);
			break;
/*
		
			CC[Elem].Kp = CC[Elem].Ka = CC[Elem].Kb = CC[Elem].Kc = CC[Elem].Kd = CC[Elem].Kf = CC[Elem].Kg = CC[Elem].Ki = CC[Elem].Kj = CC[Elem].Kl = 0.0;
			break;

		
	
			CC[Elem].Kp = CC[Elem].Ka = CC[Elem].Kb = CC[Elem].Kc = CC[Elem].Kd = CC[Elem].Kf = CC[Elem].Kg = CC[Elem].Ki = CC[Elem].Kj = CC[Elem].Kl = 0.0;
			break;
*/		
		default:
			printf("Default Case CC[%d].Type=%d",Elem, CC[Elem].Type);
			//getch();
			exit(0);
		}
		SumK = CC[Elem].Ka + CC[Elem].Kb + CC[Elem].Kc + CC[Elem].Kd + CC[Elem].Kf + CC[Elem].Kg + CC[Elem].Ki + CC[Elem].Kj + CC[Elem].Kl + CC[Elem].Kp;
		//printf("\nCC[%d].SumK = %e",Elem,SumK);
		if(fabsl(SumK) > 1e-5)
		{
			SumK=SumK;
			printf("\nCC[%d].Type = %d",Elem, CC[Elem].Type);
		}
	}

	printf("\nKs Found...\n");
	
	//getch();
	//free(Node);
}

void WritePLTFILE(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, NODE *Node)
{
	FILE *Fwrite;
	int i;

	char Filename[100];
	sprintf(Filename,"PLTPLT.plt");

	if( (Fwrite = fopen(Filename,"wt")) == NULL)
	{
		printf("\nUnable to open Plt File");
		return;
	}

	i=34;

	fprintf(Fwrite,"TITLE = %c2D Flow Past Circular Cylinder%c\n",i,i);
	
	fprintf(Fwrite, "VARIABLES = \"X\", \"Y\",  \"Type\"\n");
	fprintf(Fwrite, "ZONE N=%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED), ZONETYPE=FETRIANGLE\n",NoOfNodes,NoOfElements);

		
	for(i=0;i<NoOfNodes;i++)
	{
		fprintf(Fwrite,"%lf\n", Node[i].x);
	}
	fprintf(Fwrite,"\n\n\n");
	for(i=0;i<NoOfNodes;i++)
	{
		fprintf(Fwrite,"%lf\n", Node[i].y);
	}
	fprintf(Fwrite,"\n\n\n");
	
	for(i=0;i<NoOfElements;i++)
	{
		fprintf(Fwrite,"%d\n", CC[i].Type);
	}
	
	

	fprintf(Fwrite,"\n\n");
	for(i=0; i<NoOfElements; i++)
	{
		fprintf(Fwrite,"%d\t%d\t%d\n", CC[i].Connect[0], CC[i].Connect[1], CC[i].Connect[2] );
	}
	
	fclose(Fwrite);
	//exit(0);
}
void ComputeNodeDistance(int NoOfNodes, CELLDETAILS *CC, NODE *Node)
{
	int Nod;
	double Xc, Yc;
	double Xp, Yp;
	int count;
	for(Nod = 0; Nod<NoOfNodes;Nod++)
	{
		Xp = Node[Nod].x;
		Yp = Node[Nod].y;
		if(Nod==134)
				Nod=Nod;
		for(count=0;count<Node[Nod].ShareCount;count++)
		{
			Xc = CC[Node[Nod].ShareElements[count]].Xp;
			Yc = CC[Node[Nod].ShareElements[count]].Yp;

			Node[Nod].CC2NodeDist[count] = sqrt(pow( (Xp - Xc), 2 ) + pow( (Yp-Yc),2 ));
		}
	}
}
