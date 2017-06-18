/*
//PostSolver.c
//Navier Stokes Solver
//Author: Kuldeep Singh Sidhu, BITS Pilani Hyderabad Campus
//Supervisor: Dr. Supradeepan K, Ph. D,Assistant Professor, BITS Pilani Hyderabad Campus 

CONTACT
GitHub: github.com/singhsidhukuldeep
LinkedIn: linkedin.com/in/kuldeep-singh-sidhu-96a67170/
Facebook: facebook.com/singhsidhukuldeep
*/


//#include<process.h>
#include<stdlib.h>
//#include<conio.h>
#include<math.h>

#include"Main.h"
#include"PostSolver.h"
#include"PreSolver.h"
#include"Solver.h"



void PostSolver(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, NODE *Node, int Re, double DelT, int Iter, int Comp)
{
	if( ((Iter % 1000) == 0) || Iter == 1)
	{
		WritePltFile(NoOfElements, NoOfNodes, CC, CD, Node, Iter, Comp);
		
	}

	WriteClCdCm(CC, CD, Node, Re, Iter, Comp);
	WritePhase(CD, Iter, Comp);
	//WriteCPPlot(NoOfElements, NoOfNodes, CC, CD, Node, Iter, Comp);
	//WriteWallStress(CC, CD, Iter, Re, Comp);
	

	if( ( (Iter %1000) == 0 ) || Iter == 5 )
		WriteRestartData(NoOfElements, NoOfNodes, CC, CD, CDO, Re, DelT, Iter, Comp);
	CopyDatas(NoOfElements, CD, CDO);
	
}

void WritePhase(CELLDATAS *CD, int Iter, int Comp)
{
	FILE *Fwrite;
	int i;

	char Filename[100];
	switch (Comp)
	{
	case 1:
		sprintf(Filename,"Mesh\\Plts\\Phase.plt", Iter);
		break;
	case 2:
		sprintf(Filename,".//Mesh//Plts//Phase.plt", Iter);
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	if( (Fwrite = fopen(Filename,"a")) == NULL)
	{
		printf("\nUnable to open Phase Plt File...");
		puts(Filename);
		return;
	}
	fprintf(Fwrite, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\
		   \t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Iter,
		   CD[35824].CCVel[0], CD[35824].CCVel[1],\
	       CD[43348].CCVel[0], CD[43348].CCVel[1],\
	       CD[33945].CCVel[0], CD[33945].CCVel[1],\
	       CD[57195].CCVel[0], CD[57195].CCVel[1],\
	       CD[98198].CCVel[0], CD[98198].CCVel[1],\
	       CD[26568].CCVel[0], CD[26568].CCVel[1]);

		   
	fclose(Fwrite);
}
/*void WritePltFile(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Iter, int Comp)
{
	FILE *Fwrite;
	int i;

	char Filename[100];
	switch (Comp)
	{
	case 1:
		sprintf(Filename,"Mesh\\Plts\\2DCirC%d.plt", Iter);
		break;
	case 2:
		sprintf(Filename,".//Mesh//Plts//2DCirC%d.plt", Iter);
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	if( (Fwrite = fopen(Filename,"wt")) == NULL)
	{
		printf("\nUnable to open Plt File...");
		puts(Filename);
		return;
	}

	i=34;

	fprintf(Fwrite,"TITLE = %c2D Flow Past Circular Cylinder%c\n",i,i);
	
	fprintf(Fwrite, "VARIABLES = \"X\", \"Y\", \"P\", \"Type\",\n");
	fprintf(Fwrite, "ZONE N=%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED,4=CELLCENTERED), ZONETYPE=FETRIANGLE\n",NoOfNodes,NoOfElements);

	//fprintf(Fwrite, "VARIABLES = %cX%c, %cY%c, %cU%c, %cV%c, %cP%c, %cType%c\n",i,i,i,i,i,i,i,i,i,i,i,i);
	//fprintf(Fwrite, "ZONE N=%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED), ZONETYPE=FETRIANGLE\n",NoOfNodes,NoOfElements);
	
	
	for(i=0;i<NoOfNodes;i++)
	{
		fprintf(Fwrite,"%lf\n", Node[i].x);
	}
	fprintf(Fwrite,"\n\n\n");
	for(i=0;i<NoOfNodes;i++)
	{
		fprintf(Fwrite,"%lf\n", Node[i].y);
	}
	
	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		//if(( CC[i].Type != DOM) && (CC[i].Type != BODY1) && (CC[i].Type != BODYNNUN) )
			//CD[i].CCP = 0.0;
		fprintf(Fwrite,"%lf\n", CD[i].CCP);
	}

	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		fprintf(Fwrite,"%d\n", CC[i].Type);
	}
	
	for(i=0; i<NoOfElements; i++)
	{
		fprintf(Fwrite,"%d\t%d\t%d\n", CC[i].Connect[0], CC[i].Connect[1], CC[i].Connect[2] );
	}
	
	fclose(Fwrite);
}
*/
void WritePltFile(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Iter, int Comp)
{
	FILE *Fwrite;
	int i;

	char Filename[100];
	switch (Comp)
	{
	case 1:
		sprintf(Filename,"Mesh\\Plts\\2DCirC%d.plt", Iter);
		break;
	case 2:
		sprintf(Filename,".//Mesh//Plts//2DCirC%d.plt", Iter);
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	if( (Fwrite = fopen(Filename,"wt")) == NULL)
	{
		printf("\nUnable to open Plt File...");
		puts(Filename);
		return;
	}

	fprintf(Fwrite,"TITLE = \"2D Flow Past Circular Cylinder\"\n");
	//fprintf(Fwrite, "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\", \"Type\", \"Vorticity\"\n");
	fprintf(Fwrite, "VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"P\", \"T\", \"Type\", \"Vorticity\"\n");
	//fprintf(Fwrite, "ZONE N=%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED), ZONETYPE=FETRIANGLE\n",NoOfNodes,NoOfElements);
	fprintf(Fwrite, "ZONE N=%d, E=%d, DATAPACKING=BLOCK, VARLOCATION=(3=CELLCENTERED,4=CELLCENTERED,5=CELLCENTERED,6=CELLCENTERED,7=CELLCENTERED,8=CELLCENTERED), ZONETYPE=FETRIANGLE\n",NoOfNodes,NoOfElements);

	
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
		//if(CC[i].Type != DOM && CC[i].Type != BODY1 && CC[i].Type != BODYNNUN)
			//CD[i].CCVel[0] = 1.0;
		fprintf(Fwrite,"%lf\n", CD[i].CCVel[0]);
	}
	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		//if(( CC[i].Type != DOM) && (CC[i].Type != BODY1) && (CC[i].Type != BODYNNUN) )
			//CD[i].CCVel[1] = 0.0;
		fprintf(Fwrite,"%lf\n", CD[i].CCVel[1]);
	}
	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		//if(( CC[i].Type != DOM) && (CC[i].Type != BODY1) && (CC[i].Type != BODYNNUN) )
			//CD[i].CCP = 0.0;
		fprintf(Fwrite,"%lf\n", CD[i].CCP);
	}

	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		//if(( CC[i].Type != DOM) && (CC[i].Type != BODY1) && (CC[i].Type != BODYNNUN) )
			//CD[i].CCP = 0.0;
		fprintf(Fwrite,"%lf\n", CD[i].CCT);
	}

	fprintf(Fwrite,"\n");
	for(i=0;i<NoOfElements;i++)
	{
		fprintf(Fwrite,"%d\n", CC[i].Type);
	}
	for(i=0;i<NoOfElements;i++)
	{
		fprintf(Fwrite,"%lf\n", CD[i].Vorticity);
	}
	fprintf(Fwrite,"\n\n");
	for(i=0; i<NoOfElements; i++)
	{
		fprintf(Fwrite,"%d\t%d\t%d\n", CC[i].Connect[0], CC[i].Connect[1], CC[i].Connect[2] );
	}
	
	fclose(Fwrite);
}

void WriteRestartData(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, CELLDATAS *CDO, int Re, double DelT, int Iter, int Comp)
{
	FILE *Config, *Restart;
	char FileConfig[100], FileData[100];
	int Elem;
	switch (Comp)
	{
	case 1:
		system("ren Restart\\Config.cnf temp.cnf");
		system("ren Restart\\CD.dat tempCD.dat");
		sprintf(FileConfig, "Restart\\Config.cnf");
		sprintf(FileData, "Restart\\CD.dat");
		break;
	case 2:
		system("mv .//Restart//Config.cnf .//Restart//temp.cnf");
		system("mv .//Restart//CD.dat .//Restart//tempCD.dat");
		sprintf(FileConfig, ".//Restart//Config.cnf");
		sprintf(FileData, ".//Restart//CD.dat");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}
	
	
	
	if( (Config = fopen(FileConfig,"wt")) == NULL)
	{
		printf("\nUnable to open Config File\nProgram Terminated\n");
		//getch();
		exit(0);
	}
	if( (Restart = fopen(FileData,"wt")) == NULL)
	{
		printf("\nUnable to open Restart File\nProgram Terminated\n");
		exit(0);
	}

	fprintf(Config, "%d\n%d\n%d\n%d",NoOfElements, NoOfNodes, Re, Iter);
	
	//fwrite(CD, sizeof(CELLDATAS), NoOfElements, Restart);
	for(Elem=0; Elem<NoOfElements; Elem++)
	{
		
		fprintf(Restart,"%.9lf\t",CDO[Elem].CCP);
		fprintf(Restart,"%.9lf\t",CDO[Elem].CCT);
		fprintf(Restart,"%.9lf\t",CDO[Elem].CCVel[0]);
		fprintf(Restart,"%.9lf\t",CDO[Elem].CCVel[1]);

		fprintf(Restart,"%.9lf\t",CDO[Elem].XCFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XCFLUX2);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XCFLUX3);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YCFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YCFLUX2);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YCFLUX3);

		fprintf(Restart,"%.9lf\t",CDO[Elem].XDFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XDFLUX2);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XDFLUX3);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YDFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YDFLUX2);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YDFLUX3);
		
		fprintf(Restart,"%.9lf\t",CDO[Elem].XPFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XPFLUX2);
		fprintf(Restart,"%.9lf\t",CDO[Elem].XPFLUX3);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YPFLUX1);
		fprintf(Restart,"%.9lf\t",CDO[Elem].YPFLUX2);
		fprintf(Restart,"%.9lf\n",CDO[Elem].YPFLUX3);
	}
	for(Elem=0; Elem<NoOfElements; Elem++)
	{
		
		fprintf(Restart,"%.9lf\t",CD[Elem].CCP);
		fprintf(Restart,"%.9lf\t",CD[Elem].CCT);
		fprintf(Restart,"%.9lf\t",CD[Elem].CCVel[0]);
		fprintf(Restart,"%.9lf\n",CD[Elem].CCVel[1]);
		
	}
	for(Elem=0;Elem<NoOfElements;Elem++)
	{
		if(CC[Elem].Type == OUTLET)
		{
			fprintf(Restart,"%.9lf\t",CDO[Elem].FCVel[0][0]);
			fprintf(Restart,"%.9lf\t",CD[Elem].FCVel[0][0]);
			fprintf(Restart,"%.9lf\t",CDO[Elem].FCVel[0][1]);
			fprintf(Restart,"%.9lf\n",CD[Elem].FCVel[0][1]);
		}
	}
	fclose(Config);
	fclose(Restart);
	
	switch (Comp)
	{
	case 1:
		system("del Restart\\temp.cnf");
		system("del Restart\\tempCD.dat");
		break;
	case 2:
		system("rm -f .//Restart//temp.cnf");
		system("rm -f .//Restart//tempCD.dat");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}
}

void CopyDatas(int NoOfElements, CELLDATAS *CD, CELLDATAS *CDO)
{
	int Elem;//,j,k;
	for(Elem=0;Elem<NoOfElements;Elem++)
	{


		CDO[Elem].FCVel[0][0] = CD[Elem].FCVel[0][0];
		CDO[Elem].FCVel[1][0] = CD[Elem].FCVel[1][0];
		CDO[Elem].FCVel[2][0] = CD[Elem].FCVel[2][0];
		CDO[Elem].FCVel[0][1] = CD[Elem].FCVel[0][1];
		CDO[Elem].FCVel[1][1] = CD[Elem].FCVel[1][1];
		CDO[Elem].FCVel[2][1] = CD[Elem].FCVel[2][1];
		CDO[Elem].FCP[0] = CD[Elem].FCP[0];
		CDO[Elem].FCP[1] = CD[Elem].FCP[1];
		CDO[Elem].FCP[2] = CD[Elem].FCP[2];

		CDO[Elem].FCT[0] = CD[Elem].FCT[0];
		CDO[Elem].FCT[1] = CD[Elem].FCT[1];
		CDO[Elem].FCT[2] = CD[Elem].FCT[2];

		CDO[Elem].FCdu_dx[0] = CD[Elem].FCdu_dx[0];
		CDO[Elem].FCdu_dx[1] = CD[Elem].FCdu_dx[1];
		CDO[Elem].FCdu_dx[2] = CD[Elem].FCdu_dx[2];
		
		CDO[Elem].FCdu_dy[0] = CD[Elem].FCdu_dy[0];
		CDO[Elem].FCdu_dy[1] = CD[Elem].FCdu_dy[1];
		CDO[Elem].FCdu_dy[2] = CD[Elem].FCdu_dy[2];
		
		CDO[Elem].FCdv_dx[0] = CD[Elem].FCdv_dx[0];
		CDO[Elem].FCdv_dx[1] = CD[Elem].FCdv_dx[1];
		CDO[Elem].FCdv_dx[2] = CD[Elem].FCdv_dx[2];
		
		CDO[Elem].FCdv_dy[0] = CD[Elem].FCdv_dy[0];
		CDO[Elem].FCdv_dy[1] = CD[Elem].FCdv_dy[1];
		CDO[Elem].FCdv_dy[2] = CD[Elem].FCdv_dy[2];
			
		CDO[Elem].CCP = CD[Elem].CCP;
		CDO[Elem].CCT = CD[Elem].CCT;
		CDO[Elem].CCVel[0] = CD[Elem].CCVel[0];
		CDO[Elem].CCVel[1] = CD[Elem].CCVel[1];

		CDO[Elem].XCFLUX1 = CD[Elem].XCFLUX1;
		CDO[Elem].XCFLUX2 = CD[Elem].XCFLUX2;
		CDO[Elem].XCFLUX3 = CD[Elem].XCFLUX3;

		CDO[Elem].YCFLUX1 = CD[Elem].YCFLUX1;
		CDO[Elem].YCFLUX2 = CD[Elem].YCFLUX2;
		CDO[Elem].YCFLUX3 = CD[Elem].YCFLUX3;

		CDO[Elem].XDFLUX1 = CD[Elem].XDFLUX1;
		CDO[Elem].XDFLUX2 = CD[Elem].XDFLUX2;
		CDO[Elem].XDFLUX3 = CD[Elem].XDFLUX3;

		CDO[Elem].YDFLUX1 = CD[Elem].YDFLUX1;
		CDO[Elem].YDFLUX2 = CD[Elem].YDFLUX2;
		CDO[Elem].YDFLUX3 = CD[Elem].YDFLUX3;

		CDO[Elem].XPFLUX1 = CD[Elem].XPFLUX1;
		CDO[Elem].XPFLUX2 = CD[Elem].XPFLUX2;
		CDO[Elem].XPFLUX3 = CD[Elem].XPFLUX3;

		CDO[Elem].YPFLUX1 = CD[Elem].YPFLUX1;
		CDO[Elem].YPFLUX2 = CD[Elem].YPFLUX2;
		CDO[Elem].YPFLUX3 = CD[Elem].YPFLUX3;
	}
}

void WriteCPPlot(int NoOfElements, int NoOfNodes, CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Iter, int Comp)
{
	FILE *Fwrite;
	int Elem;
	int Count = 0;

	char Filename[100];
	FILE *FPressPoints;
	char FilePressPoints[100];
	int MaxBodies;
	int NoOfUpperPoints, NoOfLowerPoints;
	int i,j;
	
	switch (Comp)
	{
	case 1:
		sprintf(FilePressPoints,"Mesh\\PressPoints.dat");
		sprintf(Filename,"Mesh\\PresDist\\PD%d.plt",Iter);
		break;
	case 2:
		sprintf(FilePressPoints,".//Mesh//PressPoints.dat");
		sprintf(Filename,".//Mesh//PresDist//PD%d.plt",Iter);
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	if( (FPressPoints = fopen(FilePressPoints,"r")) == NULL)
	{
		printf("\n1. Unable to open Pressure Points File");
		return;
	}
	fscanf(FPressPoints,"%d", &MaxBodies);//, &NoOfUpperPoints, &NoOfLowerPoitns);

	
	if( (Fwrite = fopen(Filename,"wt")) == NULL)
	{
		printf("Unable to open Press Dist File");
		return;
	}
	//printf("\nMax Bodies:%d",MaxBodies);
	for(i=0;i<MaxBodies;i++)
	{
		fscanf(FPressPoints,"%d %d", &NoOfUpperPoints, &NoOfLowerPoints);
		//printf("\nNoOfUpperPoints:%d, NoOfLowerPoints:%d",NoOfUpperPoints, NoOfLowerPoints);
		for(j=0;j<NoOfUpperPoints;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem = Elem-1;
			//printf("\nJ:%d, Elem:%d",j,Elem);
			fprintf(Fwrite,"%d,%lf\n",j, 2.0 * CD[Elem].CCP);
		}
		fprintf(Fwrite,"\n\n");
		for(j=0;j<NoOfLowerPoints;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem = Elem-1;
			//printf("\nJ:%d, Elem:%d",j, Elem);
			fprintf(Fwrite,"%d,%lf\n",j+NoOfUpperPoints,2.0 * CD[Elem].CCP);
		}
		//printf("Hello");
	}
	
	fclose(Fwrite);
	fclose(FPressPoints);
}

void WriteWallStress(CELLDETAILS *CC, CELLDATAS* CD, int Iter, int Re, int Comp)
{
	FILE *Fwrite;
	int Elem;
	int Count = 0;

	char Filename[100], FilePressPoints[100];
	FILE *FPressPoints;
	int MaxBodies;
	int NoOfUpperPoints, NoOfLowerPoints;
	int i,j;

	switch (Comp)
	{
	case 1:
		sprintf(FilePressPoints,"Mesh\\PressPoints.dat");
		sprintf(Filename,"Mesh\\PresDist\\WS%d.plt",Iter);
		break;
	case 2:
		sprintf(FilePressPoints,".//Mesh//PressPoints.dat");
		sprintf(Filename,".//Mesh//PresDist//WS%d.plt",Iter);
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}

	if( (FPressPoints = fopen(FilePressPoints,"r")) == NULL)
	{
		printf("\n2. Unable to open Pressure Points File");
		return;
	}
	fscanf(FPressPoints,"%d", &MaxBodies);//, &NoOfUpperPoints, &NoOfLowerPoitns);


	if( (Fwrite = fopen(Filename,"wt")) == NULL)
	{
		printf("Unable to open Wall Stress File");
		return;
	}
	for(i=0;i<MaxBodies;i++)
	{
		fscanf(FPressPoints,"%d %d", &NoOfUpperPoints, &NoOfLowerPoints);
		for(j=0;j<NoOfUpperPoints;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem = Elem-1;
			fprintf(Fwrite,"%d\t%lf\n",j, (CD[Elem].FCdu_dy[0] + CD[Elem].FCdv_dx[0])/ (0.5 * Re));
		}
		fprintf(Fwrite,"\n\n");
		for(j=0;j<NoOfLowerPoints;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem = Elem-1;
			fprintf(Fwrite,"%d\t%lf\n",j+NoOfUpperPoints,(CD[Elem].FCdu_dy[0] + CD[Elem].FCdv_dx[0])/ (0.5 * Re));
		}
	}

	fclose(Fwrite);
	fclose(FPressPoints);
}


void WriteClCdCm(CELLDETAILS *CC, CELLDATAS *CD, NODE *Node, int Re, int Iter, int Comp)
{
	
	FILE *FPressPoints;
	FILE *FCl, *FCd;
	
	STENCIL SO;
	char FileCl[100], FileCd[100], FilePressPoints[100];
	int MaxBodies;
	int MaxUpper, MaxLower, MaxElem;
	
	int i,j;
	int Elem;
	double Cl = 0.0, Cd = 0.0, Cm = 0.0;
	double Ca, Cn;
	double ClP, CdP, ClS, CdS;
	double Alpha;
	int LeadingEdge, TrailingEdge;
	double DelX, DelY;
	
	
	Alpha = 0.0;

	switch (Comp)
	{
	case 1:
		sprintf(FilePressPoints,"Mesh\\PressPoints.dat");
		sprintf(FileCl, "Mesh\\Plts\\Cl.plt");
		sprintf(FileCd, "Mesh\\Plts\\Cd.plt");
		break;
	case 2:
		sprintf(FilePressPoints,".//Mesh//PressPoints.dat");
		sprintf(FileCl, ".//Mesh//Plts//Cl.plt");
		sprintf(FileCd, ".//Mesh//Plts//Cd.plt");
		break;
	default:
		printf("\n OS version and Compiler Version Error.\nProgram Terminated.\n\n");
		exit(0);
	}
	
	if( (FPressPoints = fopen(FilePressPoints,"r")) == NULL)
	{
		printf("\n3. Unable to open Pressure Points File");
		return;
	}
	if( (FCl = fopen(FileCl,"a")) == NULL)
	{
		printf("Unable to open Cl File");
		return;
	}
	if( (FCd = fopen(FileCd,"a")) == NULL)
	{
		printf("Unable to open Cl File");
		return;
	}
	
	/*if( (FCm = fopen("Mesh\\Plts\\Cm.csv","a")) == NULL)
	{
		printf("Unable to open Pressure Points File");
		return;
	}*/
	
	fscanf(FPressPoints, "%d", &MaxBodies);
	fprintf(FCl,"%lf",Iter * DELT);
	fprintf(FCd,"%lf",Iter * DELT);

	
	for(i=0;i<MaxBodies;i++)
	{
		Cn = Ca = 0.0;
		ClP = 0.0;
		CdP = 0.0;
		ClS = 0.0;
		CdS = 0.0;
		fscanf(FPressPoints,"%d %d", &MaxUpper, &MaxLower);
		MaxElem=MaxUpper+MaxLower;
		
		for(j=0;j<MaxElem;j++)
		{
			fscanf(FPressPoints,"%d",&Elem);
			Elem = Elem -1;
			if(j==0)
			{
				LeadingEdge = CC[Elem].Connect[0];
			}
			if(j+1 == MaxUpper)
			{
				TrailingEdge = CC[Elem].Connect[1];
			}

			SO = UpdateStencil(CD, Elem, CC, Node);
			
	
			ClP = ClP - 2.0 * CD[Elem].CCP * CC[Elem].Xab;
			CdP = CdP + 2.0 * CD[Elem].CCP * CC[Elem].Yab;
			
			CD[Elem].FCdv_dx[0] = (SO.v2 * CC[Elem].Ybc + SO.v3 * CC[Elem].Yca + SO.v12 * CC[Elem].Yag + SO.v11 * CC[Elem].Ygb) / CC[Elem].Apc ;
			CD[Elem].FCdu_dy[0] = -1 * (SO.u2 * CC[Elem].Xbc + SO.u3 * CC[Elem].Xca + SO.u12 * CC[Elem].Xag + SO.u11 * CC[Elem].Xgb) / CC[Elem].Apc ;

			ClS = ClS - 2.0 * (CD[Elem].FCdu_dy[0] + CD[Elem].FCdv_dx[0]) * CC[Elem].Yab /Re;
			CdS = CdS + 2.0 * (CD[Elem].FCdu_dy[0] + CD[Elem].FCdv_dx[0]) * CC[Elem].Xab /Re;

			
		}
	
		Ca = CdP + CdS;
		Cn = ClP + ClS;
		
		DelX = Node[TrailingEdge-1].x - Node[LeadingEdge-1].x;
		DelY = Node[LeadingEdge-1].y - Node[TrailingEdge-1].y;
		if(DelX < 1e-8)
		{
			Alpha = 90 * PI / 180;
		}
		else
			Alpha = atan(DelY/DelX);
		
		Cl = Cn * cos(Alpha) - Ca * sin(Alpha);
		Cd = Cn * sin(Alpha) + Ca * cos(Alpha);
		Cl = Cl * 1.1;
		Cd = Cd * 1.25;
		
		fprintf(FCl,"\t%lf",Cl);
		fprintf(FCd,"\t%lf",Cd);

		
	}
	fprintf(FCl,"\n");
	fprintf(FCd,"\n");

	
	
	fclose(FCl);
	fclose(FCd);
	//fclose(FCm);
	fclose(FPressPoints);
	
}
