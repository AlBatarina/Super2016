#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>


// Domain size.
const double A = 2.0;
const double B = 2.0;

const double A0 = -1.0;
const double B0 = -1.0;

const double A1 = 1.0;
const double B1 = 1.0;

const int SDINum = 1;

#define Test
#define Print
#define TRUE  ((int) 1)
#define FALSE ((int) 0)
#define Step 100

#define Max(A,B) ((A)>(B)?(A):(B))
#define Min(A,B) ((A)<(B)?(A):(B))
#define R2(x,y) ((x)*(x)+(y)*(y))
#define Cube(x) ((x)*(x)*(x))
#define Sqr(x) ((x)*(x))

#define hx(i,XNodes)  (XNodes[i+1]-XNodes[i])
#define hy(j,YNodes)  (YNodes[j+1]-YNodes[j])

#define LeftPart(P,i,j,XNodes,YNodes,NX,NY)\
((-(P[NX*(j)+i+1]-P[NX*(j)+i])/hx(i,XNodes)+(P[NX*(j)+i]-P[NX*(j)+i-1])/hx(i-1,XNodes))/(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))+\
(-(P[NX*(j+1)+i]-P[NX*(j)+i])/hy(j,YNodes)+(P[NX*(j)+i]-P[NX*(j-1)+i])/hy(j-1,YNodes))/(0.5*(hy(j,YNodes)+hy(j-1,YNodes))))

double Solution(double x,double y)
{
    return Sqr(1-x*x)+Sqr(1-y*y);	//log(1 + x*y);
}

double BoundaryValue(double x, double y)
{
	return Solution(x,y);
}

int RightPart(double * rhs, double * XNodes, double * YNodes, int NX, int NY)
{
    int i, j;

	for(j=0; j<NY; j++)
	   for(i=0; i<NX; i++)
	       rhs[j*NX+i] = 4*(2-3*Sqr(XNodes[i])-3*Sqr(YNodes[j]));	//(Sqr(XNodes[i])+Sqr(YNodes[j]))/Sqr(1+XNodes[i]*YNodes[j]);
	return 0;
}

int MeshGenerate(double * XNodes, double * YNodes, int N0, int N1)
{
	//const double q = 1.5;
	int i;
	double hx = A / (N0-1);
	double hy = B / (N1-1);

	for(i=0; i<N0; i++)
		XNodes[i] = A0 + i*hx;	//XNodes[i] = A*(pow(1.0+i/(N0-1.0),q)-1.0)/(pow(2.0,q)-1.0);
	for(i=0; i<N1; i++)
		YNodes[i] = B0 + i*hy;	//YNodes[i] = B*(pow(1.0+i/(N1-1.0),q)-1.0)/(pow(2.0,q)-1.0);
	XNodes[N0-1]=A1;
	YNodes[N1-1]=B1;
	return 0;
}

int IsPower(int Number)
// the function returns log_{2}(Number) if it is integer. If not it returns (-1). 
{
    unsigned int M;
    int p;
    
    if(Number <= 0)
        return(-1);
        
    M = Number; p = 0;
    while(M % 2 == 0)
    {
        ++p;
        M = M >> 1;
    }
    if((M >> 1) != 0)
        return(-1);
    else
        return(p);
}

int SplitFunction(int N0, int N1, int p)
// This is the splitting procedure of proc. number p. The integer p0
// is calculated such that abs(N0/p0 - N1/(p-p0)) --> min.
{
    float n0, n1;
    int p0, i;
    
    n0 = (float) N0; n1 = (float) N1;
    p0 = 0;
    
    for(i = 0; i < p; i++)
        if(n0 > n1)
        {
            n0 = n0 / 2.0;
            ++p0;
        }
        else
            n1 = n1 / 2.0;
    
    return(p0);
}

int Send(double * Vect, int NX, int NY, int up, int down, int left, int right, int tag, MPI_Comm Grid_Comm){
	int j;
	double * TmpVect= (double *)malloc((NY-2)*sizeof(double));						// for sending or receiving

	if (up >= 0) MPI_Send(Vect+NX+1, NX-2, MPI_DOUBLE, up, tag, Grid_Comm);
	if (down >= 0) MPI_Send(Vect+NX*(NY-2)+1, NX-2, MPI_DOUBLE, down, tag, Grid_Comm);

	if (left >= 0){
		for(j=1; j < NY-1; j++)
			TmpVect[j-1] = Vect[NX*j+1];
		MPI_Send(TmpVect, NY-2, MPI_DOUBLE, left, tag, Grid_Comm);
	}

	if (right >= 0){
		for(j=1; j < NY-1; j++)
			TmpVect[j-1] = Vect[NX*j+NX-2];
		MPI_Send(TmpVect, NY-2, MPI_DOUBLE, right, tag, Grid_Comm);
	}
	free(TmpVect);
	return 0;
}

int Receive(double * Vect, int NX, int NY, int up, int down, int left, int right, int tag, MPI_Comm Grid_Comm){
	int j;
    MPI_Status status;
	double * TmpVect= (double *)malloc((NY-2)*sizeof(double));						// for sending or receiving

	if (up >= 0) MPI_Recv(Vect+1, NX-2, MPI_DOUBLE, up, tag, Grid_Comm, &status);
	if (down >= 0) MPI_Recv(Vect+NX*(NY-1)+1, NX-2, MPI_DOUBLE, down, tag, Grid_Comm, &status);

	if (left >= 0){
		MPI_Recv(TmpVect, NY-2, MPI_DOUBLE, left, tag, Grid_Comm, &status);
		for(j=1; j < NY-1; j++)
			Vect[NX*j] = TmpVect[j-1];
	}

	if (right >= 0){
		MPI_Recv(TmpVect, NY-2, MPI_DOUBLE, right, tag, Grid_Comm, &status);
		for(j=1; j < NY-1; j++)
			Vect[NX*j+NX-1] = TmpVect[j-1];
	}	
	
	free(TmpVect);
	return 0;
}

int CGM(int CGMNum, double * XNodes, double * YNodes, int NX, int NY, int N0, int N1, MPI_Comm Grid_Comm){

	double * SolVect;						// the solution array.
	double * ResVect;						// the residual array.
	double * BasisVect;						// the vector of A-orthogonal system in CGM.
	double * RHS_Vect;						// the right hand side of Puasson equation.
	double sp, alpha, tau, NewValue, err;	// auxiliary values.
	int counter;							// the current iteration number.
	double tmp;
    int left, right, up, down;      // the neighbours of the process.
    int rank;
    int ProcNum;

	int i,j;
	char str[127];
	FILE * fp;
	int x0, xn, y0, yn;						// internal boundaries


    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Cart_shift(Grid_Comm, 0, 1, &left, &right);
    MPI_Cart_shift(Grid_Comm, 1, 1, &down, &up);

	SolVect   = (double *)malloc(NX*NY*sizeof(double));
	ResVect   = (double *)malloc(NX*NY*sizeof(double));
	RHS_Vect  = (double *)malloc(NX*NY*sizeof(double));

// Initialization of Arrays
	memset(ResVect,0,NX*NY*sizeof(double));
	memset(SolVect,1,NX*NY*sizeof(double));
	RightPart(RHS_Vect,XNodes,YNodes,NX,NY);

	x0=y0=0;
	xn=NX; yn=NY;

	if (down < 0){
		for(i=0; i<NX; i++)
			SolVect[i] = BoundaryValue(XNodes[i],B0);
		y0=1;
	} 		
	if (up < 0){
		for(i=0; i<NX; i++)
			SolVect[NX*(NY-1)+i] = BoundaryValue(XNodes[i],B1);
		yn=NY-1;
	}   
	if (left < 0){
		for(j=0; j<NY; j++)
			SolVect[NX*j] = BoundaryValue(A0,YNodes[j]);
		x0=1;
	}
	if (right < 0){
		for(j=0; j<NY; j++)
			SolVect[NX*j+(NX-1)] = BoundaryValue(A1,YNodes[j]);
		xn=NX-1;
	}

// Steep descent iterations begin ...
	#ifdef Print
		if (rank == 0) printf("\nSteep descent iterations begin ...\n");
	#endif

	for(counter=0; counter<SDINum; counter++)
	{
// The residual vector r(k) = Ax(k)-f is calculating ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				ResVect[NX*j+i] = LeftPart(SolVect,i,j,XNodes,YNodes,NX,NY)-RHS_Vect[NX*j+i];

// Send ResVect to neighbours
		Send(ResVect, NX, NY, up, down, left, right, 2, Grid_Comm);

// The value of product (r(k),r(k)) is calculating ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += ResVect[NX*j+i]*ResVect[NX*j+i]*(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))*(0.5*(hy(j,YNodes)+hy(j-1,YNodes)));
		tmp = sp;
		MPI_Allreduce(&tmp, &sp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = sp;

// Receive ResVect from neighbours
		Receive(ResVect, NX, NY, up, down, left, right, 2, Grid_Comm);

// The value of product sp = (Ar(k),r(k)) is calculating ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += LeftPart(ResVect,i,j,XNodes,YNodes,NX,NY)*ResVect[NX*j+i]*(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))*(0.5*(hy(j,YNodes)+hy(j-1,YNodes)));
		tmp = sp;
		MPI_Allreduce(&tmp, &sp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = tau/sp;

// The x(k+1) is calculating ...
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
			{
				NewValue = SolVect[NX*j+i]-tau*ResVect[NX*j+i];
				err += Sqr(NewValue-SolVect[NX*j+i]);
				SolVect[NX*j+i] = NewValue;
			}

// Send boundary SolVect
		Send(SolVect, NX, NY, up, down, left, right, 1, Grid_Comm);

		tmp = err;
		MPI_Allreduce(&tmp, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		err = sqrt(err);


        if(1)//counter%Step == 0)
        {
            if (rank == 0) printf("The Steep Descent iteration %d has been performed.\n",counter);
            
    #ifdef Print
            if (rank == 0) printf("\nThe Steep Descent iteration k = %d has been performed.\n"
					"Step \\tau(k) = %f.\nThe difference value is estimated by %.12f.\n",\
					counter, tau, err);
    #endif

/*    #ifdef Test
			err = 0.0;
			for(j=1; j < NY-1; j++)
				for(i=1; i < NX-1; i++)
					err = Max(err, fabs(Solution(XNodes[i],YNodes[j])-SolVect[NX*j+i]));
			printf("The Steep Descent iteration %d have been performed. "
					   "The residual error is %.12f\n", counter, err);
    #endif*/
        }
    }
// the end of steep descent iteration.
	if (rank == 0) printf("\nSteep descent iteration ended\n\n");

	BasisVect = ResVect;    // g(0) = r(k-1).
	ResVect = (double *)malloc(NX*NY*sizeof(double));
	memset(ResVect,0,NX*NY*sizeof(double));

// CGM iterations begin ...
// sp == (Ar(k-1),r(k-1)) == (Ag(0),g(0)), k=1.
	#ifdef Print
		if (rank == 0) printf("\nCGM iterations begin ...\n");
	#endif

	for(counter=1; counter<=CGMNum; counter++)
	{
	// Receive SolVect from neighbours
		Receive(SolVect, NX, NY, up, down, left, right, 1, Grid_Comm);

	// The residual vector r(k) is calculating ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				ResVect[NX*j+i] = LeftPart(SolVect,i,j,XNodes,YNodes,NX,NY)-RHS_Vect[NX*j+i];

	// Send ResVect to neighbours
		Send(ResVect, NX, NY, up, down, left, right, 2, Grid_Comm);

	// Receive ResVect from neighbours
		Receive(ResVect, NX, NY, up, down, left, right, 2, Grid_Comm);

	// The value of product (Ar(k),g(k-1)) is calculating ...
		alpha = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				alpha += LeftPart(ResVect,i,j,XNodes,YNodes,NX,NY)*BasisVect[NX*j+i]*(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))*(0.5*(hy(j,YNodes)+hy(j-1,YNodes)));
		tmp = alpha;
		MPI_Allreduce(&tmp, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		alpha = alpha/sp;

	// The new basis vector g(k) is being calculated ...
		for(j=0; j < NY; j++)
			for(i=0; i < NX; i++)
				BasisVect[NX*j+i] = ResVect[NX*j+i]-alpha*BasisVect[NX*j+i];

	// The value of product (r(k),g(k)) is being calculated ...
		tau = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				tau += ResVect[NX*j+i]*BasisVect[NX*j+i]*(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))*(0.5*(hy(j,YNodes)+hy(j-1,YNodes)));
		tmp = tau;
		MPI_Allreduce(&tmp, &tau, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		
	// The value of product sp = (Ag(k),g(k)) is being calculated ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += LeftPart(BasisVect,i,j,XNodes,YNodes,NX,NY)*BasisVect[NX*j+i]*(0.5*(hx(i,XNodes)+hx(i-1,XNodes)))*(0.5*(hy(j,YNodes)+hy(j-1,YNodes)));
		tmp = sp;
		MPI_Allreduce(&tmp, &sp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		tau = tau/sp;

	// The x(k+1) is being calculated ...,NX,NY
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
			{
				NewValue = SolVect[NX*j+i]-tau*BasisVect[NX*j+i];
				err += Sqr(NewValue-SolVect[NX*j+i]);
				SolVect[NX*j+i] = NewValue;
			}
// Send boundary SolVect
		Send(SolVect, NX, NY, up, down, left, right, 1, Grid_Comm);

// Calculate err
		tmp = err;
		MPI_Allreduce(&tmp, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		err = sqrt(err);

		if(counter%Step == 0 && rank == 0)
		{
			printf("The %d iteration of CGM method has been carried out.\n", counter);
            
#ifdef Print            
            printf("\nThe iteration %d of conjugate gradient method has been finished.\n"
                       "The value of \\alpha(k) = %f, \\tau(k) = %f. The difference value is %f.\n",\
                        counter, alpha, tau, err);
#endif
/*#ifdef Test
			tmp = 0.0;
			for(j=1; j < NY-1; j++)
				for(i=1; i < NX-1; i++)
					tmp = Max(err, fabs(Solution(XNodes[i],YNodes[j])-SolVect[NX*j+i]));
			printf("The %d iteration of CGM have been performed. The residual error is %.12f\n",\
                        counter, tmp);
#endif*/
		}
        if (err < 0.00001) break;
	}
// the end of CGM iterations.

// printing some results ...
	if (rank == 0) printf("\nThe %d iterations are carried out. The error of iterations is estimated by %.12f.\n", CGMNum, err);

    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);

	sprintf(str,"PuassonSerial_ECGM_%d_%dx%d.dat", ProcNum, N0, N1);
	fp = fopen(str,"a");
		for (j=0; j < NY; j+=100)
		{
			for (i=0; i < NX; i+=100)
				fprintf(fp,"%f %f %f\n", XNodes[i]+1, YNodes[j]+1, SolVect[NX*j+i]);
		}
	fclose(fp);

	MPI_Barrier(MPI_COMM_WORLD);
	free(SolVect); free(ResVect); free(BasisVect); free(RHS_Vect);
	return(0);
}

int main(int argc, char **argv)
{
	int N0, N1;                     // Mesh has N0 x N1 nodes.
    int ProcNum, rank;              // the number of processes and rank in communicator.
    int power, p0, p1;              // ProcNum = 2^(power), power splits into sum p0 + p1.
    int dims[2];                    // dims[0] = 2^p0, dims[1] = 2^p1 (--> M = dims[0]*dims[1]).
    int n0,n1, k0,k1;               // N0 = n0*dims[0] + k0, N1 = n1*dims[1] + k1.
    int Coords[2];                  // the process coordinates in the cartesian topology created for mesh.
    
    MPI_Comm Grid_Comm;             // this is a handler of a new communicator.
    const int ndims = 2;            // the number of a process topology dimensions.
    int periods[2] = {0,0};         // it is used for creating processes topology.

	int CGMNum;						// the number of CGM iterations.

	int NX, NY;						// the number of internal points on axes (ox) and (oy).
	int i0, j0;						// Initial mesh coordinates for the process
	double * XNodes, * YNodes;		// mesh node coords are stored here.
	double time0, time1;			// for report
    int left, right, up, down;      // the neighbours of the process.

    // command line analizer
	switch (argc)
	{
	case 4:{
				CGMNum = atoi(argv[3]);
				break;
		   }
	default:{
				printf("Wrong number of parameters in command line.\nUsage: <ProgName> "
                       "<Nodes number on (0x) axis> <Nodes number on (0y) axis> "
                       "<the number of conjugate gragient iterations>\nFinishing...\n");
				return(-1);
			}
	}
    
    N0 = atoi(argv[1]); N1 = atoi(argv[2]);
    
    if((N0 <= 0)||(N1 <= 0))
    {
        if(rank == 0)
            printf("The first and the second arguments (mesh numbers) should be positive.\n");
        
        MPI_Finalize();
        return(2);
    }

    // MPI Library is being activated ...
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if((power = IsPower(ProcNum)) < 0)
    {
        if(rank == 0)
            printf("The number of procs must be a power of 2.\n");
        return(3);
    }

    p0 = SplitFunction(N0, N1, power);
    p1 = power - p0;
    
    dims[0] = (unsigned int) 1 << p0;   dims[1] = (unsigned int) 1 << p1;
    n0 = N0 >> p0;                      n1 = N1 >> p1;
    k0 = N0 - dims[0]*n0;               k1 = N1 - dims[1]*n1;

#ifdef Print
    if(rank == 0)
    {
        printf("The number of processes ProcNum = 2^%d. It is split into %d x %d processes.\n"
               "The number of nodes N0 = %d, N1 = %d. Blocks B(i,j) have size:\n", power, dims[0],dims[1], N0,N1);

	if((k0 > 0)&&(k1 > 0))
	    printf("-->\t %d x %d iff i = 0 .. %d, j = 0 .. %d;\n", n0+1,n1+1, k0-1,k1-1);
        if(k1 > 0)
            printf("-->\t %d x %d iff i = %d .. %d, j = 0 .. %d;\n", n0,n1+1, k0,dims[0]-1, k1-1);
        if(k0 > 0)
            printf("-->\t %d x %d iff i = 0 .. %d, j = %d .. %d;\n", n0+1,n1, k0-1, k1,dims[1]-1);

        printf("-->\t %d x %d iff i = %d .. %d, j = %d .. %d.\n", n0,n1, k0,dims[0]-1, k1,dims[1]-1);
    }
#endif

    // the cartesian topology of processes is being created ...
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, TRUE, &Grid_Comm);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Cart_coords(Grid_Comm, rank, ndims, Coords);

    MPI_Cart_shift(Grid_Comm, 0, 1, &left, &right);
    MPI_Cart_shift(Grid_Comm, 1, 1, &down, &up);

    NX=n0; NY=n1;

    if(Coords[0] < k0)
        ++NX;
    if(Coords[1] < k1)
        ++NY;

    i0=Max(n0*Coords[0] + Min(Coords[0],k0) - 1,0);
    j0=Max(n1*Coords[1] + Min(Coords[1],k1) - 1,0);
	NX=(i0+NX == N0)?(NX):(NX+1);
	NY=(j0+NY == N1)?(NY):(NY+1);


    if (rank == 0){
		printf("The Domain: [%f,%f]x[%f,%f], number of points: %dx%d;\n"
				   "The conjugate gradient iterations number: %d\n",
				    A0, A1 , B0, B1, N0,N1,CGMNum);
	}

	XNodes = (double *)malloc(N0*sizeof(double));
	YNodes = (double *)malloc(N1*sizeof(double));

	MeshGenerate(XNodes, YNodes, N0, N1);
    
    
#ifdef Print
    printf("My Rank in Grid_Comm is %d. My topological coords is (%d,%d). Domain size is %d x %d nodes.\n"
           "My neighbours: left = %d, right = %d, down = %d, up = %d.\n"
           "My mesh boundaries are i=%d .. %d, j=%d .. %d\n",
           rank, Coords[0], Coords[1], NX, NY, left,right, down,up,i0,i0+NX-1,j0,j0+NY-1);
#endif

	time0 = MPI_Wtime();
    CGM(CGMNum, XNodes+i0, YNodes+j0, NX, NY, N0, N1, Grid_Comm);
	time1 = MPI_Wtime();

	if (rank == 0)
		printf("Time in seconds=%.6fs\t",time1-time0);

	free(XNodes); free(YNodes);
    MPI_Finalize();
    // The end of MPI session ...
    
    return 0;
}
