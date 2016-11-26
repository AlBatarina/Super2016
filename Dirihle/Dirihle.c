#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>


// Domain size.
const double A = 2.0;
const double B = 2.0;

int NX, NY;							        // the number of internal points on axes (ox) and (oy).
double * XNodes, * YNodes;						// mesh node coords are stored here.


#define Test
#define Print
#define TRUE  ((int) 1)
#define FALSE ((int) 0)
#define Step 10

#define Max(A,B) ((A)>(B)?(A):(B))
#define R2(x,y) ((x)*(x)+(y)*(y))
#define Cube(x) ((x)*(x)*(x))

#define hx(i)  (XNodes[i+1]-XNodes[i])
#define hy(j)  (YNodes[j+1]-YNodes[j])

#define LeftPart(P,i,j)\
((-(P[NX*(j)+i+1]-P[NX*(j)+i])/hx(i)+(P[NX*(j)+i]-P[NX*(j)+i-1])/hx(i-1))/(0.5*(hx(i)+hx(i-1)))+\
(-(P[NX*(j+1)+i]-P[NX*(j)+i])/hy(j)+(P[NX*(j)+i]-P[NX*(j-1)+i])/hy(j-1))/(0.5*(hy(j)+hy(j-1))))


#ifdef Test
    double Solution(double x,double y)
    {
        return 2.0/(1.0+R2(x,y));
    }
#endif

double BoundaryValue(double x, double y)
{
    #ifdef Test
        return Solution(x,y);
    #else
        // replace this default value by yours.
        return 0;
    #endif
}

int RightPart(double * rhs)
{
    int i, j;
    double kappa2 = (16.0/R2(A,B));

    #ifdef Test
        for(j=0; j<NY; j++)
            for(i=0; i<NX; i++)
                rhs[j*NX+i] = 8.0*(1.0-R2(XNodes[i],YNodes[j]))/Cube(1.0+R2(XNodes[i],YNodes[j]));
        return 0;
    #else
        memset(rhs,0,NX*NY*sizeof(double));
    
    // place your code here.
    
        return 0;
    #endif
}

int MeshGenerate(int NX, int NY)
{
	const double q = 1.5;
	int i;

	for(i=0; i<NX; i++)
		XNodes[i] = A*(pow(1.0+i/(NX-1.0),q)-1.0)/(pow(2.0,q)-1.0);
	for(i=0; i<NY; i++)
		YNodes[i] = B*(pow(1.0+i/(NY-1.0),q)-1.0)/(pow(2.0,q)-1.0);
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

int CGM(int NX, int NY, double A, double B){

	double * SolVect;					// the solution array.
	double * ResVect;					// the residual array.
	double * BasisVect;					// the vector of A-orthogonal system in CGM.
	double * RHS_Vect;					// the right hand side of Puasson equation.
	double sp, alpha, tau, NewValue, err;			// auxiliary values.
	int counter;						// the current iteration number.

	int i,j;
	char str[127];
	FILE * fp;

	printf("The Domain: [0,%f]x[0,%f], number of points: N[0,A] = %d, N[0,B] = %d;\n"
			   "The conjugate gradient iterations number: %d\n",
			    A,B, NX,NY,CGMNum);

	XNodes = (double *)malloc(NX*sizeof(double));
	YNodes = (double *)malloc(NY*sizeof(double));
	SolVect   = (double *)malloc(NX*NY*sizeof(double));
	ResVect   = (double *)malloc(NX*NY*sizeof(double));
	RHS_Vect  = (double *)malloc(NX*NY*sizeof(double));

// Initialization of Arrays
	MeshGenerate(NX, NY);
	memset(ResVect,0,NX*NY*sizeof(double));
	RightPart(RHS_Vect);
    
	for(i=0; i<NX; i++)
	{
		SolVect[i] = BoundaryValue(XNodes[i],0.0);
		SolVect[NX*(NY-1)+i] = BoundaryValue(XNodes[i],B);
	}
	for(j=0; j<NY; j++)
	{
		SolVect[NX*j] = BoundaryValue(0.0,YNodes[j]);
		SolVect[NX*j+(NX-1)] = BoundaryValue(A,YNodes[j]);
	}

// Iterations ...
	#ifdef Test
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				err = Max(err, fabs(Solution(XNodes[i],YNodes[j])-SolVect[NX*j+i]));
		fprintf(fp,"\nNo iterations have been performed. The residual error is %.12f\n", err);
	#endif

		BasisVect = ResVect;    // g(0) = r(k-1).
	ResVect = (double *)malloc(NX*NY*sizeof(double));
	memset(ResVect,0,NX*NY*sizeof(double));

// CGM iterations begin ...
// sp == (Ar(k-1),r(k-1)) == (Ag(0),g(0)), k=1.
	#ifdef Print
		printf("\nCGM iterations begin ...\n");
	#endif

	for(counter=1; counter<=CGMNum; counter++)
	{
	// The residual vector r(k) is calculating ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				ResVect[NX*j+i] = LeftPart(SolVect,i,j)-RHS_Vect[NX*j+i];

	// The value of product (Ar(k),g(k-1)) is calculating ...
		alpha = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				alpha += LeftPart(ResVect,i,j)*BasisVect[NX*j+i]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		alpha = alpha/sp;

	// The new basis vector g(k) is being calculated ...
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				BasisVect[NX*j+i] = ResVect[NX*j+i]-alpha*BasisVect[NX*j+i];

	// The value of product (r(k),g(k)) is being calculated ...
		tau = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				tau += ResVect[NX*j+i]*BasisVect[NX*j+i]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		
	// The value of product sp = (Ag(k),g(k)) is being calculated ...
		sp = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
				sp += LeftPart(BasisVect,i,j)*BasisVect[NX*j+i]*(0.5*(hx(i)+hx(i-1)))*(0.5*(hy(j)+hy(j-1)));
		tau = tau/sp;

	// The x(k+1) is being calculated ...
		err = 0.0;
		for(j=1; j < NY-1; j++)
			for(i=1; i < NX-1; i++)
			{
				NewValue = SolVect[NX*j+i]-tau*BasisVect[NX*j+i];
				err = Max(err, fabs(NewValue-SolVect[NX*j+i]));
				SolVect[NX*j+i] = NewValue;
			}

		if(counter%Step == 0)
		{
			printf("The %d iteration of CGM method has been carried out.\n", counter);
            
#ifdef Print            
            fprintf(fp,"\nThe iteration %d of conjugate gradient method has been finished.\n"
                       "The value of \\alpha(k) = %f, \\tau(k) = %f. The difference value is %f.\n",\
                        counter, alpha, tau, err);
#endif

#ifdef Test
			err = 0.0;
			for(j=1; j < NY-1; j++)
				for(i=1; i < NX-1; i++)
					err = Max(err, fabs(Solution(XNodes[i],YNodes[j])-SolVect[NX*j+i]));
			fprintf(fp,"The %d iteration of CGM have been performed. The residual error is %.12f\n",\
                        counter, err);
#endif
		}
	}
// the end of CGM iterations.

// printing some results ...
	printf("\nThe %d iterations are carried out. The error of iterations is estimated by %.12f.\n",
                CGMNum, err);

	sprintf(str,"PuassonSerial_ECGM_%dx%d.dat", NX, NY);
	fp = fopen(str,"w");
		fprintf(fp,"# This is the conjugate gradient method for descrete Puasson equation.\n"
				"# A = %f, B = %f, N[0,A] = %d, N[0,B] = %d, SDINum = %d, CGMNum = %d.\n"
				"# One can draw it by gnuplot by the command: splot 'MyPath\\FileName.dat' with lines\n",\
				A, B, NX, NY, SDINum, CGMNum);
		for (j=0; j < NY; j++)
		{
			for (i=0; i < NX; i++)
				fprintf(fp,"\n%f %f %f", XNodes[i], YNodes[j], SolVect[NX*j+i]);
			fprintf(fp,"\n");
		}
	fclose(fp);

	free(XNodes); free(YNodes);
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
    int left, right, up, down;      // the neighbours of the process.

	int CGMNum;						// the number of CGM iterations.

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
    
    if((power = IsPower(ProcNum)) < 0)
    {
        if(rank == 0)
            printf("The number of procs must be a power of 2.\n");
        return(3);
    }
    
    // MPI Library is being activated ...
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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
    MPI_Comm_rank(Grid_Comm, &rank);
    MPI_Cart_coords(Grid_Comm, rank, ndims, Coords);
    
    if(Coords[0] < k0)
        ++n0;
    if(Coords[1] < k1)
        ++n1;
    
    MPI_Cart_shift(Grid_Comm, 0, 1, &left, &right);
    MPI_Cart_shift(Grid_Comm, 1, 1, &down, &up);
    
#ifdef Print
    printf("My Rank in Grid_Comm is %d. My topological coords is (%d,%d). Domain size is %d x %d nodes.\n"
           "My neighbours: left = %d, right = %d, down = %d, up = %d.\n",
           rank, Coords[0], Coords[1], n0, n1, left,right, down,up);
#endif
    
    CGM(n0,n1,A,B,);

    MPI_Finalize();
    // The end of MPI session ...
    
    return 0;
}
