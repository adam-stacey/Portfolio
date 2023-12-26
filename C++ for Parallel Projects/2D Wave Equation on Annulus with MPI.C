//==============================================================
//
// MPI VERSION OF PROJECT CODE
//
// SOLVE WAVE EQUATION ON ANNULUS WITH VARIABLE COEFFICIENT
//
// u_tt = del * (c(theta) del u) + F
//
//
//
//===============================================================

// Import libraries and define type
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <mpi.h>
using std::string;
using std::max;

#include <ctime>

#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN


typedef double real; // can be switched between double and float


// getCPU() : Return the current wall-clock time in seconds
#include "getCPU.h"

// include commands tp parse command line arguments
#include "parseCommand.h"

// include command to create local grids
#include "getLocalIndexBounds.h"

int getDimX(int np)
{
  
  if(np == 1) return 1;

  if(np == 2) return 1;

  if (np == 4) return 2;

  if(np == 8) return 2;

  if(np == 16) return 4;

  if(np == 32) return 4;

  if(np == 64) return 8;

  if(np == 128) return 8;

  if(np == 256) return 16;

  if(np == 512) return 16;

  if(np == 1024) return 32;

  return 0;

}

int main(int argc, char *argv[])
{

    // initialize MPI and get processor and rank info
    MPI_Init(&argc,&argv);
    int myRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    int np; // assume an even number of processors (or 1)
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    // check to see if np is even or 1
    if(np % 2 != 0 && np != 1)
    {
      if(myRank == 0)
        printf("Number of Processors must be even, aborting\n");

      abort();
    }

    // get number of processors in each direction

    int dimX = getDimX(np);
    int dimY = np/dimX;

    int pi2 = myRank/dimX; // processor index 2
    int pi1 = myRank - dimX*(pi2); // processor index 1

    // define pi
    const real pi = 4.*atan2(1.,1.);

    // describe what boundary conditions are possible 
    enum BoundaryConditionsEnum
    {
        periodic= -1,
        dirichlet= 1,
        neumann= 2,
        parallelGhost = 3
    };

    // number of spatial dimensions
    const int numberOfDimensions=2;

    int debug = 0;
    real ra = 1.0, rb = 1.5; // radius of annulus is 1/2 at inner ring, 1 at outer ring
    real wa = -pi, wb = pi; // theta varies from -pi to pi

    real tFinal = 0.1;
    int nr = 100, nw = nr; // number of grid points in each dimesion

    int writeSolution = 0; // option to write final solution to csv: 0 for no, 1 for yes
    int writeGrids = 0; // option to write X and Y grids to csv's
    string solutionFileName = "serialSolution.csv";

    // parse commands from command line
    string line;
    bool echo; // turn echo on/off
    if(np == 1)
      echo = true;
    else
      echo = false;
    
    for( int i=1; i<argc; i++ )
    {
        line=argv[i];
        if( parseCommand( line,"-nr=",nr,echo) )
        {
            nw=nr;
        }
        else if( parseCommand( line,"-debug=",debug, echo) ){}
        else if( parseCommand( line, "-tFinal=",tFinal, echo) ){}
        else if( parseCommand( line, "-writeSolution=",writeSolution, echo) ) {}
        else if( parseCommand( line, "-writeGrids=",writeGrids, echo) ) {}
        else if( parseCommand( line, "-solutionFileName=",solutionFileName, echo) ) {}
    }


    FILE *debugFile=NULL;
    if(debug)
    {
        // open a debug file on each processor (in the debug folder)
        char debugFileName[80];
        sprintf(debugFileName,"debug/debugFileProjectMPIProc%dOf%d.debug",myRank,np);
        debugFile = fopen(debugFileName,"w");
        fprintf(debugFile,"---------- Project MPI Code--------------\n");
        fprintf(debugFile," np=%d, myRank=%d\n",np,myRank);
    }

     // define solution macro, if not specified at compile time
    // to specify at compile time, do something like g++ -DSOLUTIONTYPE=0 ...
    #ifndef SOLUTIONTYPE
        #define SOLUTIONTYPE 0 // default is manufactured solution
    #endif

    // define global grid

    // dim 1 is radius
    const int numGhost = 1;
    const int n1a = 0;
    const int n1b = n1a + nr;
    const int nd1a = n1a - numGhost;
    const int nd1b = n1b + numGhost;
    const int nd1 = nd1b - nd1a + 1;

    const int n2a = 0;
    const int n2b = n2a + nw;
    const int nd2a = n2a - numGhost;
    const int nd2b = n2b + numGhost;
    const int nd2 = nd2b - nd2a + 1;

    // get local grid  
    int n1a_l;
    int n1b_l;
    int n2a_l;
    int n2b_l;
    int nd1_l;
    int nd2_l;
    getLocalIndexBounds (pi2, dimY, nw, nd2_l, n2a_l,n2b_l);
    getLocalIndexBounds (pi1, dimX, nr, nd1_l, n1a_l, n1b_l);
    
    int numParallelGhost = numGhost;
    int nd1a_l = n1a_l - numParallelGhost;
    int nd1b_l = n1b_l + numParallelGhost;
    nd1_l = nd1b_l - nd1a_l + 1;
    int nd2a_l = n2a_l - numParallelGhost;
    int nd2b_l = n2b_l + numParallelGhost;
    nd2_l = nd2b_l - nd2a_l + 1;
    
    if(debug)
    {
      fprintf(debugFile,"Global indicies:  nd1a: %d, n1a: %d,  n1b: %d, nd1b: %d\n",nd1a, n1a, n1b, nd1b);
      fprintf(debugFile,"Global ghosts:  nd2a:%d, n2a: %d  n2b: %d, nd2b: %d\n",nd2a, n2a, n2b, nd2b);
      fprintf(debugFile,"Local indices: nd1a_l: %d, n1a_l: %d n1b_l:%d, nd1b_l: %d\n",nd1a_l,n1a_l,n1b_l,nd1b_l);
      fprintf(debugFile,"Local indices: nd2a_l: %d, n2a_l: %d n2b_l:%d, nd2b_l: %d\n",nd2a_l,n2a_l,n2b_l,nd2b_l);
    }
    

    // set up for boundary conditions and solution type   
    int *boundaryCondition_p = new int [2*numberOfDimensions];
    #define boundaryCondition(side,axis) boundaryCondition_p[(side)+2*(axis)]

    // define solution, variable coefficient, forcing term, and boundary conditions
    const real a0 = 0.1, a1 = 0.5, a2 = 0.1;
    const real b0 = 0.25, b1 = 0.4, b2 = 0.1;
    const real c0 = 0.3, c1 = 0.6, c2 = 0.2;
    const real p1 = 1.0;

    #if SOLUTIONTYPE == 0
        #define UTRUE(r,w,t)  ( (a0 + (t)*(a1 + a2 *(t) ) ) * (b0 + (r)*(b1 + b2*(r) ) ) * (c0 + (w)*(c1 + c2*(w) ) ) ) 
        #define UTRUE_R(r,w,t) ( (a0 + (t)*(a1 + a2 *(t) ) ) * ( b1 + 2.*b2*(r) ) * (c0 + (w)*(c1 + c2*(w) ) ) )
        #define UTRUE_RR(r,w,t) ( ( (a0 + (t)*(a1 + a2 *(t) ) ) * ( 2.*b2 ) * (c0 + (w)*(c1 + c2*(w) ) ) ) )
        #define UTRUE_T(r,w,t) ( (a1 + 2.*a2*(t) ) * (b0 + (r)*(b1 + b2*(r) ) ) * (c0 + (w)*(c1 + c2*(w) ) ) ) 
        #define UTRUE_TT(r,w,t) ( (2.*a2) * (b0 + (r)*(b1 + b2*(r) ) ) * (c0 + (w)*(c1 + c2*(w) ) ) )
        #define UTRUE_W(r,w,t) ( (a0 + (t)*(a1 + a2 *(t) ) ) * (b0 + (r)*(b1 + b2*(r) ) ) * (c1 + 2.*c2*(w) ) )
        #define UTRUE_WW(r,w,t) ( (a0 + (t)*(a1 + a2 *(t) ) ) * (b0 + (r)*(b1 + b2*(r) ) ) * (2.*c2) )
        #define LAPLACE_R(r,w,t) ( COEF(w)*UTRUE_RR(r,w,t) + COEF(w)*(p1/(r))*UTRUE_R(r,w,t) ) // components of laplacian wrt radius
        #define LAPLACE_W(r,w,t) ( p1/((r)*(r)) ) * ( COEFP(w)*UTRUE_W(r,w,t) + COEF(w)*UTRUE_WW(r,w,t) ) // component of laplacian wrt to theta
        #define FORCE(r,w,t) ( UTRUE_TT(r,w,t) - LAPLACE_R(r,w,t) - LAPLACE_W(r,w,t) )  // forcing term for manufactured solution

        // set solution name
        const char solutionName[] = "manufactured";

	// set boundary conditions
	boundaryCondition(0,0) = n1a_l == 0 ? dirichlet : parallelGhost; // left
        boundaryCondition(1,0) = n1b_l == nr ? dirichlet : parallelGhost; //right
        boundaryCondition(0,1) = n2a_l == 0 ? dirichlet : parallelGhost; // bottom
        boundaryCondition(1,1) = n2b_l == nw ? dirichlet : parallelGhost; // top

    #elif SOLUTIONTYPE == 1 
        #define UTRUE(r,w,t)  ( 0.0 )
        #define UTRUE_R(r,w,t) ( 0.0 )
        #define UTRUE_RR(r,w,t) ( 0.0 )
        #define UTRUE_T(r,w,t) ( 0.0 )
        #define UTRUE_TT(r,w,t) ( 0.0 )
        #define UTRUE_W(r,w,t) ( 0.0 )
        #define UTRUE_WW(r,w,t) ( 0.0)
     
        #define FORCE(r, w,t)					       \
        ({                                                             \
            real answer;                                               \
            real c = 0.5*(ra + rb);				       \
            if ( abs((r) - c) < dx[0] && abs((w) - 0.5*pi) < dx[1] )   \
            {							       \
                answer = 1.;                                           \
            } else {                                                   \
                answer = 0.;                                           \
            }                                                          \
            answer;                                                    \
        })

	// set solution name
	const char solutionName[] = "wave";

	// set boundary conditions
	boundaryCondition(0,0) = n1a_l == 0 ? dirichlet : parallelGhost; // left
        boundaryCondition(1,0) = n1b_l == nr ? dirichlet : parallelGhost; //right
        boundaryCondition(0,1) = n2a_l == 0 ? periodic : parallelGhost; // bottom
        boundaryCondition(1,1) = n2b_l == nw ? periodic : parallelGhost; // top
	
    #endif

    // define macro for variable coefficient
    #ifndef COEFTYPE
      #define COEFTYPE 0 // 0 for constant coef of 1.0, 1 for periodic variable coef
    #endif


    #if COEFTYPE == 1
        #define COEF(w) (sin(p1*(w))*sin(p1*(w)) + 1.0) // variable coef
        #define COEFP(w) (2.*sin(p1*(w))*cos(p1*(w)) ) // derivative of variable coef

	const char coefName[] = "variable";
	
    #elif COEFTYPE == 0
        #define COEF(w) (1.0) // variable coef
        #define COEFP(w) (0.0 ) // derivative of variable coef

	const char coefName[] = "constant";

    #endif

    if(debug)
      {
        fprintf(debugFile, "Solution Type: %s\n",solutionName);
        fprintf(debugFile, "Coef Type: %s\n",coefName);
        fprintf(debugFile, "Left BC: %d, Right BC: %d, Bottom BC: %d, Top BC: %d\n",
                boundaryCondition(0,0),boundaryCondition(1,0),boundaryCondition(0,1), boundaryCondition(1,1) );
       }


    // define spatial grid
    real dx[2];
    dx[0] = (rb-ra)/nr; // grid spacing in r direction
    dx[1] = (wb-wa)/nw; // grid spacing in w direction
    
    // define local grids
    real *r_l = new real [nd1_l*nd2_l];
    real *w_l = new real[nd1_l*nd2_l];
    #define R(i,j) r_l[((i) - nd1a_l) + nd1_l*((j) - nd2a_l)]
    #define W(i,j) w_l[((i) - nd1a_l) + nd1_l*((j) - nd2a_l)]

    for(int i = nd1a_l; i<= nd1b_l; i++)
    for(int j = nd2a_l; j <= nd2b_l; j++)
    {
	   R(i,j) =  ra + (i - n1a)*dx[0];
	   W(i,j) =  wa + (j - n2a)*dx[1];
    }
    

    // define time step

    // first, compute maximum value of coefficient
    
    #if COEFTYPE == 1
        real cMax = 2.0;
    #elif COEFTYPE == 0
        real cMax = 1.0;
    #endif
   
      
    // compute time step
    real cfl = 0.9;
    real dr2 = dx[0]*dx[0];
    real dw2 = dx[1]*dx[1];
    real dt = cfl * sqrt( 0.5 * ( 1.0/cMax) * (1.0/ ( 1.0/dr2 + 1.0/dw2 ) ) );
    int numSteps=ceil(tFinal/dt);
    dt = tFinal/numSteps; // adjust dt to reach the final time

    // solution storage
    
    real *u_p[3];
    u_p[0] = new real [nd1_l*nd2_l];
    u_p[1] = new real [nd1_l*nd2_l];
    u_p[2] = new real [nd1_l*nd2_l];

    #define um(i1,i2) u_p[past][((i1) - nd1a_l) + nd1_l*((i2) - nd2a_l)]
    #define uc(i1,i2) u_p[cur ][((i1) - nd1a_l) + nd1_l*((i2) - nd2a_l)]
    #define un(i1,i2) u_p[next][((i1) - nd1a_l) + nd1_l*((i2) - nd2a_l)]

    // apply initial conditions and copy to GPU
    
    int past = 0;
    int cur = 1;
    real dt2 = dt*dt;
    for(int i2 = nd2a_l; i2 <= nd2b_l; i2++)
    for(int i1 = nd1a_l; i1 <= nd1b_l; i1++)
    {
       real r = R(i1,i2), w = W(i1,i2);
       real C = COEF(w), Cp = COEFP(w);
       real overR = 1.0/r, overR2 = 1.0/(r*r);

       um(i1,i2) = UTRUE(r,w,0.0);
       uc(i1,i2) = um(i1,i2) + dt*UTRUE_T(r,w,0.0) + 0.5 * dt2 * ( C*UTRUE_RR(r,w,0.0) +
			   C*overR*UTRUE_R(r,w,0.0) + overR2 * ( Cp*UTRUE_W(r,w,0.0) +
			   C*UTRUE_WW(r,w,0.0) ) + FORCE(r,w,0.0) );
    }

    
    // print out preliminary information
    if (!myRank)
    {
      
      printf("-----------------Solve the Wave Equation on an Annulus--------------------------\n");
      printf(" Project MPI Code: Number of Processors: %d\n",np);
      printf(" Solution Type: %d, Solution Name: %s\n",SOLUTIONTYPE,solutionName);
      printf(" Coefficient Type: %d, Coefficient Name: %s\n",COEFTYPE,coefName);
      printf(" Radius: %6.2f to %6.2f, Theta: %6.4f to %6.4f, tFinal = %6.4f\n",ra,rb,wa,wb,tFinal);
      printf(" Nr=%d, Nw=%d, numSteps = %d, debug = %d\n",nr,nw,numSteps,debug);
      printf(" nd1a = %d, n1a = %d, n1b = %d, nd1b = %d, n1a = %d\n",nd1a,n1a,n1b,nd1b,nd1);
      printf(" nd2a = %d, n2a = %d, n2b = %d, nd2b = %d, n2a = %d\n",nd2a,n2a,n2b,nd2b,nd2);
    }
    
    //----------------------TIME STEP-----------------------------
    
    // define necessary constants
    const real overDr2 = 1.0/dr2;
    const real overDw2 = 1.0/dw2;
    const real overDr = 0.5/dx[0];
    const real overDw = 0.5/dx[1];

    // define time
    real t = dt;

    // record initial cpu time
    real cpu0 = getCPU();

    // Allocate Buffers to send/recv parallel ghost points
    real *lrBuffer = new real [nd2_l];
    real *tbBuffer = new real [nd1_l];
    
    
    for(int n = 2; n <= numSteps; n++)
      {

	// update arrays
        int past = (n + 1) % 3;
        int cur = (n + 2) % 3;
        int next = (n + 3) % 3;
	// update interior scheme
	for(int i2 = n2a_l; i2 <= n2b_l; i2++)
	  for(int i1 = n1a_l; i1 <= n1b_l; i1++)
	    {
	      
	      // compute coefficients
	      real r = R(i1,i2), w = W(i1,i2); // radius and theta for this point
	      real overR = 1.0/(r);
	      real overR2 = 1.0/(r*r); // 1/r and 1/r^2 for given point
	      real C = COEF(w), Cp = COEFP(w); // variable coef and its derivative evaluated at given point

	      // apply update
	      un(i1,i2) = 2.*uc(i1,i2) - um(i1,i2) + dt2 * ( C * overDr2 *( uc(i1 + 1,i2) - 2.*uc(i1,i2) + uc(i1 - 1,i2) )
						       + C * overDr * overR *( uc(i1 + 1,i2) - uc(i1 - 1,i2) )
						       + overR2 * ( Cp * overDw * ( uc(i1,i2 + 1) - uc(i1, i2 - 1) )
						       + C * overDw2 * ( uc(i1, i2 + 1) - 2.*uc(i1,i2) + uc(i1,i2 - 1) ) ) ) + dt2*FORCE(r,w,t);
	      
	    }
	
	// update time
	t = n*dt;
    
	// update boundary terms
	for(int axis = 0; axis < numberOfDimensions; axis++)
	for(int side = 0; side < 2; side++)
	  {
            
            int is = 1 - 2*side; // 1 on left/bottom, -1 on right/top

            // Apply Dirichlet Boundary Conditions
	    if(boundaryCondition(side,axis) == dirichlet)
	      {
	       
                // left and right
		if(axis == 0)
		  {
		    int i1 = side == 0 ? n1a_l : n1b_l; // index of boundary 
		    int i1g = side == 0 ? nd1a_l : nd1b_l; // index of ghost
		    for(int i2 = nd2a_l; i2 <= nd2b_l; i2++)
		       {
			 un(i1,i2) = UTRUE(R(i1,i2),W(i1,i2),t);
			 un(i1g,i2) = 3.*un(i1,i2) - 3.*un(i1+is,i2) + un(i1+2*is,i2); // extrap ghost
		       }

                    if(pi1 == 0 && dimX > 1 && !side)
                    {
                      for(int i2 = nd2a_l; i2 <= nd2b_l; i2++) {lrBuffer[i2 - nd2a_l] = un(n1b_l,i2);}
                      MPI_Send(lrBuffer,nd2_l,MPI_DOUBLE,myRank + 1,myRank,MPI_COMM_WORLD);
                      
                    }

                    if(pi1 == dimX - 1 && dimX > 1 && side)
                    {
                      for(int i2 = nd2a_l; i2 <= nd2b_l; i2++) {lrBuffer[i2 - nd2a_l] = un(n1a_l,i2);}
                      MPI_Send(lrBuffer,nd2_l,MPI_DOUBLE,myRank - 1,myRank,MPI_COMM_WORLD);
                    }
                    
		  }
                
		// bottom and top
		if(axis == 1)
		  {
		    int i2 = side == 0 ? n2a_l : n2b_l; // boundary
		    int i2g = side == 0? nd2a_l : nd2b_l; // ghost
		    for(int i1 = nd1a_l; i1 <= nd1b_l; i1++)
		       {
			  un(i1,i2) = UTRUE(R(i1,i2),W(i1,i2),t);
			  un(i1,i2g) = 3.*un(i1,i2) - 3.*un(i1,i2+is) + un(i1,i2+2*is); // extrap ghost

		       }

                    if(pi2 == 0 && dimY > 1 && !side)
                    {
                      for(int i1 = nd1a_l; i1 <= nd1b_l; i1++) {tbBuffer[i1 - nd1a_l] = un(i1,n2b_l);}
                      MPI_Send(tbBuffer,nd1_l,MPI_DOUBLE,myRank + dimX,myRank,MPI_COMM_WORLD);

                    }

                    if(pi2 == dimY - 1 && dimY > 1 && side)
                    {
                      for(int i1 = nd1a_l; i1 <= nd1b_l; i1++) {tbBuffer[i1 - nd1a_l] = un(i1,n2a_l);}
                      MPI_Send(tbBuffer,nd1_l,MPI_DOUBLE,myRank - dimX,myRank,MPI_COMM_WORLD);
                    }

		  }
                

	      } // end of dirichlet condition case

	    // periodic boundary conditions
	    if(boundaryCondition(side,axis) == periodic)
	      {
		// periodic conditons only applied to theta
		if(axis == 1 && side == 0)
		  {  
		    for(int i1 = nd1a_l; i1 <= nd1b_l; i1++)
		      {
			un(i1,n2b_l) = un(i1,n2a_l);
			un(i1,nd2b_l) = un(i1,n2a_l + 1);
			un(i1,nd2a_l) = un(i1,n2b_l - 1);
		      }
                  }
                
                    if(pi2 == 0 && dimY > 1 && !side)
                    {
                      for(int i1 = nd1a_l; i1 <= nd1b_l; i1++) {tbBuffer[i1 - nd1a_l] = un(i1,n2b_l);}
                      MPI_Send(tbBuffer,nd1_l,MPI_DOUBLE,myRank + dimX,myRank,MPI_COMM_WORLD);
                    }

                    if(pi2 == dimY - 1 && dimY > 1 && side)
                    {
                      for(int i1 = nd1a_l; i1 <= nd1b_l; i1++) {tbBuffer[i1 - nd1a_l] = un(i1,n2a_l);}
                      MPI_Send(tbBuffer,nd1_l,MPI_DOUBLE,myRank - dimX,myRank,MPI_COMM_WORLD);
                    }
              } // end of periodic boundary conditions
            
            
            //parallel ghost boundary points
            if(boundaryCondition(side,axis) == parallelGhost)
            {
              // left and right
              if(axis == 0)
              {
                    int recvDir = myRank - is;
                    int sendCond = boundaryCondition(side + is,axis) == parallelGhost ? 1 : 0;
                    int i1g = side == 0 ? nd1a_l : nd1b_l; // index of ghost
                    int leftIndex = nd2a_l;
                    int rightIndex = nd2b_l;
                    
                    MPI_Recv(lrBuffer,nd2_l,MPI_DOUBLE,recvDir,recvDir,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                    // save information from buffer
                   
                    for(int i2 = leftIndex; i2 <= rightIndex; i2++) {un(i1g,i2) = lrBuffer[i2 - leftIndex];}
                    
                    // send along data if necessary

                    if(sendCond)
                    {
                      int sendDir = myRank + is;
                      int tag = myRank;
                      int sendLeftRight = side == 0 ? n1b_l : n1a_l;
                      for(int i2 = leftIndex; i2 <= rightIndex; i2++){ lrBuffer[i2 - leftIndex] = un(sendLeftRight,i2);}
                      MPI_Send(lrBuffer,nd2_l,MPI_DOUBLE,sendDir,tag,MPI_COMM_WORLD);  
                    }
                    

              } // end of left/right parallel ghosts

              // top/bottom
              if(axis == 1)
              {
                    int recvDir = myRank - is*dimX;
                    int sendCond = boundaryCondition(side + is,axis) == parallelGhost ? 1 : 0;
                    int i2g = side == 0 ? nd2a_l : nd2b_l; // index of ghost
                    int leftIndex = nd1a_l;
                    int rightIndex = nd1b_l;

                    MPI_Recv(tbBuffer,nd1_l,MPI_DOUBLE,recvDir,recvDir,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                    // save information from buffer

                    for(int i1 = leftIndex; i1 <= rightIndex; i1++) {un(i1,i2g) = tbBuffer[i1 - leftIndex];}

                    // send along data if necessary

                    if(sendCond)
                    {
                      int sendDir = myRank + is*dimX;
                      int tag = myRank;
                      int sendTopBottom = side == 0 ? n2b_l : n2a_l;
                      for(int i1 = leftIndex; i1 <= rightIndex; i1++){ tbBuffer[i1 - leftIndex] = un(i1,sendTopBottom);}
                      MPI_Send(tbBuffer,nd1_l,MPI_DOUBLE,sendDir,tag,MPI_COMM_WORLD);
                    }
                    

              } // end of top/bottom parallel ghosts
              


            } // end of parallel ghost conditions
            

              

	  } // end of boundary condition loop


	// check time step error
	if(debug)
	{
		real maxErr = 0.0;
		for(int j = n2a_l; j <= n2b_l; j++)
		for(int i = n1a_l; i <= n1b_l; i++)
		{
		  real err = abs(un(i,j) - UTRUE(R(i,j),W(i,j),t));
		  maxErr = max(err,maxErr);
		}
		fprintf(debugFile, "Step: %d, time = %f, maxErr = %.2e\n",n,t,maxErr);
	}
        
       

      } // end of timestep
    

    // record final cpu time and get global time
    real cpuTime = getCPU() - cpu0;
    real globalTime;
    MPI_Allreduce(&cpuTime,&globalTime,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

    // collect global solution on rank 0

    // define global solution array
    real *globalU_p = new real [nd1*nd2];
    #define globalSolution(i,j) globalU_p[((i) - nd1a) + nd1*((j) - nd2a)]

    // collect local solution on each processor
    real *localSolution = new real[nd1_l*nd2_l];
    int next = (numSteps + 3) % 3;
    for(int i2 = nd2a_l; i2 <= nd2b_l; i2++)
      for(int i1 = nd1a_l; i1 <= nd1b_l; i1++)
        localSolution[((i1) - nd1a_l) + nd1_l*((i2) - nd2a_l)] = un(i1,i2);
    

  
    // save solution from rank 0
    if(!myRank)
    {
      for(int j = nd2a_l; j <= nd2b_l;j++)
      for(int i = nd1a_l; i <= nd1b_l; i++)
          globalSolution(i,j) = localSolution[((i) - nd1a_l) + nd1_l*((j) - nd2a_l)];
      
    }

    if(np > 1)
    {
      for(int rank = 1; rank < np; rank++)
      {
        
        int pi2_s = rank/dimX; // processor index 2
        int pi1_s = rank - dimX*(pi2_s); // processor index 1
        
        // get local indices of sending rank
        getLocalIndexBounds (pi2_s, dimY, nw, nd2_l, n2a_l,n2b_l);
        getLocalIndexBounds (pi1_s, dimX, nr, nd1_l, n1a_l,n1b_l);
        int nd1a_l = n1a_l - numGhost;
        int nd2a_l = n2a_l - numGhost;
        int nd1b_l = n1b_l + numGhost;
        int nd2b_l = n2b_l + numGhost;
        int nd1_l = nd1b_l - nd1a_l + 1;
        int nd2_l = nd2b_l - nd2a_l + 1;
        
        
        int sendSize = nd1_l*nd2_l;
        
        
        int topIndex = n2b_l;
        int bottomIndex = n2a_l;
        int rightIndex = n1b_l;
        int leftIndex = n1a_l;

        if(pi2_s == dimY - 1) {++topIndex;}
        if(pi2_s == 0) {--bottomIndex;}
        if(pi1_s == dimX - 1) {++rightIndex;}
        if(pi1_s == 0) {--leftIndex;} 

        real *temp = new real [sendSize];
        
        
        if(myRank == rank)
        {
          MPI_Send(localSolution,sendSize,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
        }

        if(myRank == 0)
        {
          MPI_Recv(temp,sendSize,MPI_DOUBLE,rank,rank,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

          for(int i2 = bottomIndex; i2 <= topIndex; i2++)
            for(int i1 = leftIndex; i1 <= rightIndex; i1++)
            {
              int index = (i1 - nd1a_l) + nd1_l*(i2 - nd2a_l);
              
              globalSolution(i1,i2) = temp[index ];

            }
        }
        
        
     
        delete [] temp;
        
      }// end of rank loop
    } // end of if (np > 1)

    // compute global error
    real globalError = 0.;
    real globalMaxNorm = 0.;
    if(!myRank)
    {
      
      for(int i2 = n2a; i2 <= n2b; i2++)
        for(int i1 = n1a; i1 <= n1b; i1++)
        {
          real r = ra + (i1 - n1a)*dx[0];
          real w = wa + (i2 - n2a)*dx[1];
        
          real err = abs(globalSolution(i1,i2) - UTRUE(r,w,tFinal));
          globalError = max(err,globalError);
          globalMaxNorm = max(globalSolution(i1,i2),globalMaxNorm);
          
          
        }
    }
    
    globalError /= max(globalMaxNorm,REAL_MIN); // relative error

    if(!myRank)
    {
      printf("\n========================================================\n");
      printf("Maximum Relative Error = %9.2e, Max Norm = %9.2e, CPU(s) = %9.2e\n",globalError,globalMaxNorm,globalTime);
      printf("========================================================\n");
    }
    

    if(debug)
    {
      fflush(debugFile);
      fclose(debugFile);
    }

    
    // write solution and global grids if requested
    if(writeGrids && !myRank)
      {

	// open csv files for the X and Y grids
	FILE *gridX;
	FILE *gridY;

	gridX = fopen("X_Grid.csv","w+");
	gridY = fopen("Y_Grid.csv","w+");

	// write to X and Y files
	for(int i2 = nd2a; i2 <= nd2b; i2++)
	  {
	    for(int i1 = nd1a; i1 <= nd1b; i1++)
	      {
                real r = ra + (i1 - n1a)*dx[0];
                real w = wa + (i2 - n2a)*dx[1];
                
		fprintf(gridX,"%f,",r );
		fprintf(gridY,"%f,",w );
	      }
	    fprintf(gridX,"\n");
	    fprintf(gridY,"\n");
	  }

	// close files
	fclose(gridX);
	fclose(gridY);
	printf("Wrote grids!\n");
			
      }

    if(writeSolution && !myRank)
      {

	// open csv
	FILE *fpt;
	fpt = fopen(solutionFileName.c_str(),"w+");

	// write solution data
	for(int i2 = nd2a; i2 <= nd2b; i2++)
	  {
	    for(int i1 = nd1a; i1 <= nd1b; i1++)
	      {
		fprintf(fpt,"%f,",globalSolution(i1,i2));
	      }
	    fprintf(fpt,"\n");
	  }

	// close file
	fclose(fpt);
	printf("Wrote solution at final time t = %.2f\n",tFinal);

      }

    
    
    // clean up memory
    delete [] boundaryCondition_p;
    delete [] r_l;
    delete [] w_l;
    delete [] u_p[0];
    delete [] u_p[1];
    delete [] u_p[2];
    delete [] globalU_p;
    delete [] localSolution;
    delete [] lrBuffer;
    delete [] tbBuffer;
    

    // close down MPI
    MPI_Finalize();


    return 0;

}
