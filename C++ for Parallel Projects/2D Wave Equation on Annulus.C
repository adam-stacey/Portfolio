//==============================================================
//
// SERIAL VERSION OF PROJECT CODE
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
using std::string;
using std::max;

#include <ctime>

#define REAL_EPSILON DBL_EPSILON
#define REAL_MIN DBL_MIN


typedef double real; // can be switched between double and float


// getCPU() : Return the current wall-clock time in seconds
inline double getCPU()
{
    return ( 1.0*std::clock() )/CLOCKS_PER_SEC ;
}

// include commands tp parse command line arguments
#include "parseCommand.h"


int main(int argc, char *argv[])
{


    // define pi
    const real pi = 4.*atan2(1.,1.);

    // describe what boundary conditions are possible 
    enum BoundaryConditionsEnum
    {
        periodic= -1,
        dirichlet= 1,
        neumann= 2,
	ghost = 3
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
    bool echo = true; // turn echo on/off
    
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

     // define solution macro, if not specified at compile time
    // to specify at compile time, do something like g++ -DSOLUTIONTYPE=0 ...
    #ifndef SOLUTIONTYPE
        #define SOLUTIONTYPE 0 // default is manufactured solution
    #endif

    // define grid indices
    
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
	boundaryCondition(0,0) = dirichlet; // left
        boundaryCondition(1,0) = dirichlet; //right
        boundaryCondition(0,1) = dirichlet; // bottom
        boundaryCondition(1,1) = dirichlet; // top

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
	boundaryCondition(0,0) = dirichlet; // left
        boundaryCondition(1,0) = dirichlet; //right
        boundaryCondition(0,1) = periodic; // bottom
        boundaryCondition(1,1) = periodic; // top
	
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


    // define spatial grid
    real dx[2];
    dx[0] = (rb-ra)/nr; // grid spacing in r direction
    dx[1] = (wb-wa)/nw; // grid spacing in w direction

    real *r_p = new real [nd1*nd2];
    real *w_p = new real[nd1*nd2];
    #define R(i,j) r_p[((i) - nd1a) + nd1*((j) - nd2a)]
    #define W(i,j) w_p[((i) - nd1a) + nd1*((j) - nd2a)]

    for(int i = nd1a; i<= nd1b; i++)
    for(int j = nd2a; j <= nd2b; j++)
    {
	   R(i,j) =  ra + (i - n1a)*dx[0];
	   W(i,j) =  wa + (j - n2a)*dx[1];
    }


    if (debug > 1)
    {

	   printf("Check R Grid:\n");
	   for(int j = nd2b; j >= nd2a; j--)
	   for(int i = nd1a; i <= nd1b; i++)
	   {
		   printf("%.2f ",R(i,j));
		   if(i == nd1b)
			printf("\n");
	   }
	   printf("\nCheck W Grid:\n");
	   for(int j = nd2b; j >= nd2a; j--)
	   for(int i = nd1a; i <= nd1b; i++)
	   {
		   printf("%.2f ",W(i,j));
		   if(i == nd1b)
			printf("\n");
	   }
	   printf("\n");
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
    u_p[0] = new real [nd1*nd2];
    u_p[1] = new real [nd1*nd2];
    u_p[2] = new real [nd1*nd2];

    #define um(i1,i2) u_p[past][((i1) - nd1a) + nd1*((i2) - nd2a)]
    #define uc(i1,i2) u_p[cur ][((i1) - nd1a) + nd1*((i2) - nd2a)]
    #define un(i1,i2) u_p[next][((i1) - nd1a) + nd1*((i2) - nd2a)]

    // apply initial conditions and copy to GPU
    int past = 0;
    int cur = 1;
    real dt2 = dt*dt;
    for(int i = nd1a; i <= nd1b; i++)
    for(int j = nd2a; j <= nd2b; j++)
    {
       real r = R(i,j), w = W(i,j);
       real C = COEF(w), Cp = COEFP(w);
       real overR = 1.0/r, overR2 = 1.0/(r*r);

       um(i,j) = UTRUE(r,w,0.0);
       uc(i,j) = um(i,j) + dt*UTRUE_T(r,w,0.0) + 0.5 * dt2 * ( C*UTRUE_RR(r,w,0.0) +
			   C*overR*UTRUE_R(r,w,0.0) + overR2 * ( Cp*UTRUE_W(r,w,0.0) +
			   C*UTRUE_WW(r,w,0.0) ) + FORCE(r,w,0.0) );
    }

    if(debug > 1)
    {
	   printf("After Initial Conditions, Error in UM = :\n");
	   for(int j = nd2b; j >= nd2a; j--)
	   for(int i = nd1a; i <= nd1b; i++)
	   {
	      printf("%.14f ",um(i,j) - UTRUE(R(i,j),W(i,j),0.0));
	       if(i == nd1b)
		  printf("\n");
	    }
	    printf("\n");
	    printf("After Initial Conditions, Error in UC = :\n");
            for(int j = nd2b; j >= nd2a; j--)
            for(int i = nd1a; i <= nd1b; i++)
            {
	       printf("%.14f ",uc(i,j)- UTRUE(R(i,j),W(i,j),dt));
               if(i == nd1b)
                  printf("\n");
            }
            printf("\n");
    }

    // print out preliminary information
    printf("-----------------Solve the Wave Equation on an Annulus--------------------------\n");
    printf(" Project Serial Code: No Parallelization!\n");
    printf(" Solution Type: %d, Solution Name: %s\n",SOLUTIONTYPE,solutionName);
    printf(" Coefficient Type: %d, Coefficient Name: %s\n",COEFTYPE,coefName);
    printf(" Radius: %6.2f to %6.2f, Theta: %6.4f to %6.4f, tFinal = %6.4f\n",ra,rb,wa,wb,tFinal);
    printf(" Nr=%d, Nw=%d, numSteps = %d, debug = %d\n",nr,nw,numSteps,debug);
    printf(" nd1a = %d, n1a = %d, n1b = %d, nd1b = %d, n1a = %d\n",nd1a,n1a,n1b,nd1b,nd1);
    printf(" nd2a = %d, n2a = %d, n2b = %d, nd2b = %d, n2a = %d\n",nd2a,n2a,n2b,nd2b,nd2);
    printf(" Left BC: %d, Right BC: %d, Bottom BC: %d, Top BC: %d\n",boundaryCondition(0,0), boundaryCondition(0,1),
	   boundaryCondition(1,0), boundaryCondition(1,1) );
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
    for(int n = 2; n <= numSteps; n++)
      {

	// update arrays
        int past = (n + 1) % 3;
        int cur = (n + 2) % 3;
        int next = (n + 3) % 3;
	// update interior scheme
	for(int i2 = n2a; i2 <= n2b; i2++)
	  for(int i1 = n1a; i1 <= n1b; i1++)
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
	    // Apply Dirichlet Boundary Conditions
	    if(boundaryCondition(side,axis) == dirichlet)
	      {
		int is = 1 - 2*side; // 1 on left/bottom, -1 on right/top
		// left and right
		if(axis == 0)
		  {
		    int i1 = side == 0 ? n1a : n1b; // index of boundary 
		    int i1g = side == 0 ? nd1a : nd1b; // index of ghost
		    for(int i2 = nd2a; i2 <= nd2b; i2++)
		       {
			 un(i1,i2) = UTRUE(R(i1,i2),W(i1,i2),t);
			 un(i1g,i2) = 3.*un(i1,i2) - 3.*un(i1+is,i2) + un(i1+2*is,i2); // extrap ghost
		       }

		  }
		// bottom and top
		if(axis == 1)
		  {
		    int i2 = side == 0 ? n2a : n2b; // boundary
		    int i2g = side == 0? nd2a : nd2b; // ghost
		    for(int i1 = nd1a; i1 <= nd1b; i1++)
		       {
			  un(i1,i2) = UTRUE(R(i1,i2),W(i1,i2),t);
			  un(i1,i2g) = 3.*un(i1,i2) - 3.*un(i1,i2+is) + un(i1,i2+2*is); // extrap ghost

		       }

		  }

	      } // end of dirichlet condition case

	    // periodic boundary conditions
	    if(boundaryCondition(side,axis) == periodic)
	      {
		// periodic conditons only applied to theta
		if(axis == 1 && side == 1)
		  {
		    for(int i1 = nd1a; i1 <= nd1b; i1++)
		      {
			un(i1,n2b) = un(i1,n2a);
			un(i1,nd2b) = un(i1,n2a + 1);
			un(i1,nd2a) = un(i1,n2b - 1);
		      }
		  }
	      } // end of periodic condtitions

	  } // end of boundary condition loop


	// check time step error
	if(debug)
	{
		real maxErr = 0.0;
		for(int j = n2a; j <= n2b; j++)
		for(int i = n1a; i <= n1b; i++)
		{
		  real err = abs(un(i,j) - UTRUE(R(i,j),W(i,j),t));
		  maxErr = max(err,maxErr);
		}
		printf("Step: %d, time = %f, maxErr = %.2e\n",n,t,maxErr);
	}
       

      } // end of timestep

    // record final cpu time
    real cpuTime = getCPU() - cpu0;

    // compute final error
    int next = (numSteps + 3) % 3;
    real cpuMaxErr = 0.;
    real cpuMaxNorm = 0.;
    for(int j = n2a; j <= n2b; j++)
    for(int i = n1a; i <= n1b; i++)
    {
	    real err = abs(un(i,j) - UTRUE(R(i,j),W(i,j),tFinal));
	    cpuMaxErr = max(err,cpuMaxErr);
	    cpuMaxNorm = max(cpuMaxNorm,un(i,j));

    }

    cpuMaxErr /= max(cpuMaxNorm,REAL_MIN); // relative error

    printf("\n========================================================\n");
    printf("Maximum Relative Error = %9.2e, Max Norm = %9.2e, CPU(s) = %9.2e\n",cpuMaxErr,cpuMaxNorm,cpuTime);
    printf("========================================================\n");


    // convert R and W to X and Y and prepare to export
    
    if(writeGrids)
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
		fprintf(gridX,"%f,",R(i1,i2)*cos(W(i1,i2) ) );
		fprintf(gridY,"%f,",R(i1,i2)*sin(W(i1,i2) ) );
	      }
	    fprintf(gridX,"\n");
	    fprintf(gridY,"\n");
	  }

	// close files
	fclose(gridX);
	fclose(gridY);
	printf("Wrote grids!\n");
			
      }

    if(writeSolution)
      {

	// open csv
	FILE *fpt;
	fpt = fopen(solutionFileName.c_str(),"w+");

	// write solution data
	for(int i2 = nd2a; i2 <= nd2b; i2++)
	  {
	    for(int i1 = nd1a; i1 <= nd1b; i1++)
	      {
		fprintf(fpt,"%f,",un(i1,i2));
	      }
	    fprintf(fpt,"\n");
	  }

	// close file
	fclose(fpt);
	printf("Wrote solution at final time t = %.2f\n",tFinal);

      }
      

    
       
    
    // clean up memory
    delete [] boundaryCondition_p;
    delete [] r_p;
    delete [] w_p;
    delete [] u_p[0];
    delete [] u_p[1];
    delete [] u_p[2];
    
    return 0;

}
