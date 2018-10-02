#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define PI 3.14159265358979323846

/*Variables*/
char parfile[100];
double **psi, **psi0, **u, **u0, **mu, **mu0;
double **comp;
double cleqm, cseqm, mueqm;
double pc, ml, ms, clinitial;
double rt, beta, lambda;
double **gradw2x, **gradw2y;
double **gradpsix, **gradpsiy;
double **a, **theta;
double **w, **tau, **wderivative;
double tau0, w0;
double K, D;
double alpha, gama, ue, delta;
double epsi4, ja, theta0;
double amplitude;
double dx, dy, dt;
double idx, idy, idt, iD;
int nx, ny, maxiter, iter = 0;
int nxcg, nycg;
int nxhalf, nyhalf;
int nxcghalf, nycghalf;
double ri;

//coarse grid
double **s, **m, **h, **h0, **ucg;
double lb, c;
int n, mwalker;
double hc;
int countwalkerupdates=1;

//linked list
struct list {
	struct list *next;
	double x, y;
        int nt;
};
struct list **backbone;

//control walker pointers
struct list *walker, *bufferwalker, *poswalker;

int imc=0; //choose to solve using mc algo or not
int itemp=0; //choose to solve for thermal field
int icomp=1; //choose to solve for solutal field

/*Functions*/
void setfiles(int args, char *argv[]); //check for proper input parameters from command line
double **create2Ddouble(int nx, int ny); //create 2D array of type double
void readinput(void); //read from input file
void initialize(void); //initialize and allocate memory for system wide variables
void stepahead(void); //step ahead in time
void updatecg(void); //update coarse grid
void updatewalker(int choice, double x, double y, int k, int l); //create or delete walkers
struct list **createpointerarray(void);
void initializewalker(); //initialize positions for walkers
void updatewalkerpos(void); //update positions for walkers
void updatewalkerlist(void); //update walker list depending on the position of the coarse grid
void freewalkermemory(void);
void updatecgu(void); //update u field on coarse grid
void updatecgh(void); //update h field on coarse grid
void updatecgm(void); //update m field on coarse grid
void createbackbone(void); //create backbone and lists for updation
void boundary(void); //set boundary conditions
void boundarymc(void); //set boundary conditions for mc simulation

// deterministic equation solvers
void computepsi(void); //compute on whole domain
void computeu(void);
void computecomp(void);
void applysymmetry(void); //cubic symmetry

// stochastic equation solvers
void computemc(void); //compute on fine grid only
void applysymmetrymc(void); //cubic symmetry in mc

void printoutput(int r); //print output into text file
void rawoutput(int r); //print output into binary file
void freememory(void); //free allocated memory
