#include "pf.h"

int main(int args, char *argv[])
{
	int r;
	setfiles(args, argv);
	readinput();
	
        //Parameter calculation
        hc = n*n*delta/mwalker;
        nxcg = nx/n;
        nycg = ny/n;
        nxhalf = nx/2;
        nyhalf = ny/2;
        nxcghalf = nxcg/2;
        nycghalf = nycg/2;
        idx = 1.0/dx;
        idy = 1.0/dy;
        idt = 1.0/dt;
        iD = 1.0/D;
        rt = 6.0;
        
        pc = cseqm/cleqm;
        mueqm = rt*log(cleqm);
        /*mueqm = cleqm;*/
        beta = -rt*log(pc);
        
   printf("\n Maximum number of iterations = %d", maxiter);
   printf("\n Grid size = %d x %d", nx, ny);
   printf("\n Spacial and temporal steps = %lf %lf, %lf", dx, dy, dt);
   printf("\n Intial radius of solid = %lf", ri);
   
   printf("\n \n Phase field:");
   printf("\n Anisotropy parameter = %lf", epsi4);
   printf("\n Anisotropy mode number = %lf, Ref. angle = %lf", ja, theta0);	
   printf("\n Amplitude of disturbance = %lf", amplitude);
   printf("\n Surface parameter = %lf, Characteristic time = %lf", w0, tau0);
   
   printf("\n \n Thermal field:");
   printf("\n Dimensionless latent heat = %lf, Thermal Diffusivity = %lf",K, D);
   printf("\n Equilibrium temperature = %lf, undercooling = %lf", ue, delta);
   printf("\n Coupling parameters for thermal field: alpha = %lf, gamma = %lf", alpha, gama);
   
   printf("\n \n Solutal field:");
   printf("\n Equilibrium: Liq composition = %lf, Solid composition = %lf, Chemical potential = %lf", cleqm, cseqm, mueqm);
   printf("\n Partition coeff = %lf, Diffusion coeff in liquid = %lf, Diffusion coeff in solid = %lf ", pc, ml, ms);
   printf("\n Initial concentration of liquid = %lf", clinitial);
   printf("\n Coupling parameters for solutal field = %lf", lambda);
   
   printf("\n \n Monte Carlo:");
   printf("\n Random walker distance coupling parameter (c) = %lf, L_buffer = %lf, Critical H = %lf", c, lb, hc);
   printf("\n Coarse grid: spacial step size = %d, Initial no. of walkers = %d", n, mwalker);
   
   initialize();
   
   if(imc==1)
   {
        //allocate memory for global pointers
        walker = (struct list *) malloc( sizeof(struct list) );
        bufferwalker = (struct list *) malloc( sizeof(struct list) );
        poswalker = (struct list *) malloc( sizeof(struct list) );
        
	updatecg();
        initializewalker(); 
        createbackbone(); //calculate m here itself
   }
	printoutput(0);	
	for(r=0; r<maxiter; r++)
	{	
		printf("\n Iteration: %d", r+1);
		stepahead();

                if(imc==1)
                {
                    updatecg();
                    updatewalkerlist();
                    updatecgu();
                    updatecgh();
                    computemc();
                    applysymmetrymc();
                    boundary();
                    boundarymc();
                    updatewalkerpos(); //update m with it
                    freewalkermemory();
                    printf("\n Number of times walkers were updated = %d",countwalkerupdates);
                }
                else
                {
                    computepsi();
                    if(itemp==1)
                        computeu();
                    if(icomp==1)
                        computecomp();
                    applysymmetry();
                    boundary();
                }
                printf("\n");
		if((r+1)%1000==0)
			printoutput(r+1);		
	}
	
	freememory();
	printf("\n All done \n");
	return 0;
}

void setfiles(int args, char *argv[])
{
   int i;

   if(args < 3) {
    printf("Usage: %s -b parameter_filename\n",argv[0]);
    exit(0);
   } 
   else {
    for(i=1;i<args;i=i+2) 
    {
      if(strcmp(argv[i],"-b") == 0) 
      {
	strncpy(parfile,argv[i+1], (size_t) 64);
      } else 
      {
        printf("Unknown option : %s \n", argv[i]);
        exit(0);
      }
    }
   }
      
  return;
}

void restore(void)
{
/*
    printf("\n Reading from an initial configuration");
    FILE *fptr;
    char initfile[100];
    
    strcpy(initfile,"output20000.dat");
    printf("\n Opening intial configuration file %s", initfile);
    
    if((fptr == fopen(initfile,"r")) == NULL)
        printf("Error! opening file");
    
*/   
    return;
}

double **create2Ddouble(int nx, int ny)
{

  double *space;
  double **arr2d;
  int i, j;

  space = (double *) malloc (nx * ny * sizeof (double));

  arr2d = (double **) malloc (nx * sizeof (double *));
  for (i = 0; i < nx; i++) {
      arr2d[i] = space + i * ny;
    }

  for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	  arr2d[i][j] = 0.0;
	}
    }
  return (arr2d);
}

struct list **createpointerarray(void)
{
    backbone = (struct list **)malloc((maxiter+1)*sizeof(struct list *));
    
    int i;
    for(i=0;i<maxiter;i++)
    {
        backbone[i] = (struct list *)malloc(sizeof(struct list));
        backbone[i]->next = (struct list *)malloc(sizeof(struct list));
        backbone[i]->next = NULL;
        backbone[i]->x = 0.0;
        backbone[i]->y = 0.0;
        backbone[i]->nt = 0;
    }
    return backbone;
}

void readinput(void)
{ 
   printf("\n Reading input");
   int somenum;
   FILE *fptr;
   printf("\n Opening parameter file %s \n", parfile);
   
   if ((fptr = fopen(parfile,"r")) == NULL){
       printf("Error! opening file");
   }

   fscanf(fptr,"%d", &maxiter);
   fscanf(fptr,"%d %d", &nx, &ny);
   fscanf(fptr,"%lf %lf %lf", &dx, &dy, &dt);
   fscanf(fptr,"%lf %lf %lf %lf", &w0, &tau0, &K, &D);
   fscanf(fptr,"%lf %lf %lf %lf", &alpha, &gama, &ue, &delta);
   fscanf(fptr,"%lf", &epsi4);
   fscanf(fptr,"%lf %lf", &ja, &theta0);
   fscanf(fptr,"%lf", &amplitude);
   fscanf(fptr,"%lf %lf", &c, &lb);
   fscanf(fptr,"%d %d", &n, &mwalker);
   fscanf(fptr,"%lf", &ri);
   fscanf(fptr,"%lf %lf %lf %lf", &cleqm, &cseqm, &ml, &ms); 
   fscanf(fptr,"%lf %lf", &clinitial, &lambda);

   fclose(fptr); 
  
   return;
}

void initialize(void)
{
	printf("\n Initializing variables");
        
	psi = create2Ddouble(nx, ny);
	psi0 = create2Ddouble(nx, ny);
	gradpsix = create2Ddouble(nx, ny);
	gradpsiy = create2Ddouble(nx, ny);
	u = create2Ddouble(nx, ny);
	u0 = create2Ddouble(nx, ny);
	a = create2Ddouble(nx, ny);
	theta = create2Ddouble(nx, ny);
	w = create2Ddouble(nx, ny);
	tau = create2Ddouble(nx, ny);
        wderivative = create2Ddouble(nx,ny);
	gradw2x = create2Ddouble(nx, ny);
	gradw2y = create2Ddouble(nx, ny);
	s = create2Ddouble(nxcg, nycg);
	m = create2Ddouble(nxcg, nycg);
	h = create2Ddouble(nxcg, nycg);
	h0 = create2Ddouble(nxcg, nycg);
	ucg = create2Ddouble(nxcg, nycg);
        mu = create2Ddouble(nx, ny);
        mu0 = create2Ddouble(nx, ny);
        comp = create2Ddouble(nx, ny);
	
	int i, j;
	double rr;

	for(i=0; i<nx; i++)
	{
		for(j=0; j<ny; j++)
		{
			u[i][j] = 0.0;
                        //mu[i][j] = rt*log(clinitial);
			rr = fabs(sqrt((i-nxhalf)*(i-nxhalf)+(j-nyhalf)*(j-nyhalf)));
			
			if(rr<ri)
                        {
				psi[i][j] = 1.0;
                                mu[i][j] = rt*log(pc*clinitial);
                        }
			else
                        {
				psi[i][j] = 0.0;
                                mu[i][j] = rt*log(clinitial);
                        }
		}
	}
	return;
}

void initializewalker(void)
{
    printf("\n Initializing walker positions");
    int k,l,kk,ll;
    int i, check = 0;
    double x, y;
    
    for(k=nxcghalf;k<nxcg;k++)
    {
        for(l=0;l<nycghalf+1;l++)
        {
            if(s[k][l]==2.0)
            {
                m[k][l] = (double)mwalker; //initialize m
                x = drand48()+(double)k;
                y = drand48()+(double)l;
                
                updatewalker(0,x,y,k,l);
                break;
            }
        }
        if(check == 0)
            break;
    }
    
    int count = 0;
    for(k=k;k<nxcg;k++)
    {
        for(l=l;l<nycghalf+1;l++)
        {
            if(s[k][l]==2.0)
            {
                m[k][l] = (double)mwalker; //initialize m
                for(i=0;i<mwalker;i++)
                {
                    x = drand48()+(double)k;
                    y = drand48()+(double)l;
                    
                    updatewalker(1,x,y,k,l);
                    count+=1;
                }
            }
        }
    }
    printf("\n Initial number of walkers = %d", count);
    return;
}

void boundary(void)
{
	int i, j;	
	printf("\n Applying boundary conditions");
	//no flux boundary conditions
	for(i=0; i<nx; i++)
	{
		psi[i][0] = psi[i][1];
		psi[i][ny-1] = psi[i][ny-2];
		u[i][0] = u[i][1];
		u[i][ny-1] = u[i][ny-2];
                mu[i][0] = mu[i][1];
		mu[i][ny-1] = mu[i][ny-2];
	}
	for(j=0; i<ny; j++)
	{
		u[0][j] = u[1][j];
		u[nx-1][j] = u[nx-2][j];
		psi[0][j] = psi[1][j];
		psi[nx-1][j] = psi[nx-2][j];
                mu[0][j] = mu[1][j];
		mu[nx-1][j] = mu[nx-2][j];
	}
	
	return;
}

void boundarymc(void)
{
    return;
}

void computepsi(void)
{
	printf("\n Computing phase field");
	int i, j;
	double t1, t2, t3, t4, t5;
	double laplacianpsi;
	double tt1, tt2;
	double randnum;
        double mthermal, msolutal;
        
        double igradpsix;
	
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
			gradpsix[i][j] = (psi0[i+1][j]-psi0[i-1][j])*0.5*idx;
			gradpsiy[i][j] = (psi0[i][j+1]-psi0[i][j-1])*0.5*idy; 		
		}
	}

	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
			//anisotropy function
                    igradpsix = 1.0/gradpsix[i][j];
			if(gradpsix[i][j]==0.0)
			{
				if(gradpsiy[i][j] >= 0.0)
					theta[i][j] = PI*0.5;
				else
					theta[i][j] = -PI*0.5;
			}
			else if(gradpsix[i][j]<0.0)
			{
				theta[i][j] = PI+atan(gradpsiy[i][j]*igradpsix);
			}
			else
			{
				if(gradpsiy[i][j] < 0.0)
					theta[i][j] = 2.0*PI+atan(gradpsiy[i][j]*igradpsix);
				else
					theta[i][j] = atan(gradpsiy[i][j]*igradpsix);		
			} 
			
			a[i][j] = 1.0 + epsi4*cos(ja*(theta[i][j]-theta0));
			//tau calculation
			tau[i][j] = tau0*a[i][j];	
			//w calculation
			w[i][j] = w0*a[i][j];
		}
	}	
	
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
			wderivative[i][j] = -w0*epsi4*ja*sin(ja*(theta[i][j]-theta0));
			gradw2x[i][j] = (w[i+1][j]*w[i+1][j]-w[i-1][j]*w[i-1][j])*0.5*idx;
			gradw2y[i][j] = (w[i][j+1]*w[i][j+1]-w[i][j-1]*w[i][j-1])*0.5*idy;
		}
	}
	
	double iPI = 1.0/PI;
	for(i=nxhalf; i<nx-1; i++)
	{
            /*
            if(i==nxhalf)
                printf("\n");
            */
		for(j=1; j<nyhalf+1; j++)
		{
			randnum = drand48()-0.5;
                        
			//psi calculation
			t1 = (w[i][j+1]*wderivative[i][j+1]*gradpsix[i][j+1]-w[i][j-1]*wderivative[i][j-1]*gradpsix[i][j-1])*0.5*idy;
			t2 = -(w[i+1][j]*wderivative[i+1][j]*gradpsiy[i+1][j]-w[i-1][j]*wderivative[i-1][j]*gradpsiy[i-1][j])*0.5*idx;
                        
                        if(itemp==1)
                            mthermal = alpha*atan(gama*(ue-u[i][j]))*iPI;
                        else 
                            mthermal = 0.0;
                        
                        if(icomp==1)
                        {
                            msolutal = (cleqm-cseqm)*(mueqm-mu[i][j])*lambda;
                        }
                        else
                            msolutal = 0.0;
                        /*
                        if(i==nxhalf)
                            printf(" %lf",mueqm-mu[i][j]);
                        */
                        t3 = psi0[i][j]*(1.0-psi0[i][j])*(psi0[i][j]-0.5+mthermal+msolutal+amplitude*randnum); //thermal & solutal field coupling
                        
			laplacianpsi = (2.0*(psi0[i+1][j]+psi0[i-1][j]+psi0[i][j+1]+psi0[i][j-1])+psi0[i+1][j+1]+psi0[i-1][j-1]+psi0[i-1][j+1]+psi0[i+1][j-1]-12.0*psi0[i][j])*idx*idx*0.33333333;
			t4 = w[i][j]*w[i][j]*laplacianpsi;
			t5 = gradw2x[i][j]*gradpsix[i][j]+gradw2y[i][j]*gradpsiy[i][j];
                        
			psi[i][j] = psi0[i][j] + dt*(t1+t2+t4+t5+t3)/tau[i][j];
		}
	}
	return;	
}

void computeu(void)
{
	printf("\n Computing thermal field");
	int i, j;
	double t1, t2;

	for(i=nxhalf; i<nx-1; i++)
	{
            if(i==nxhalf)
                printf("\n");
		for(j=1; j<nyhalf+1; j++)
		{
			//u calculation
			t1 = (2.0*(u0[i+1][j]+u0[i-1][j]+u0[i][j+1]+u0[i][j-1])+u0[i+1][j+1]+u0[i-1][j-1]+u0[i-1][j+1]+u0[i+1][j-1]-12.0*u0[i][j])*idx*idx*0.33333333;
			t2 = K*(psi[i][j]-psi0[i][j]);
			u[i][j] = u0[i][j] + dt*t1 + t2;
                        
                        if(i==nxhalf)
                            printf(" %lf",dt*t1+t2);
		}
	}
	return;	
}

void computecomp(void)
{
    printf("\n Computing solutal field");
    int i,j;
    double lhs;
    double t1, t2;
    double cliquid, csolid, mc;
    
    	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
                    cliquid = exp(mu0[i][j]/rt);
                    csolid = pc*cliquid;
                    comp[i][j] = psi[i][j]*csolid + (1.0-psi[i][j])*cliquid;
                    
                    lhs = csolid*psi[i][j]/rt+cliquid*(1.0-psi[i][j])/rt;
                    
                    
                    mc = psi[i][j]*ms+(1.0-psi[i][j])*ml;
                    
                    
                    t1 = mc*(2.0*(mu0[i+1][j]+mu0[i-1][j]+mu0[i][j+1]+mu0[i][j-1])+mu0[i+1][j+1]+mu0[i-1][j-1]+mu0[i-1][j+1]+mu0[i+1][j-1]-12.0*mu0[i][j])*idx*idx*0.33333333;
                    
                    t2 = (csolid-cliquid)*(psi[i][j]-psi0[i][j]);
                    
                    mu[i][j] = mu0[i][j] + (dt*t1-t2)/lhs;
                    /*
                    cliquid = mu0[i][j];
                    csolid = pc*cliquid;
                    comp[i][j] = psi[i][j]*csolid + (1.0-psi[i][j])*cliquid;
                    
                    lhs = pc*psi[i][j]+(1.0-psi[i][j]);
                    mc = psi[i][j]*ms+(1.0-psi[i][j])*ml;
                    t1 = mc*(2.0*(mu0[i+1][j]+mu0[i-1][j]+mu0[i][j+1]+mu0[i][j-1])+mu0[i+1][j+1]+mu0[i-1][j-1]+mu0[i-1][j+1]+mu0[i+1][j-1]-12.0*mu0[i][j])*idelx*idelx*0.33333333;
                    t2 = (csolid-cliquid)*(psi[i][j]-psi0[i][j]);
                    mu[i][j] = mu0[i][j] + (dt*t1 - t2)/lhs;
                    */
		}
	}
    return;
}

void stepahead(void)
{
	printf("\n Stepping ahead");
        iter = iter+1;
	int i, j;
	for(i=0; i<nx; i++)
	{
		for(j=0; j<ny; j++)
		{
			u0[i][j] = u[i][j];
			psi0[i][j] = psi[i][j];
                        mu0[i][j] = mu[i][j];
		}
	}
	
	for(i=0;i<nxcg;i++)
	{
		for(j=0;j<nycg;j++)
			h0[i][j] = h[i][j];
	}	
	return;
}


void computemc(void)
{
	printf("\n Computing using MC");
	int i, j;
	double t1, t2, t3, t4, t5;
	double laplacianpsi;
	double tt1, tt2;
	double randnum;
        double mthermal, msolutal;
        
        double igradpsix;
        
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
                    if(s[i/n][j/n] == -2.0 || s[i/n][j/n] == -1.0)
                    {
			gradpsix[i][j] = (psi0[i+1][j]-psi0[i-1][j])*0.5*idx;
			gradpsiy[i][j] = (psi0[i][j+1]-psi0[i][j-1])*0.5*idy; 
                    }
		}
	}
        
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
                    if(s[i/n][j/n] == -2.0 || s[i/n][j/n] == -1.0)
                    {
			//anisotropy function
                    igradpsix = 1.0/gradpsix[i][j];
			if(gradpsix[i][j]==0.0)
			{
				if(gradpsiy[i][j] >= 0.0)
					theta[i][j] = PI*0.5;
				else
					theta[i][j] = -PI*0.5;
			}
			else if(gradpsix[i][j]<0.0)
			{
				theta[i][j] = PI+atan(gradpsiy[i][j]*igradpsix);
			}
			else
			{
				if(gradpsiy[i][j] < 0.0)
					theta[i][j] = 2.0*PI+atan(gradpsiy[i][j]*igradpsix);
				else
					theta[i][j] = atan(gradpsiy[i][j]*igradpsix);		
			} 
			
			a[i][j] = 1.0 + epsi4*cos(ja*(theta[i][j]-theta0));
			//tau calculation
			tau[i][j] = tau0*a[i][j];	
			//w calculation
			w[i][j] = w0*a[i][j];
                    }
		}
	}	
	
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
                    if(s[i/n][j/n] == -2.0 || s[i/n][j/n] == -1.0)
                    {
			wderivative[i][j] = -w0*epsi4*ja*sin(ja*(theta[i][j]-theta0));
			gradw2x[i][j] = (w[i+1][j]*w[i+1][j]-w[i-1][j]*w[i-1][j])*0.5*idx;
			gradw2y[i][j] = (w[i][j+1]*w[i][j+1]-w[i][j-1]*w[i][j-1])*0.5*idy;
                    }
		}
	}
	
	double iPI = 1.0/PI;
	for(i=nxhalf; i<nx-1; i++)
	{
		for(j=1; j<nyhalf+1; j++)
		{
                    if(s[i/n][j/n] == -2.0 || s[i/n][j/n] == -1.0)
                    {
			randnum = drand48()-0.5;
			mthermal = alpha*atan(gama*(ue-u0[i][j]))*iPI;
			//psi calculation
			t1 = (w[i][j+1]*wderivative[i][j+1]*gradpsix[i][j+1]-w[i][j-1]*wderivative[i][j-1]*gradpsix[i][j-1])*0.5*idy;
			t2 = -(w[i+1][j]*wderivative[i+1][j]*gradpsiy[i+1][j]-w[i-1][j]*wderivative[i-1][j]*gradpsiy[i-1][j])*0.5*idx;
			t3 = psi0[i][j]*(1.0-psi0[i][j])*(psi0[i][j]-0.5+mthermal+amplitude*randnum);
			//t3 = psi0[i][j]*(1.0-psi0[i][j])*(psi0[i][j]-0.5+m);
			laplacianpsi = (2.0*(psi0[i+1][j]+psi0[i-1][j]+psi0[i][j+1]+psi0[i][j-1])+psi0[i+1][j+1]+psi0[i-1][j-1]+psi0[i-1][j+1]+psi0[i+1][j-1]-12.0*psi0[i][j])*idx*idx*0.33333333;
			t4 = w[i][j]*w[i][j]*laplacianpsi;
			t5 = gradw2x[i][j]*gradpsix[i][j]+gradw2y[i][j]*gradpsiy[i][j];
			psi[i][j] = psi0[i][j] + dt*(t1+t2+t3+t4+t5)/tau[i][j];
			
			//u calculation
			tt1 = (2.0*(u0[i+1][j]+u0[i-1][j]+u0[i][j+1]+u0[i][j-1])+u0[i+1][j+1]+u0[i-1][j-1]+u0[i-1][j+1]+u0[i+1][j-1]-12.0*u0[i][j])*idx*idx*0.33333333;
			tt2 = K*(psi[i][j]-psi0[i][j]);
			u[i][j] = u0[i][j] + dt*tt1 + tt2;
                    }
		}
	}
	
	
	return;	
}


void applysymmetry(void)
{
    int i,j,k,l;
    
    printf("\n Applying symmetry");
    //anticlockwise
    //quadrant 1
    for(i=0;i<nxhalf+1;i++)
    {
        for(j=0;j<nyhalf+1;j++)
        {
            psi[i][j] = psi[nx-1-i][j];
            u[i][j] = u[nx-1-i][j];
            mu[i][j] = mu[nx-1-i][j];
            comp[i][j] = comp[nx-1-i][j];
        }
    }
    
    //quadrant 3
    for(i=nxhalf;i<nx-1;i++)
    {
        for(j=nyhalf;j<ny-1;j++)
        {
            psi[i][j] = psi[nx-1-i][ny-1-j];
            u[i][j] = u[nx-1-i][ny-1-j];
            mu[i][j] = mu[nx-1-i][ny-1-j];
            comp[i][j] = comp[nx-1-i][ny-1-j];
        }
    }
    
    //quadrant 4
    for(i=0;i<nxhalf+1;i++)
    {
        for(j=nyhalf;j<ny-1;j++)
        {
            psi[i][j] = psi[i][ny-1-j];
            u[i][j] = u[i][ny-1-j];
            mu[i][j] = mu[i][ny-1-j];
            comp[i][j] = comp[i][ny-1-j];
        }
    }
    return;
}

void applysymmetrymc(void)
{
    int i,j,k,l;
    
    printf("\n Applying symmetry in MC");
    //anticlockwise
    //quadrant 1
    for(i=0;i<nxhalf+1;i++)
    {
        for(j=0;j<nyhalf+1;j++)
        {
            psi[i][j] = psi[nx-1-i][j];
            u[i][j] = u[nx-1-i][j];
        }
    }
    
    for(k=0;k<nxcghalf+1;k++)
    {
        for(l=0;l<nycghalf+1;l++)
        {
            s[k][l] = s[nxcg-1-k][l];
            h[k][l] = h[nxcg-1-k][l];
            m[k][l] = m[nxcg-1-k][l];
            ucg[k][l] = ucg[nxcg-1-k][l];
        }
    }
    
    //quadrant 3
    for(i=nxhalf;i<nx-1;i++)
    {
        for(j=nyhalf;j<ny-1;j++)
        {
            psi[i][j] = psi[nx-1-i][ny-1-j];
            u[i][j] = u[nx-1-i][ny-1-j];
        }
    }
    
    for(k=nxcghalf;k<nxcg-1;k++)
    {
        for(l=nycghalf;l<nycg-1;l++)
        {
            s[k][l] = s[nxcg-1-k][nycg-1-l];
            h[k][l] = h[nxcg-1-k][nycg-1-l];
            m[k][l] = m[nxcg-1-k][nycg-1-l];
            ucg[k][l] = ucg[nxcg-1-k][nycg-1-l];
        }
    }
    
    //quadrant 4
    for(i=0;i<nxhalf+1;i++)
    {
        for(j=nyhalf;j<ny-1;j++)
        {
            psi[i][j] = psi[i][ny-1-j];
            u[i][j] = u[i][ny-1-j];
        }
    }
    
    for(k=0;k<nxcghalf+1;k++)
    {
        for(l=nycghalf;l<nycg-1;l++)
        {
            s[k][l] = s[k][nycg-1-l];
            h[k][l] = h[k][nycg-1-l];
            m[k][l] = m[k][nycg-1-l];
            ucg[k][l] = ucg[k][nycg-1-l];
        }
    }
    return;
}

void updatecg(void)
{
    printf("\n Updating coarse grid");
    
	int i,j,k,l,kk,ll;
	int check=0;
        double d;

        //random initialization of liquid cells
	for(k=0;k<nxcg;k++)
	{
		for(l=0;l<nycg;l++)
			s[k][l] = 2.0; //random initialization
	}
        
	//solid cell
	for(k=0;k<nxcg;k++)
	{
		for(l=0;l<nycg;l++)
		{
			for(i=n*k;i<n*(k+1);i++)
			{
				for(j=n*l;j<n*(l+1);j++)
				{
					if(psi[i][j]>0.5)
                                        {
						check = 1;
                                        }
				}
			}
			if(check==1)
                        {
                            s[k][l] = -2.0;
                        }
			
			check = 0;
		}
	}
	
	//buffer cell
        for(k=0;k<nxcg;k++)
        {
            for(l=0;l<nycg;l++)
            {
                if(s[k][l] == -2.0) //solid
                {
                    for(kk=0;kk<nxcg;kk++)
                    {
                        for(ll=0;ll<nycg;ll++)
                        {
                            if(kk!=k && ll!=l)
                            {
                                d = fabs(sqrt((k-kk)*(k-kk)+(l-ll)*(l-ll)));
                                if(d<lb && s[kk][ll]!=-2.0)
                                    s[kk][ll] = -1.0; 
                            }
                        }
                    }
                }
            }
        }
        
	//conversion cell
	for(k=1;k<nxcg-1;k++)
	{
            if(s[k][0]==2.0 || s[k][ny-1]==2.0 || s[0][k]==2.0 || s[nxcg-1][k]==0.0)
            {
		if(s[k-1][0]==-1.0||s[k+1][0]==-1.0||s[k-1][1]==-1.0||s[k+1][1]==-1.0||s[k][1]==-1.0)
			s[k][0] = 0.0;
		else if(s[k-1][nycg-1]==-1.0||s[k+1][nycg-1]==-1.0||s[k-1][nycg-2]==-1.0||s[k+1][nycg-2]==-1.0||s[k][nycg-2]==-1.0 && s[k][nycg-1]==2.0)
			s[k][nycg-1] = 0.0;
		else if(s[0][k-1]==-1.0||s[0][k+1]==-1.0||s[1][k-1]==-1.0||s[1][k+1]==-1.0||s[1][k]==-1.0 && s[0][k]==2.0)
			s[0][k] = 0.0;
		else if(s[nxcg-1][k-1]==-1.0||s[nxcg-1][k+1]==-1.0||s[nxcg-2][k-1]==-1.0||s[nxcg-2][k+1]==-1.0||s[nxcg-2][k]==-1.0&& s[nxcg-1][k]==2.0)
			s[nxcg-1][k] = 0.0;
            }
            
            for(l=1;l<nycg-1;l++)
            {
                if(s[k][l] == 2.0)
                {
                    if(s[k-1][l]==-1.0 || s[k+1][l] ==-1.0 || s[k][l-1] ==-1.0 || s[k][l+1] ==-1.0 || s[k+1][l+1] ==-1.0 || s[k+1][l-1] ==-1.0 || s[k-1][l+1] ==-1.0 || s[k-1][l-1] ==-1.0 && s[k][l]==2.0)
			s[k][l] = 0.0;
                }
            }
	}
	
	if(s[0][0]==2.0 || s[0][ny-1]==2.0 || s[nxcg-1][0]==2.0 || s[nxcg-1][nycg-1]==2.0)
        {
            if(s[0][1]==-1.0||s[1][0]==-1.0||s[1][1]==-1.0 && s[0][0]==2.0)
		s[0][0] = 0.0;
            if(s[0][nycg-2]==-1.0||s[1][nycg-1]==-1.0||s[1][nycg-2]==-1.0 && s[0][nycg-1]==2.0)
		s[0][nycg-1] = 0.0;
            if(s[nxcg-2][0]==-1.0||s[nxcg-1][1]==-1.0||s[nxcg-2][1]==-1.0 && s[nxcg-1][0]==2.0)
		s[nxcg-1][0] = 0.0;
            if(s[nxcg-1][nycg-2]==-1.0||s[nxcg-2][nycg-1]==-1.0||s[nxcg-2][nycg-2]==-1.0 && s[nxcg-1][nycg-1]==2.0)
		s[nxcg-1][nycg-1] = 0.0;
        }
	
	return;
}

void updatecgu(void)
{
    printf("\n Updating temperature field on coarse grid");
	int k,l;
	int i,j;
	double imwalker = 1.0/(double)mwalker;
        
	for(k=nxcghalf;k<nxcg;k++)
	{
		for(l=0;l<nycghalf+1;l++)
		{
			if(s[k][l]==2.0)
				ucg[k][l] = -delta*(1.0-m[k][l]*imwalker);

			else if(s[k][l]==0.0)
                        {
                            for(i=n*k;i<n*(k+1);i++)
                            {
                                for(j=n*l;j<n*(l+1);j++)	
					u[i][j] = ucg[k][l];
                            }
                        }
                        else if(s[k][l]==-1.0)
                        {
                            ucg[k][l] = delta*(m[k][l]-mwalker+h[k][l]/hc)*imwalker;
                            for(i=n*k;i<n*(k+1);i++)
                            {
                                for(j=n*l;j<n*(l+1);j++)	
					u[i][j] = ucg[k][l];
                            }
                        }
		}
	}
	return;
}

void updatecgh(void)
{
    printf("\n Updating enthalpy field on coarse grid");
	int k,l;
	int i,j;
	double deltau;
	
	for(k=nxcghalf;k<nxcg;k++)
	{
		for(l=1;l<nycghalf+1;l++)
		{
			if(s[k][l]==0.0)
			{
				if(s[k-1][l]==-1.0)
				{
					i = (k-1)*n+1;
					for(j=(l-1)*n+1;j<l*n;j++)
						deltau = u[i-1][j]-u[i][j];
				}
				else if(s[k+1][l]==-1.0)
				{
					i = k*n-1;
					for(j=(l-1)*n+1;j<l*n;j++)
						deltau = u[i+1][j]-u[i][j];
				}
				else if(s[k][l-1]==-1.0)
				{
					j = (l-1)*n+1;
					for(i=(k-1)*n+1;i<k*n;i++)
						deltau = u[i][j-1]-u[i][j];	
				}
				else if(s[k][l+1]==-1.0)
				{
					j = l*n-1;
					for(i=(k-1)*n+1;i<k*n;i++)
						deltau = u[i][j+1]-u[i][j];
				}
				
				h[k][l] = h0[k][l] + dt*(deltau)*idx*idx;
				if(h[k][l]>hc)
                                {
                                    h[k][l] -= hc; 
                                    updatewalker(2,0.0,0.0,k,l); //create walker with random position
                                }
                                else if(h[k][l]<-hc)
                                {
                                    h[k][l] += hc;
                                    updatewalker(-1,0.0,0.0,k,l); //delete walker
                                }
			}
			else
                            h[k][l] = 0.0;
		}
		
	}
	
	return;
}

void updatewalker(int choice, double x, double y, int k, int l)
{
    int i,j,h;
    int kk,ll;
    double length, dcb, dcbmin = 1000.0;
    int checkcb = 0;
    
    struct list *tempwalker, *newwalker;
    
    if(choice == 2)
    {
        //add walkers to the existing list
        printf("\n Adding walkers to the list");
        m[k][l] += 1.0;
        
        newwalker = (struct list *) malloc( sizeof(struct list) );

        if(backbone[iter]->next!=NULL)
        {
            bufferwalker = backbone[iter]->next;
            backbone[iter]->next = newwalker;
            newwalker->next = bufferwalker;
        }
        else
        {
            backbone[iter]->next = newwalker;
            newwalker->next = NULL;
        }
        
        newwalker->x = drand48()+(double)k;
        newwalker->y = drand48()+(double)l;
        
        for(kk=nxcghalf;kk<nxcg;kk++)
        {
            for(ll=0;ll<nycghalf+1;ll++)
            {
                if(s[kk][ll]==0.0)
                {
                    dcb = fabs(sqrt((k-kk)*(k-kk)+(l-ll)*(l-ll)));
                    if(dcb<dcbmin && dcb>0.0)
                        dcbmin = dcb;
                }
            }
        }
        
        length = c*dcbmin;
        newwalker->nt = floor(length*length*0.5*idt*iD);
        
        if(newwalker->nt==0)
            newwalker->nt = 1;
                    
        for(h=1;h<maxiter/newwalker->nt;h++)
        {
            if(backbone[h*newwalker->nt]->next == NULL)
            {
                tempwalker = (struct list *) malloc( sizeof(struct list) );
                tempwalker->x = newwalker->x;
                tempwalker->y = newwalker->y;
                tempwalker->nt = newwalker->nt;
                tempwalker->next = NULL;
                backbone[h*newwalker->nt]->next = tempwalker;
            }
            else
            {
                
                tempwalker = (struct list *) malloc( sizeof(struct list) );
                tempwalker->x = newwalker->x;
                tempwalker->y = newwalker->y;
                tempwalker->nt = newwalker->nt;
                
                bufferwalker = backbone[h*newwalker->nt]->next;
                backbone[h*newwalker->nt]->next = tempwalker;
                backbone[h*newwalker->nt]->next->next = bufferwalker;
                //tempwalker->next = bufferwalker;
            }
        }
        newwalker = NULL;
    }
    else if(choice == -1)
    {
        //remove a walker
        printf("\n Removing walkers from list");
        countwalkerupdates+=1;
        m[k][l] -= 1.0;
        /*
        tempwalker = (struct list *) malloc( sizeof(struct list) );
        tempwalker = backbone[iter]->next;
        while(tempwalker!=NULL && tempwalker->next != walker)
        {
            tempwalker = tempwalker->next;
        }
        
        bufferwalker = walker;
        tempwalker->next = walker->next;
        free(bufferwalker);
        free(tempwalker);
        */
    }
    else if(choice == 1 || choice == 0)
    {
        printf("\n Creating walker list");
        newwalker = (struct list *) malloc( sizeof(struct list) );
        
        if(choice == 1)
        {
            bufferwalker = backbone[0]->next;
            backbone[0]->next = newwalker;
            newwalker->next = bufferwalker;
        }
        else 
        {
            backbone = createpointerarray();
            backbone[0]->next = newwalker;
            newwalker->next = NULL;
        }
        newwalker->x = x;
        newwalker->y = y;
        
        for(kk=nxcghalf;kk<nxcg;kk++)
        {
            for(ll=0;ll<nycghalf+1;ll++)
            {
                if(s[kk][ll]==0.0)
                {
                    dcb = fabs(sqrt((k-kk)*(k-kk)+(l-ll)*(l-ll)));
                    if(dcb<dcbmin && dcb>0.0)
                        dcbmin = dcb;
                }
            }
        }

        length = c*dcbmin;
        newwalker->nt = floor(length*length*0.5*idt*iD);

        if(newwalker->nt==0)
            newwalker->nt = 1;
    }
    
    return;
}


void updatewalkerlist(void)
{    
    printf("\n Updating walker list based on position of coarse grid");

    int k,l;
    int count=0;
    walker = backbone[iter]->next;

    if(walker!=NULL)
    {
        do
        {
            k = (int)floor(walker->x);
            l = (int)floor(walker->y);

            if(s[k][l] == -2.0 || s[k][l] == -1.0)
            {
                m[k][l] -= 1.0;
                printf("\n Removing walkers from fine grid region");
                count+=1;
                
                if(walker != backbone[iter]->next)
                {
                    //node to be removed is not a part of the backbone
                    poswalker = backbone[iter]->next;
                    while(poswalker!=NULL && poswalker->next != walker)
                    {
                        poswalker = poswalker->next;
                    }
                    
                    bufferwalker = walker;
                    poswalker->next = walker->next;
                    free(bufferwalker);
                    walker = poswalker;
                }
                else
                {
                    //node to be removed is in the backbone
                    bufferwalker = walker;
                    backbone[iter]->next = walker->next;
                    walker = backbone[iter]->next;
                }
            }
            
            walker = walker->next;
        }while(walker!=NULL);
    }
    
    printf("\n Number of times walkers were removed from fine grid region = %d", count);
    return;
}

void createbackbone(void)
{
    printf("\n Creating backbone list");
    int k,l;
    int i,j,h;
    int kk,ll;
    double length, dcb, dcbmin = 1000.0;
    struct list *tempwalker;
    walker = backbone[0]->next;
    
    do{
        
        for(h=1;h<maxiter/walker->nt;h++)
        {
            if(backbone[h*walker->nt]->next == NULL)
            {
                tempwalker = (struct list *) malloc( sizeof(struct list) );
                tempwalker->x = walker->x;
                tempwalker->y = walker->y;
                tempwalker->nt = walker->nt;
                tempwalker->next = NULL;
                backbone[h*walker->nt]->next = tempwalker;
            }
            else
            {
                
                tempwalker = (struct list *) malloc( sizeof(struct list) );
                tempwalker->x = walker->x;
                tempwalker->y = walker->y;
                tempwalker->nt = walker->nt;
                
                bufferwalker = backbone[h*walker->nt]->next;
                backbone[h*walker->nt]->next = tempwalker;
                backbone[h*walker->nt]->next->next = bufferwalker;
            }
        }
        walker = walker->next;
    }while(walker != NULL);
    
    return;
}

void updatewalkerpos(void)
{
    printf("\n Updating walker position");
    
    double randnum, length;
    int k_old, l_old;

    walker = backbone[iter]->next;

    if(walker!=NULL)
    {
        do
        {
            randnum = drand48();
            k_old = floor(walker->x);
            l_old = floor(walker->y);
            length = fabs(sqrt(2.0*D*dt*(double)(walker->nt)));

            walker->x += length*randnum;
            randnum = drand48();
            walker->y += length*randnum;

            if((int)floor(walker->x)!=k_old || (int)floor(walker->y)!=l_old)
            {
                //update m
                m[k_old][l_old] -= 1.0;
                m[(int)floor(walker->x)][(int)floor(walker->y)] += 1.0;
            }
            
            countwalkerupdates+=1;
            walker = walker->next;
        }while(walker!=NULL);
    }
    
    return;
}

void printoutput(int r)
{
	int i, j;
  	FILE *fp;
  	char filename[100];
  
  	sprintf (filename, "output%d.dat", r);
  	printf("\n Starting output %s \n",filename);

  	fp=fopen(filename,"w");
  	fprintf(fp, "psi=[");
  
 	for ( i = 0; i<nx; i++)
  	{
    	for(j=0; j<ny; j++)
    		fprintf(fp, "%lf  ", psi[i][j]);
    
    	fprintf(fp,"\n");
  	}
  	fprintf(fp, "];\n");
  	
  	fprintf(fp, "u=[");
  
 	for ( i = 0; i<nx; i++)
  	{
    	for(j=0; j<ny; j++)
    		fprintf(fp, "%lf  ", u[i][j]);
    
    	fprintf(fp,"\n");
  	}
  	fprintf(fp, "];\n");
        
        fprintf(fp, "comp=[");
  
 	for ( i = 0; i<nx; i++)
  	{
    	for(j=0; j<ny; j++)
    		fprintf(fp, "%lf  ", comp[i][j]);
    
    	fprintf(fp,"\n");
  	}
  	fprintf(fp, "];\n");
        
        fprintf(fp, "mu=[");
  
 	for ( i = 0; i<nx; i++)
  	{
    	for(j=0; j<ny; j++)
    		fprintf(fp, "%lf  ", mu[i][j]);
    
    	fprintf(fp,"\n");
  	}
  	fprintf(fp, "];\n");
        
  	fprintf(fp, "s=[");
  
 	for ( i = 0; i<nxcg; i++)
  	{
    	for(j=0; j<nycg; j++)
    		fprintf(fp, "%lf  ", s[i][j]);
    
    	fprintf(fp,"\n");
  	}
  	fprintf(fp, "];\n");
        
  	fprintf(fp,"\nxset(\"colormap\",jetcolormap(64));");
	fprintf(fp,"\ncolorbar(0,1);");
	fprintf(fp,"\nSgrayplot(1:%d,1:%d,psi);",nx,ny);
        //fprintf(fp,"\ncontour2d(1:%d,1:%d,psi,3);",nx,ny);
  	fclose(fp);

/*  	
  	printf("\n Psi = \n");
  	for(i=0; i<nx; i++)
  		printf(" %lf", psi[i][ny/2]);
  	printf("\n u = \n");
  	for(i=0; i<nx; i++)
  		printf(" %lf", u[i][ny/2]);
  	printf("\n");
*/
	return;
}

void rawoutput(int r)
{
    
    return;
}

void freewalkermemory(void)
{
    printf("\n Freeing used walker memory");
    walker = backbone[iter]->next;
    
    if(walker!=NULL)
    {
        do
        {
            bufferwalker = walker;
            walker = walker->next;
            free(bufferwalker);
        }while(walker!=NULL);
    }
    return;
}

void freememory(void)
{
	printf("\n Freeing memory");
	free(psi);
	free(psi0);
	free(gradpsix);
	free(gradpsiy);
	free(u);
	free(u0);
	free(a);
	free(theta);
	free(w);
	free(tau);
	free(wderivative);
	free(gradw2x);
        free(gradw2y);
        free(walker);
        free(bufferwalker);
        free(poswalker);
        free(mu);
        free(mu0);
        free(comp);
	
	return;
}
