/* INTEGRATOR RK4

 The code compute Positions and Velocities of Sun and Earth starting from Initial condintions and using the RK4 method.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define me 5.9736E24
#define ms 1.9891E30
#define a 149598261000.0
#define e 0.0167112303531389
#define G 6.67428E-11
#define I1 8.008E37
#define I3 8.034E37
#define Tday0 86164.1
#define theta0 23.45
#define PI 3.1415926535897932384626433832795028841971693

double Tday = Tday0*1;


void AngularMomentum(double *y, double *Ls, double *Le)
{
  //Angular momentum SUN
  Ls[1]=ms*(y[4]*y[11]-1.0*y[5]*y[10]);
  Ls[2]=ms*(y[5]*y[9]-1.0*y[3]*y[11]);
  Ls[3]=ms*(y[3]*y[10]-1.0*y[4]*y[9]);
  
 // ngular Momentum Earth
  Le[1]=me*(y[1]*y[8]-1.0*y[2]*y[7]);
  Le[2]=me*(y[2]*y[6]-1.0*y[0]*y[8]);
  Le[3]=me*(y[0]*y[7]-1.0*y[1]*y[6]);
}


void SeW(double *y, double *Spine)
{
  double Stemp[3]= {Spine[0],Spine[1],Spine[2]}; // to save the values in the body coordinate system
  
  Spine[0] = Stemp[0]*(cos(y[14])*cos(y[12])-1.0*cos(y[13])*sin(y[12])*sin(y[14]))+Stemp[1]*(-1.0*sin(y[14])*cos(y[12])-1.0*cos(y[13])*sin(y[12])*cos(y[14]))+Stemp[2]*(sin(y[13])*sin(y[12])); //world
  
  Spine[1] = Stemp[0]*(cos(y[14])*sin(y[12])+cos(y[13])*cos(y[12])*sin(y[14]))+Stemp[1]*(-1.0*sin(y[14])*sin(y[12])+cos(y[12])*cos(y[13])*cos(y[14]))+Stemp[2]*(-1.0*sin(y[13])*cos(y[12])) ; // world
  
  Spine[2] = Stemp[0]*(sin(y[13])*sin(y[14]))+Stemp[1]*(sin(y[13])*cos(y[14]))+Stemp[2]*(cos(y[13])) ; // world
 

}

void OmegaFunction(double *y,double *omeg, double *Spine)
{
  // 16 theta dot .. 13 theta
  omeg[1]=y[16]*cos(y[14])+y[15]*sin(y[13])*sin(y[14]);
  omeg[2]=-y[16]*sin(y[14])+1.0*y[15]*sin(y[13])*cos(y[14]);
  omeg[3]=y[15]*cos(y[13])+y[17];
 
  Spine[1]=I1*omeg[1]; //in the body system
  Spine[2]=I1*omeg[2]; //body
  Spine[3]=I3*omeg[3]; //body
  
}

void gravforce(double *y, int n, double t, double *der)
{
  for(int i=0; i<6; i++) der[i]=y[6+i];//for position of sun and earth, the derivative is the velocity
  
  double dv[3]={y[3]-y[0], y[4]-y[1], y[5]-y[2]};//dv vector which is R_s - R_e
  double dis=pow(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2],0.5);
  double f[3];
  double d3=((dv[0]*sin(y[12])-dv[1]*cos(y[12]))*sin(y[13])+dv[2]*cos(y[13]));//used in the following
  f[0]=-(G*ms*me/(dis*dis*dis))*dv[0]+((3.0*G*ms)/(pow(dis,5.0)))*(I1-I3)*((dv[0]/2.0)-((2.5*d3*d3*dv[0])/(dis*dis))+d3*sin(y[12])*sin(y[13]));
  f[1]=-(G*ms*me/(dis*dis*dis))*dv[1]+((3.0*G*ms)/(pow(dis,5.0)))*(I1-I3)*((dv[1]/2.0)-((2.5*d3*d3*dv[1])/(dis*dis))-d3*cos(y[12])*sin(y[13]));
  f[2]=-(G*ms*me/(dis*dis*dis))*dv[2]+((3.0*G*ms)/(pow(dis,5.0)))*(I1-I3)*((dv[2]/2.0)-((2.5*d3*d3*dv[2])/(dis*dis))+d3*cos(y[13]));
  
  for(int i = 0; i<3; ++i) der[6+i]=-f[i]/me;//for velocity of sun and earth, the derivative is the acceleration
  for(int i = 0; i<3; ++i) der[9+i]=f[i]/ms;
  
  for(int i=0; i<3; i++) der[12+i]=y[15+i];//the derivative for the angle is the angular velocity
  double c0=(2.0*PI)/(Tday);//used in the following
  der[15]=-2.0*y[15]*y[16]*cos(y[13]);
  der[15]=der[15]+(I3/I1)*c0*y[16];
  der[15]=der[15]+((3.0*G*ms*(I1-I3)*d3)/(pow(dis,5.0)*I1))*(dv[0]*cos(y[12])+dv[1]*sin(y[12]));
  der[15]=der[15]/(sin(y[13]));
  
  der[16]=y[15]*y[15]*sin(y[13])*cos(y[13]);
  der[16]=der[16]-(I3/I1)*c0*y[15]*sin(y[13]);
  der[16]=der[16]+((3.0*G*ms*(I1-I3)*d3)/(pow(dis,5.0)*I1))*((dv[0]*sin(y[12])-dv[1]*cos(y[12]))*cos(y[13])-dv[2]*sin(y[13]));
  
  der[17]=0;//the derivative for the angular velocity in psi is not need, since dot{psi} can be computed by c0-dot{phi}*cos(theta)
}


void integrate(double* yo, int n, double* yn, void f(double*, int, double, double*), double t, double dt)//integrate by RK4
{
  double *k1,*k2,*k3,*k4,*yt;
  k1=(double*)malloc(n*sizeof(double));
  k2=(double*)malloc(n*sizeof(double));
  k3=(double*)malloc(n*sizeof(double));
  k4=(double*)malloc(n*sizeof(double));
  yt=(double*)malloc(n*sizeof(double));
  
  f(yo, n, t, k1);
  for(int i=0; i<n-1; i++) yt[i]=yo[i]+dt*k1[i]/2.0;
  f(yt, n, t+dt/2.0, k2);
  for(int i=0; i<n-1; i++) yt[i]=yo[i]+dt*k2[i]/2.0;
  f(yt, n, t+dt/2.0, k3);
  for(int i=0; i<n-1; i++) yt[i]=yo[i]+dt*k3[i];
  f(yt, n, t+dt, k4);
  
  for(int i=0; i<n-1; ++i) yn[i]=yo[i]+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
  yn[17]=(2*PI/Tday)-yn[15]*cos(yn[13]);
  
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(yt);
}

int main(int argc, char** argv)
{
  clock_t start, finish;
  double  duration;
  start = clock();
  
  int n = 18;
  double *yold, *ynew, *Omega, *SpinEARTH, *LEarth, *LSun;
  yold=(double*)malloc(n*sizeof(double));
  ynew=(double*)malloc(n*sizeof(double));
  Omega=(double*)malloc(3*sizeof(double));
  SpinEARTH=(double*)malloc(3*sizeof(double));
  LEarth=(double*)malloc(3*sizeof(double));
  LSun=(double*)malloc(3*sizeof(double));
  
  FILE *result=NULL;
  result=fopen("result.txt","w"); //result file
  
  int nstep = atof(argv[1]);//the first input is number of step
  int numberyear = atof(argv[2]);//the second is number of year
  
  //double Tyear = 31552800.0;
  double Tyear = 2*PI*pow(a,1.5)/(pow(G*(me+ms),0.5));
  printf("time year: %.20lf \n", Tyear);
  double time = numberyear * Tyear;
  double dt = time/nstep;
  
  double r0=(a*(1+e)*ms)/(ms+me);//initial position for earth
  double KE=(G*ms*me/(a))*(1/(1+e)-0.5);
  double v0=pow((2*KE)/(me*(1+me/(ms))),0.5);//initial velocity for earth
  
  double re0[3]={r0,0,0};
  double ve0[3]={0,v0,0};
  double rs0[3]={(-me/(ms))*r0,0,0};//initial position of sun
  double vs0[3]={0,(-me/(ms))*v0,0};//initial velocity of sun
  double angle0[3]={0,(theta0*PI)/180.0,0};//euler angle (phi,theta,psi)
  double angle_v0[3]={0,0,(2.0*PI)/(Tday)};//angular velociy in euler angel
  
  for(int i=0; i<3; i++) yold[0+i]=re0[i];//position x y z of earth
  for(int i=0; i<3; i++) yold[3+i]=rs0[i];//position x y z of sun
  for(int i=0; i<3; i++) yold[6+i]=ve0[i];//velocity of earth
  for(int i=0; i<3; i++) yold[9+i]=vs0[i];//velocity of sun
  for(int i=0; i<3; i++) yold[12+i]=angle0[i];//euler angle
  for(int i=0; i<3; i++) yold[15+i]=angle_v0[i];//Euler angle derivative
  
  double ttt = 0;///ttt will be the calculated period of precession
  double d_phi=10;
  for(int step = 0; step < nstep; ++step)
  {
    double t=step*dt;
    integrate(yold, n, ynew, gravforce, t, dt);
    
    // Printing data on a File:
   double diff[3] ;
   double differenza;
    
    
    if(step%300000==0 && n!=0)
    {
                                                                                                                                                        
      OmegaFunction(ynew,Omega,SpinEARTH);
      SeW(ynew,SpinEARTH);
      AngularMomentum(ynew,LSun,LEarth);
    
      
      
      for(int i=0; i<3; ++i)
      {
        diff[i]=(SpinEARTH[1]+LSun[i]+LEarth[i])*10E-20 ;
              }
  
      differenza=sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
    //printf("| %.20e  ", differenza);
    //4printf("\n");
      //fprintf(result, "%.20f ", t+dt);
      //for(int i=0; i<n; ++i) fprintf(result, "%.20e ", ynew[i]); //Printing on file the result stored in ynew
      //fprintf(result," \n");
      //for(int i=0; i<n; ++i) printf( "%.20e  ", ynew[i]); //Printing on file the result stored in ynew
     // printf(" \n");
      
    }

    
    
    if(fabs(fabs(ynew[12])-2*PI) < d_phi)
    {
      d_phi=fabs(fabs(ynew[12])-2*PI);

      ttt=t;
    }
    
    for(int i=0; i<n; ++i) yold[i] = ynew[i];
  }
  printf("period of precession: ");
  printf("%.2f",ttt/Tyear);
  printf("\n");
  
  double ReTyear = pow(ynew[0]*ynew[0]+ynew[1]*ynew[1]+ynew[2]+ynew[2],0.5);
  double error = (ReTyear - r0)/(r0);
  
  printf("the original position of earth: ");
  for(int i=0; i<3; ++i) printf("%e \t", re0[i]);
  printf("\n");
  
  printf("the final position of earth: ");
  for(int i=0; i<3; ++i) printf("%e \t", ynew[i]);
  printf("\n");
  
  printf("the time step dt: %lf \n", dt);
  
  printf("The error for the position after computation: ");
  printf("%.35f \n", error);
  
  free(yold);
  free(ynew);
  fclose(result);
  
  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;//comput program running time
  printf( "Time to do integration is: ");
  printf( "%18.10f seconds\n", duration );
  
  return 0;
}
