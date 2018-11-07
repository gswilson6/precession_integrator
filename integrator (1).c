/* Wilson Stewart
   November 11, 2015

   This program numerically integrates
   the force on the Earth and the force
   on the Sun equations to determine
   the path taken in the Earth-Sun
   system by the Earth and Sun.

*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define a  149598261000.0
#define e  0.0167112303531389
#define G  6.67428*pow(10,-11)
#define ms 1.9891*pow(10,30)
#define me 5.9736*pow(10,24)

void gravforce(double *y, int n, double t, double *der)
{
  /*This section of code defines the gravitational
    force acting in the Earth-Sun system.
  */

  for(int i=0; i<6; i++) der[i]=y[6+i];

  double dv[3]={y[3]-y[0], y[4]-y[1], y[5]-y[2]};//distance vector between the earth and sun
  double d3=pow(dv[0]*dv[0]+dv[1]*dv[1]+dv[2]*dv[2],1.5);//magnitude of the distance vector dv
  double f[3]; for(int i = 0; i<3; ++i) f[i]=G*ms*me*dv[i]/d3;//force vector

  for(int i = 0; i<3; ++i) der[6+i]=f[i]/me;//defining the acceleration of the earth system
  for(int i = 0; i<3; ++i) der[9+i]=-f[i]/ms;//defining the acceleration of the sun system
}

void integrate_eu(double* yo, int n, double* yn, void f(double*, int, double, double*), double t, double dt)
{
  /*This section of code runs the Euler integration
    method to define the changing positions of the Earth
    and Sun.
  */

  double *der = (double*)malloc(n*sizeof(double));//allocating memory for the *der array
  f(yo, n, t, der);//defining the force between 2 objects
  for(int i = 0; i<n; ++i) yn[i] = yo[i]+ dt*der[i];//gives us the new y after integration using Euler method
}

void integrate(double* yo, int n, double* yn, void f(double*, int, double, double*), double t, double dt)
{

  /*this section of code runs the Runge-Kutte 4 integrator
    to define the changing positions of the Earth and Sun.
    This process is more precise than the Euler integration
    method
  */


  double *k1,*k2,*k3,*k4,*yt;
  k1=(double*)malloc(n*sizeof(double));
  k2=(double*)malloc(n*sizeof(double));
  k3=(double*)malloc(n*sizeof(double));
  k4=(double*)malloc(n*sizeof(double));
  yt=(double*)malloc(n*sizeof(double));

  f(yo, n, t, k1);
  for(int i=0; i<n; i++) yt[i]=yo[i]+dt*k1[i]/2;
  f(yt, n, t+dt/2, k2);
  for(int i=0; i<n; i++) yt[i]=yo[i]+dt*k2[i]/2;
  f(yt, n, t+dt/2, k3);
  for(int i=0; i<n; i++) yt[i]=yo[i]+dt*k3[i];
  f(yt, n, t+dt, k4);

  for(int i=0; i<n; ++i) yn[i]=yo[i]+dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;

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

  int n = 12;
  double *yold, *ynew;
  yold=(double*)malloc(n*sizeof(double));
  ynew=(double*)malloc(n*sizeof(double));
  double dt = 3600;//one step is equal to one hour in time on Earth
  int nstep = 24*365;//sets our time steps so we are able to see one whole year (1 orbit)

  double initEarth_r = a*(1+e)*(ms/(ms+me));//defines the initial position of the Earth
  double KE_Earth = (G*ms*me/(a))*(1/(1+e)-0.5);//defines the kinetic energy of the Earth at an initial time
  double initEarth_v = sqrt(2*KE_Earth/(me*(1+me/(ms))));//defines the initial condition for the velocity of the Earth

  double initSun_r = -(me/(ms))*initEarth_r;//defines initial condition for the position of the Sun
  double initSun_v = -(me/(ms))*initEarth_v;//defines initial condition for the velocity of the Sun

  double re0[3]={initEarth_r,0,0};//sets initial condition for position of Earth
  double rs0[3]={initSun_r,0,0};//sets initial condition for position of Sun
  double ve0[3]={0,initEarth_v,0};//sets initial condition for velocity of Earth
  double vs0[3]={initSun_v,0,0};//sets initial condition for velocity of Sun


  //Storing our initial conditions
  for(int i=0; i<3; i++) yold[0+i]=re0[i];
  for(int i=0; i<3; i++) yold[3+i]=rs0[i];
  for(int i=0; i<3; i++) yold[6+i]=ve0[i];
  for(int i=0; i<3; i++) yold[9+i]=vs0[i];


  for(int step = 0; step < nstep; ++step)
  {
    double t=step*dt;//defines our time used in the integrate function
    integrate(yold, n, ynew, gravforce, t, dt);//calls the integration function
    printf("%.2f ", t+dt);
    for(int i=0; i<n; ++i) printf("%e ", ynew[i]);//prints the ynew in an array for mathematica to read
    printf("\n");
    for(int i=0; i<n; ++i) yold[i] = ynew[i];//defines the new y
  }

  finish = clock();
  duration = (double)(finish - start) / CLOCKS_PER_SEC;//comput program running time
  printf( "Time to do integration is: ");
  printf( "%18.10f seconds\n", duration );

  free(yold);
  free(ynew);
  return 0;
}
