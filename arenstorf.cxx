#include <cmath>
#include <fstream>
#include <cstdlib>
#include <iostream>

using namespace std;

void f(const double* const y, double* const k, const double& nu);
double maximum(double const a, double const b);

int main(){
  const int N=4;
  double* swap;
  

  //Mass ratio
  const double nu=0.012277471;
  //Tol for stepsize control
  double tola=0.0001;
  double max;
  double* y = new double[N];
  double k1[N]; 
  double k2[N];  
  double k3[N];
  double k4[N];
  double k5[N];
  double k6[N];
  double k7[N];
  double* temp = new double[N];
  double chi[N];
  //Time for one orbit 
  double To =20;
  //stating stepsize
  double dt = 1e-5;
  //time 
  double t = 0;
  
  ofstream out("data.txt");
  //initial conditions 
  y[0]=0.994;
  y[1]=0;
  y[2]=0;
  y[3]=-2.00158510637908;
       
  //calc k_i, i=1,...,7
 while(t<To)
 {
    out<<t<<"\t"<<dt<<"\t"<<y[0]<<"\t"<<y[1]<<"\t"<<y[2]<<"\t"<<y[3]<<endl;
    //calc k1
    f(y,k1,nu);
    
    
    //calc k2
    for(int i=0;i<N;i++)
      temp[i]=y[i]+dt*(1./5)*k1[i];
    
    
    f(temp,k2,nu);
    

    //calc k3
    for(int j=0;j<N;j++)
      temp[j]=y[j]+dt*(3./40*k1[j]+9./40*k2[j]);
  
    f(temp,k3,nu);

    
    //calc k4
    for(int k=0;k<N;k++)
      temp[k]=y[k]+dt*(44./45*k1[k]-56./15*k2[k]+32./9*k3[k]);
  
    f(temp,k4,nu);

    //calc k5
    for(int l=0;l<N;l++)
      temp[l]=y[l]+dt*(19372./6561*k1[l]-25360./2187*k2[l]+64448./6561*k3[l]-212./729*k4[l]);
    
    f(temp,k5,nu);
    
    //calc k6
    for(int m=0;m<N;m++)
      temp[m]=y[m]+dt*(9017./3168*k1[m]-355./33*k2[m]+46732./5247*k3[m]+49./176*k4[m]-5103./18656*k5[m]);
    
    f(temp,k6,nu);
    
    //calc k7
    for(int q=0;q<N;q++)
      temp[q]=y[q]+dt*(35./384*k1[q]+500./1113*k3[q]+125./192*k4[q]-2187./6784*k5[q]+11./84*k6[q]);
    
    f(temp,k7,nu);
    
    
    //now temp is actualy y_n+1 of RK4
    
    //calc y of RK5
    for(int ii=0;ii<N;ii++)
      y[ii]=y[ii]+dt*(5179./57600*k1[ii]+7571./16695*k3[ii]+393./640*k4[ii]-92097./339200*k5[ii]+187./2100*k6[ii]+1./40*k7[ii]);
    
       
      //next step
    t+=dt;
    //stepsize control for the next step
    for(int jj=0;jj<N;jj++)
      chi[jj]=abs(temp[jj]-y[jj]);
    
    max=maximum(chi[0],chi[2]);
    
    //new stepsize with quality faktor of 0.9
    double dtt=dt;
    dt=0.9*dt*pow(tola/max,1./5);
    //if bigger then 10⁻3 then use old value
    if (dt>1e-3)
    dt=dtt;
    //use result of RK4 (temp) and calc next step. We can just use that y and temp are just pointer to arrays. So if the addres to the array y is exchanged with the addres of temp (y=temp), y point direktly to y_n+1 cal by RK4
    swap=y;
    
    
    y=temp;
    temp=swap;
       
    
  }
 
  out.close();
  delete[] y;
  delete[] temp;
  return 0;
}


void f(const double* const y, double* const k, const double& nu){
  //y[0]=x,y[1]=x',y[2]=y,y[3]=y'
  //parametervektor
  double* const p = new double [2];
  //r
  p[0]=sqrt(pow((y[0]+nu),2)+pow(y[2],2));
  //s
  p[1]=sqrt(pow((y[0]-1+nu),2)+pow(y[2],2));
  
  
  k[0]=y[1];
  k[1]=y[0]+2*y[3]-((1-nu)*(y[0]+nu))/(pow(p[0],3))-(nu*(y[0]-1+nu))/(pow(p[1],3));
  k[2]=y[3];
  k[3]=y[2]-2*y[1]-((1-nu)*y[2])/pow(p[0],3)-(nu*y[2])/pow(p[1],3);
  delete[] p;
}


double maximum(double const a, double const b){
  return(a<b?b:a);  
}