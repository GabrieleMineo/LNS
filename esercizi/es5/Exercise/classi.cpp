

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "classi.h"

using namespace std;


int Fattoriale(int n){
    if (n == 0)
       return 1;
    return n * Fattoriale(n - 1);
}

//coeff newton
double Bin(int j, int k){
	return (double)Fattoriale(j)/(double)(Fattoriale(k)*Fattoriale(j-k));
}


 
Random :: Random(){
 int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();
 
 ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
return;
}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

/////////////////////////////////////////////////////////////////////

MT :: MT (int configurazione = 0, int if_gauss=0){
	count = 0;	
	n_step_equilibrazione = 50; 
	gauss_step=if_gauss;
	if(configurazione==0){
		orb.set_n(1);	
		orb.set_l(0);	
		orb.set_m(0);	
		a = 1.;	
	}else{
		orb.set_n(2);	
		orb.set_l(1);	
		orb.set_m(0);	
		a = 2.8;
	}
	x[0]=0.3;
	x[1]=-0.3;
	x[2]=0.1;
}


void MT :: Equilibrate(){
	for(int i = 0; i < n_step_equilibrazione; i++) Stepper();
	count = 0;
}

void MT :: Stepper(){
	double p, pnew;
	double x_try[3]; 
	double x_inc[3];	
	if(gauss_step==1){				//incrementi gaussiani
		x_inc[0] = rnd.Gauss(0., a/sqrt(3.));
		x_inc[1] = rnd.Gauss(0., a/sqrt(3.));
		x_inc[2] = rnd.Gauss(0., a/sqrt(3.));
	}else{									//uniforme
		double theta = acos(1.-2.*rnd.Rannyu());		
		double phi = rnd.Rannyu(0., 2.*M_PI);
		x_inc[0] = a*sin(theta)*cos(phi);
		x_inc[1] = a*sin(theta)*sin(phi);
		x_inc[2] = a*cos(theta);
	}

	p = orb.Probability(x);			// attuale				

	for (int i = 0; i < 3; i ++) x_try[i] = x[i] + x_inc[i] ;		
	pnew = orb.Probability(x_try);			// prob della nuova configurazione

	if (AR(pnew, p)){			// accept/reject
		for (int i = 0; i < 3; i ++) x[i] = x_try[i];			 
  	count++;
	}

}

bool MT :: AR(double x, double y){
	bool a;	
	if( x > y || rnd.Rannyu() < x/y  ) a = true;
	else a = false;
	return a;
}

double MT :: GetCoord(int n){	
	if(n<3){	
		return x[n];
	}else return -6*1e6;
}

double MT :: GetRadius(void){
 double r = 0.;
 for(int i = 0; i < 3; i++) r += x[i]*x[i];
 return sqrt(r);
}

Psi :: Psi (unsigned int primo_o_secondo){
	if (primo_o_secondo==0){
		n=1;
		l=0;
		m=0;
	}
	else{
		n=2;
		l=1;
		m=0;
	}
}

double Psi :: Norm(void){
 int j = n - l - 1;
 return pow(2.0/(double)n, 3) * Fattoriale(j)/(double)(2.0 * n * Fattoriale(n + l));
}

double Psi :: Laguerre(double x){
	int j = n - l - 1;
	int a = 2*l + 1;
	double y = 0.;
	for(int i = 0; i < j+1; i++){
		y += pow((-1), i) * Bin(j + a, j - i) * pow(x, i) / (double)Fattoriale(i); 
	}

	return y*y;
}

double Psi :: Harmonic(double x){
	double h=0.;	
	if(l == 0) h=0.25/M_PI;	
	else       h=0.25 * 3./M_PI * x * x;
	return h;		
}

// probabilità come modulo quadro della funzione 
double Psi :: Probability(double r[3]){
 double R = 0;
 for(int i = 0; i < 3; i++) R += r[i]*r[i];
 R = sqrt(R);
 double cos_theta = r[2]/R;
 return Norm() * Laguerre(2*R/(double)n) * Harmonic(cos_theta) *  pow((2*R/(double)n), 2*l) * exp(- 2*R/(double)n);
}

