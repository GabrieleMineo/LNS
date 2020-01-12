/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"
#include <random>

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exponential(double lambda){
	double u = Rannyu();
	double x = -log(u)/lambda;
	return x;

}

double Random :: Cauchy(double x0, double gamma){
	double u = Rannyu();
	return x0 + gamma*tan(M_PI*(u-0.5));
}
double Random :: Gauss(){
	return Gauss(0.,1.);
}
double Random :: Exponential(){
	return Exponential(1.);
}
double Random :: Cauchy(){
	return Cauchy(0.,1.);
}

double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

vector<double> Random :: random_vector(int dist_type, int numero_elementi, double primo_parametro, double secondo_parametro){
	vector<double> v;
	double rv=0.;	
	for (unsigned int i = 0; i < numero_elementi; i += 1){
			if(dist_type==UNIFORM) rv=Rannyu(primo_parametro,secondo_parametro);
			else if(dist_type==GAUSS) rv=Gauss(primo_parametro,secondo_parametro);
			else if(dist_type==EXPONENTIAL) rv=Exponential(secondo_parametro);
			else  rv=Cauchy(primo_parametro,secondo_parametro);
		v.push_back(rv);
	}
	return v;

}
vector<int> Random :: i_random_vector(int dist_type, int numero_elementi, double primo_parametro, double secondo_parametro){
	vector<int> v;
	int rv=0;	
	for (unsigned int i = 0; i < numero_elementi; i += 1){
			if(dist_type==UNIFORM) rv = (int)floor(Rannyu(primo_parametro,secondo_parametro));
			else if(dist_type==GAUSS) rv = (int)floor( Gauss(primo_parametro,secondo_parametro) );
			else if(dist_type==EXPONENTIAL) rv = (int)floor( Exponential(secondo_parametro));
			else  rv = (int)floor( Cauchy(primo_parametro,secondo_parametro) );
		v.push_back(rv);
	}
	return v;

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

void Statistica :: get_mean_stdev(double& media_out, double& devstd_out){

		if(numero_blocchi==1){
			vector<double> v = v_m;
			double sum = std::accumulate(v.begin(), v.end(), 0.0);
			media = sum / v.size();

			double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
			stdev = std::sqrt(sq_sum / v.size() - media * media);
		}else{
			
			vec_media_devstd();
			media=prog_media.back();
		 	stdev=prog_error.back();
		}
		media_out = media;
		devstd_out = stdev;

}


void Statistica :: vec_media_devstd(){

	vector<double> av,av2;	

	for (unsigned int i = 0; i < numero_blocchi; i += 1)
	{
		double varianza_accu=0.,media_accu=0.;
		for (unsigned int j = 0; j < larghezza_blocco; j += 1)
		{
			media_accu += v_m[j+i*larghezza_blocco];
		}
		media_accu /= (double)larghezza_blocco;
		av.push_back(media_accu);
		av2.push_back(media_accu*media_accu);

	}

  for (unsigned int i = 0; i < numero_blocchi; i += 1)
  {
		double prog_m=0., prog_q=0.;
  	for (unsigned int j = 0; j < i+1; j += 1)
  	{
  		prog_m+=av[j];
			prog_q+=av2[j];			
  	}
		prog_m/=(double)(i+1.);
    prog_q/=(double)(i+1.);
		prog_media.push_back(prog_m);
	  prog_quadr.push_back(prog_q);
		prog_error.push_back(sqrt( (prog_q - prog_m*prog_m) /(double)i ));
  }

}




/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
