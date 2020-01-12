

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iterator>
#include "random.h"


using namespace std;
vector<double> pricer_vanilla(Random* rnd, bool Call, bool eulero, 
															double S, double T, double t, double s, 
															double r, int nstep, double K);

vector<double> pricer_vanilla(vector<double> v, bool Call=true, bool eulero=false, 
															double S=100., double T=1., double t=0., double s=0.25, 
															double r=0.1, int nstep=100, double K = 100.);
double euro_pv(double S, double K, bool call);


const int N_PATH=1E6;

int main (int argc, char *argv[]){

 	 ofstream third("third.dat");

   Random rnd;
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
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;


	double S=100., T=1., t=0., s=0.25, r=0.1, K=100.;
	bool Call=true, eulero=false;
	int nstep=10;
	


//	vector<double> w = rnd.random_vector(GAUSS, N_PATH, 0., 1.);
	vector<double> p_call_eulero    = pricer_vanilla(&rnd, Call, true, S, T, t, s, r, nstep, K);
	vector<double> p_call_analitica = pricer_vanilla(&rnd, Call, eulero, S, T, t, s, r, nstep, K);
	vector<double> p_put_analitica  = pricer_vanilla(&rnd, false, false, S, T, t, s, r, nstep, K);
	vector<double> p_put_eulero     = pricer_vanilla(&rnd, false, true, S, T, t, s, r, nstep, K);

	for (int i = 0; i < N_PATH; i += 1) 
			third << p_call_analitica[i] << "\t"<< p_call_eulero[i] << "\t"<< p_put_analitica[i] << "\t"<<p_put_eulero[i] << "\t"<< endl;
	
	 third.close();	

   return 0;
}



vector<double> pricer_vanilla(vector<double> w, bool Call, bool eulero, 
															double S, double T, double t, double s, 
															double r, int nstep, double K){
		int Nsim=w.size();
		double S_T=S;		
		vector<double> payoffs;

		for (int j = 0; j < w.size(); j++){
			S_T=S*exp(   (r-0.5*s*s)*T + s*sqrt(T)*w[j] ); 
			payoffs.push_back(exp(-r)*euro_pv( S_T, K, Call ));   
		}
		return payoffs;
				
}

vector<double> pricer_vanilla(Random* rnd, bool Call, bool eulero, 
															double S, double T, double t, double s, 
															double r, int nstep, double K){


		double S_T; 
		double dt= T/(double)nstep;		
		cout <<dt<<endl;
		vector<double> payoffs;

			for (int j = 0; j < N_PATH; j++){
				S_T=S;
				if(eulero==false) S_T*=exp(   (r-0.5*s*s)*T + s*sqrt(T)*rnd->Gauss() ); 
				else for (unsigned int i = 0; i < nstep; i += 1) S_T *= exp(  (r - 0.5*s*s)*dt + s * sqrt(dt) * rnd->Gauss()  );
				payoffs.push_back(exp(-r)*euro_pv( S_T, K, Call ));   
			}
			
		return payoffs;
				
}

double euro_pv(double S, double K, bool call){
	double pay;
	if(call) pay=max(S-K,0.);
	else pay=max(K-S,0.);
	return pay;
}


