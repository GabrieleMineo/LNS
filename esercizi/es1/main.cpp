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
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <iterator>
#include "random.h"

const int M=1E6;		//Numero totale punti generati
const int N_BLOCK=100;		//Numero dei blocchi
const int length_block=M/N_BLOCK;		//Numero dei punti in un blocco

using namespace std;


double chi_squared_test(int M, Random* rnd);
void rv_sum_vector(Random rnd, vector<double>& rvs, int N_somma, int N_realizzazioni,  int dist_type);
vector<double>  buffon_experiment(Random rnd, int ntot=1E5, int blocks=2);
double integranda(double x);
double inversa(double x);
double taylor(double x);
vector<double> mc_integral(Random rnd, int ntot=1E5, bool importance_sampling=false);


int main (int argc, char *argv[]){

 	 ofstream first_prog("out/first_prog.dat");
	 ofstream second_prog("out/second_prog.dat");

	 //chi squared file
	 ofstream chi("out/chi.dat");
	 //central limit theorem files
	 ofstream uniff("out/unif_clt.txt");
	 ofstream gaussf("out/gauss_clt.txt");
	 ofstream exponentialf("out/exponential_clt.txt");
	 ofstream cauchyf("out/cauchy_clt.txt");


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



   rnd.SaveSeed();

//////////////////////	
//Esercizi 1.01 e 1.02
	 int M = 1E4;
	 double ave_1=0., var_1=0. , ave_2=0., var_2=0.;
//genero un vettore di numeri random disbruiti uniformemente

	 vector<double> xv=rnd.random_vector(UNIFORM,M);
	 for (unsigned int i = 0; i < N_BLOCK; i += 1){
		double media=0.;
		double var=0.;
	 	for (unsigned int j = 0; j < length_block; j += 1){
	 		media+=xv[i]/length_block;
			var+= (xv[i]-0.5)*(xv[i]-0.5)/length_block;
	 	}
		ave_1 += media;
		var_1 += media*media;
		ave_2 += var;
		var_2 += var*var;
		first_prog << ave_1/(i+1) << "\t" << var_1/(i+1) <<  "\t"  << pow(ave_1/(i+1),2) <<endl;
		second_prog<< ave_2/(i+1) << "\t" << var_2/(i+1) <<  "\t" <<  pow(ave_2/(i+1),2)<<endl;		

	 }

	first_prog.close();
	second_prog.close();
	 

///////////////////////
// Esercizio chi-quadro
	 for (unsigned int i = 1; i < 100; i += 1)
	 		chi << i <<"\t"<<chi_squared_test(i, &rnd)<< endl;
//////////////////////	


//Esercizio Teorema Limite Centrale

	 vector<double> somma_rv_unif, somma_rv_gauss, somma_rv_exp, somma_rv_cauchy;
	 int N_somma=100, N_realizzazioni=1E4;
	 rv_sum_vector(rnd, somma_rv_unif, N_somma, N_realizzazioni,  UNIFORM);
	 rv_sum_vector(rnd, somma_rv_gauss, N_somma, N_realizzazioni,  GAUSS);
	 rv_sum_vector(rnd, somma_rv_exp, N_somma, N_realizzazioni,  EXPONENTIAL);
	 rv_sum_vector(rnd, somma_rv_cauchy, N_somma, N_realizzazioni,  CAUCHY);
	
	for (unsigned int i = 0; i < N_realizzazioni; i += 1)
	{
		uniff << somma_rv_unif[i] <<endl;
		gaussf << somma_rv_gauss[i] <<endl;
		exponentialf << somma_rv_exp[i] <<endl;
		cauchyf << somma_rv_cauchy[i] <<endl;
	}
	
	chi.close();
	uniff.close();
	gaussf.close();
	exponentialf.close();
	cauchyf.close();	


//esercizio ago di buffon
	ofstream pi_file("out/pi_file.dat");
	vector<double> buffon = buffon_experiment(rnd);
  Statistica b(buffon, 100);
 
  vector<double> prog_pi=b.get_progressi_media();
  vector<double> prog_epi=b.get_progressi_error();
	for (unsigned int i = 2; i < prog_pi.size(); i += 1)
		pi_file << i << "\t" << prog_pi[i] << "\t" << prog_epi[i] <<endl;
	pi_file.close();
////////////////////////////

   return 0;
}




//funzioni utili


double integranda(double x){
	return 0.5*M_PI*cos(0.5*M_PI*x);
}
double inversa(double x){


	return 2./M_PI * asin(x);

// return 2.* ( pow((sqrt(9.* CO*CO*x*x - 8.) - 3.*CO* x),(1./3.)) /M_PI + 2./(M_PI * pow((sqrt(9.*CO*CO*x*x - 8.) - 3.*CO* x),( 1./3.)  )  ) );
}
double taylor(double x){
	double CO = 1.;
	return 0.5*M_PI*(1.-M_PI*M_PI*x*x/8.)/CO;
}


vector<double> mc_integral(Random rnd, int ntot, bool importance_sampling){
	vector<double> x = rnd.random_vector(UNIFORM, ntot, 0.,1.);
	vector<double> I;
	for (unsigned int i = 0; i < ntot; i += 1)
	{
		if(importance_sampling==true)I.push_back( integranda( inversa( x[i]  ) ) );
		else 	I.push_back(integranda(x[i]));
	}	
	return I;

}



vector<double> buffon_experiment(Random rnd, int ntot, int blocks ){

		vector<double> centers = rnd.random_vector(UNIFORM,ntot,0.,1.);
		vector<double> angles = rnd.random_vector(UNIFORM,ntot,0.,0.5*M_PI);
		vector<double> results;
		double l=0.99;
		int nhit=0;

		for (unsigned int i = 0; i < ntot; i += 1)
		{
				double y=0.5*l*sin(angles[i]);
				if(centers[i]<=y) nhit++;
				results.push_back(l*i/nhit);
		}
				    
		return results;

}


double chi_squared_test(int M, Random* rnd){

	int nsim = 1E4;
	double chi=0.;
	double m = 1./(double)M;
	double expected = m*(double)nsim;
	for (unsigned int i = 0; i < M; i += 1)
	{
		double sx = i*m, dx = (i+1.)*m;
		int observed = 0;
		for (unsigned int j = 0; j < nsim; j += 1)
		{
			double ni = rnd->Rannyu(0.,1.);
			if(ni<=dx && ni>=sx) observed++;
		}
		chi+=(observed - expected)*(observed - expected) / expected;

	}

	return chi;
}


void rv_sum_vector(Random rnd, vector<double>& rvs, int N_somma, int N_realizzazioni,  int dist_type){
	
	double rv_somma, rv;
	for (unsigned int j = 0; j < N_realizzazioni; j += 1)
	{
		rv_somma = 0.;
		for (unsigned int i = 0; i < N_somma; i += 1)
		{
			if(dist_type==UNIFORM) rv=rnd.Rannyu();
			else if(dist_type==GAUSS) rv=rnd.Gauss();
			else if(dist_type==EXPONENTIAL) rv=rnd.Exponential();
			else  rv=rnd.Cauchy();
			rv_somma+=rv;
		}
		rvs.push_back(rv_somma/(double)N_somma);
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
