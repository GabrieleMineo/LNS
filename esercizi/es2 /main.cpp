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


using namespace std;


double integranda(double x);
double inversa(double x);
double taylor(double x);
vector<double> mc_integral(Random rnd, int ntot=1E6, bool importance_sampling=false);


int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;

	 vector<double> mean_results;
	 vector<double> mean_uncertainties;
   vector<double> variance_results;
	 vector<double> variance_uncertainties;
	 vector<double> somma_rv_unif, somma_rv_gauss,somma_rv_exp,somma_rv_cauchy;

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

///integrale mc
	ofstream i_file("out/i_file.dat");
	vector<double> integral = mc_integral(rnd);
	vector<double> integral_importance_sampling = mc_integral(rnd,1E6, true);
  Statistica c(integral, 1000);
	Statistica c_is(integral_importance_sampling , 1000);
  vector<double> prog_i=c.get_progressi_media();
  vector<double> prog_ei=c.get_progressi_error();
  vector<double> prog_i_is=c_is.get_progressi_media();
  vector<double> prog_ei_is=c_is.get_progressi_error();
	for ( int i = 2; i < prog_i.size(); i += 1)
		i_file << i << "\t" << prog_i[i] << "\t" << prog_ei[i]<< "\t" << prog_i_is[i] << "\t" << prog_ei_is[i] <<endl;
	i_file.close();
//importance_sampling


///////////////////////////
//walker
	ofstream w_file("out/walker_file.dat");
	ofstream w_path("out/walker_path.dat");
	ofstream wc_file("out/c_walker_file.dat");
	ofstream wc_path("out/c_walker_path.dat");
	int numero_mosse=100;
	int numero_sim=1E4;
	
	vector<int> mosse = rnd.i_random_vector(UNIFORM,numero_sim*numero_mosse,0.,6.);
	vector<double> distanze_medie;
	vector<double> distanze_errori;
 	vector<walker> walkers;
	for ( int i = 0; i < numero_sim; i += 1)
	{
		vector<int> sottomosse(mosse.begin()+i*numero_mosse, mosse.begin()+(i+1)*numero_mosse);
		walkers.push_back(walker(sottomosse));
	}

	for ( int i = 0; i < numero_mosse; i += 1){
		vector<double> distanze_attuali;
		for ( int j = 0; j < numero_sim; j += 1)	
		{
			walkers[j].aggiorna();
			distanze_attuali.push_back(walkers[j].get_distance());
		}
		w_path << i << "\t"<< walkers[0].get_x()<< "\t"<< walkers[0].get_y()<< "\t"<< walkers[0].get_z()<<endl;
		vector<double> v = distanze_attuali;
		double sum = std::accumulate(v.begin(), v.end(), 0.0);
		double mean = sum / v.size();

		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / v.size() - mean * mean)/sqrt(v.size());
		w_file << (double)1.*i << "\t" << sqrt(mean) <<"\t"<<stdev<<endl;

	}
///continuous walker


	vector<double> theta = rnd.random_vector(UNIFORM,numero_sim*numero_mosse,0.,2.*M_PI);
	vector<double> phi = rnd.random_vector(UNIFORM,numero_sim*numero_mosse,0.,M_PI);
	vector<double> c_distanze_medie;
	vector<double> c_distanze_errori;
 	vector<walker_c> c_walkers;
	for ( int i = 0; i < numero_sim; i += 1)
	{
		vector<double> sottotheta(theta.begin()+i*numero_mosse, theta.begin()+(i+1)*numero_mosse);
		vector<double> sottophi(phi.begin()+i*numero_mosse, phi.begin()+(i+1)*numero_mosse);
		c_walkers.push_back(walker_c(sottotheta,sottophi));
	}

	for ( int i = 0; i < numero_mosse; i += 1){
		vector<double> distanze_attuali;
		for ( int j = 0; j < numero_sim; j += 1)	
		{
			c_walkers[j].aggiorna();
			distanze_attuali.push_back(c_walkers[j].get_distance());
		}
		wc_path << i << "\t"<< c_walkers[0].get_x()<< "\t"<< c_walkers[0].get_y()<< "\t"<< c_walkers[0].get_z()<<endl;

		vector<double> v = distanze_attuali;
		double sum = std::accumulate(v.begin(), v.end(), 0.0);
		double mean = sum / v.size();

		double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
		double stdev = std::sqrt(sq_sum / v.size() - mean * mean)/sqrt(v.size());
		wc_file << (double)1.*i << "\t" << sqrt(mean) <<"\t"<<stdev<<endl;

	}

	w_file.close();
	wc_file.close();

   return 0;
}




//funzioni utili


double integranda(double x){
	return 0.5*M_PI*cos(0.5*M_PI*x);
}
double inversa(double x){
	return 1. - sqrt(1. - x);
}



vector<double> mc_integral(Random rnd, int ntot, bool importance_sampling){
	vector<double> x = rnd.random_vector(UNIFORM, ntot, 0.,1.);
	vector<double> I;
	for ( int i = 0; i < ntot; i += 1)
	{
		if(importance_sampling==true)I.push_back( 0.5 * integranda( inversa( x[i]  ) ) / ( 1. - inversa(x[i]) ) );
		else 	I.push_back(integranda(x[i]));
	}	
	return I;

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
