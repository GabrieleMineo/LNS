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
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
	Input(); //Inizialization
	for(int j = 0; j < nbins; j++){
	  temp +=-1.5/(double)nbins; //T ladder
		beta=1./temp;
		Start_From_Last_Config();
	 
		if (j == 0 && is_termalizzazione==1) Thermalization(); 

		for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
		{
		  Reset(iblk);   //Reset block averages
		  for(int istep=1; istep <= nstep; ++istep)
		  {
		    Move(metro);
		    Measure();
		    Accumulate(); //Update block averages
		  }
		  Averages(iblk);   //Update results for current block
		}
		Output_T();	// risultati finali per la temperatura studiata
		cout <<" T = " << temp <<  " acceptance rate = " << accepted/attempted <<endl;

	}
  ConfFinal(); //Write final configuration

  return 0;
}



void Input()
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.out");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");
  double t;
  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;
	//MIA AGGIUNTA PER TERMALIZZARE	
	ReadInput >> is_termalizzazione;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i){
		if(rnd.Rannyu() >= 0.5) s[i] = 1;
	  else s[i] = -1;
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;

 }




void Move(int metro)
{
  int o;
  double p, p_old,p_new, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    if(metro==1) //Metropolis
    {
// INCLUDE YOUR CODE HERE
			sm = s[o];
			p_old = exp( - Boltzmann(sm, o)/temp); //calcolo energia vecchia conf
			p_new = exp( - Boltzmann(-sm, o)/temp);//calcolo energia  conf con flip
			if((p_new > p_old) || (rnd.Rannyu() < p_new/p_old)){ //mt test
				s[o] = -s[o];				//flippo
				accepted++;
			}
			attempted++;
    }
    else //Gibbs sampling
    {
// INCLUDE YOUR CODE HERE
			double flip = 0.;
			if (rnd.Rannyu() >= 0.5) flip = +1.;
			else flip = -1.;
			p = 1./(1. + exp(-2.*beta*flip*(J*(s[Pbc(o-1)] + s[Pbc(o+1)]) + h) ));	
			if (rnd.Rannyu() < p) s[o] = flip;
			else s[o] = -flip; 
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
// INCLUDE YOUR CODE HERE
     m += s[i];
  }
  walker[iu] = u;
// INCLUDE YOUR CODE HERE
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{
   
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk)  
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    //cout << "Block number " << iblk << endl;
    //if(accepted/attempted < 1e-2 ) cout << "Acceptance rate " << accepted/attempted << "  for T = "<< temp << endl << endl;    //stampa solo quando diventa molto piccolo
    Ene.open("out/output.ene.0",ios::app);
    Heat.open("out/output.Heat.0",ios::app); //Heat capacity
    Mag.open("out/output.Mag.0",ios::app);
    Chi.open("out/output.Chi.0",ios::app);

    stima_u = (blk_av[iu]/blk_norm)/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);

// INCLUDE YOUR CODE HERE

    stima_c = beta*beta*(blk_av[ic]/blk_norm - (blk_av[iu]/blk_norm)*(blk_av[iu]/blk_norm))/(double)nspin; 
    stima_m = (blk_av[im]/blk_norm)/(double)nspin;
    stima_x = beta*(blk_av[ix]/blk_norm - (blk_av[im]/blk_norm)*(blk_av[im]/blk_norm)  )/(double)nspin; 

    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);

		  Ene << iblk <<  "\t" << stima_u << "\t" << glob_av[iu]/(double)iblk << "\t" << err_u << endl;
		  Heat << iblk <<  "\t" << stima_c << "\t" << glob_av[ic]/(double)iblk << "\t" << err_c << endl;
		  Mag << iblk <<  "\t" << stima_m << "\t" << glob_av[im]/(double)iblk << "\t" << err_m << endl;
		  Chi << iblk <<  "\t" << stima_x << "\t" << glob_av[ix]/(double)iblk << "\t" << err_x << endl;
    Chi.close();
    Ene.close();
    Mag.close();
    Heat.close();

    
}




void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}


void Output_T(){
			const int wd=12;
			ofstream Ene, Heat, Mag, Chi, acc;
			Ene.open("out/mt_E.dat",ios::app);
		  Heat.open("out/mt_c.dat",ios::app); //Heat capacity
		  Mag.open("out/mt_m.dat",ios::app);
		  Chi.open("out/mt_chi.dat",ios::app);
		  acc.open("out/mt_acceptance.dat",ios::app);
			Ene << "\t" << temp <<  "\t" << stima_u << "\t" << glob_av[iu]/(double)nblk << "\t" << err_u << endl;
		  Heat << "\t" << temp <<  "\t" << stima_c << "\t" << glob_av[ic]/(double)nblk << "\t" << err_c << endl;
		  Mag << "\t" << temp <<  "\t" << stima_m << "\t" << glob_av[im]/(double)nblk << "\t" << err_m << endl;
		  Chi << "\t" << temp <<  "\t" << stima_x << "\t" << glob_av[ix]/(double)nblk << "\t" << err_x << endl;
		  acc << "\t" << temp <<  "\t" << accepted/attempted << endl;

			Ene.close();
		  Chi.close();
		  Heat.close();
		  Mag.close();

}
void Thermalization(){ 
	 int term_steps=1E5;
	 for(int k = 0; k < term_steps; k++){
			Move(metro);
	 }
}

void Start_From_Last_Config(){
  	ifstream input;
  	input.open("config.final");
    if(input) for (int i=0; i<nspin; ++i)  input >> s[i];
		Measure();
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
