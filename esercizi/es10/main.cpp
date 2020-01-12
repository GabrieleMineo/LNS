#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TSM.h"
#include "random.h"
#include "mpi.h"

using namespace std;

const int N_TOT=1E3;
const double T_target=1E-3;

void read_input(int& configurazione,int& N,int& Nmosse, double& temp, vector<double>& t_i);
void generate_input(Random* rnd, vector<posizione>& cities, vector<int>& sequence, int quadrato_cerchio);

int main(int argc, char* argv[]) {

	MPI::Init(argc,argv);		
	
	int size = MPI::COMM_WORLD.Get_size();	
	int rank = MPI::COMM_WORLD.Get_rank();	
//Generatore random.
	Random rnd;		//crea oggetto Random.
   	int seed[4];
	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
   	   for(int i=0; i<rank+1; i++) {
   	   	Primes >> p1 >> p2 ;
   	   }
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
	   

	////////////////////////////////////////////

	risultato r; //qui storerò una sorta di mappa lunghezza/rank con quella lunghezza, mi serve per il minloc di mpi
	ofstream bestLength, finalConfig, CitiesPosition, BestConf, Data;
	int N, configurazione, N_mosse;
	double temp, phi, x, y, L, Lp;

	vector<double> temp_ladder;
	read_input(configurazione, N, N_mosse, temp, temp_ladder);

	vector<posizione> cities(N);	//punti
	vector<int> sequence(N);	//ordine visita.


	BestConf.open("out/conf_best.dat"); //qua ci scrive rank 0 riga per riga i vettori di sequenza di visita, mi serve per animare
	if(rank==0) {
		bestLength.open("out/lunghezzeBest.dat");
		Data.open("out/progressi.dat");
		generate_input(&rnd, cities, sequence, configurazione);		//genera l'input
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0; i<cities.size(); i++) {
			MPI_Bcast(&(cities.at(i).x),1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(cities.at(i).y),1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&(sequence.at(i)),1, MPI_INT, 0, MPI_COMM_WORLD);
		}
	MPI_Barrier(MPI_COMM_WORLD);
	
	Path p(N,sequence,cities);
	int lanci=0, accettati=0;		//Per tassi di accettazione
	risultato out;
	
	for (unsigned int i = 0; i < N_TOT; i += 1){
		accettati=0.;
		for (unsigned int j = 0; j < N_mosse; j += 1){
			L=p.GetLength(cities);
			Path proposta=Scambio(&rnd, p);	
			Lp=proposta.GetLength(cities);		
			double prob=min(1., exp((L-Lp)/temp_ladder[i]));		

			if(rnd.Rannyu(0.,1.) < prob){
				p=proposta;
				accettati++;
			}
		}
	//ad ogni temperatura il migliore stampa la configurazione
		r.lunghezza=p.GetLength(cities);
		r.rank_i=rank;
		if(rank==0){
			Data << temp_ladder[i] << "\t" << r.lunghezza << "\t" <<(double)accettati/(double)N_mosse <<endl;
			for (unsigned int k = 0; k < p.GetN(); k += 1) BestConf<< p.GetElement(k)<<"\t";
			BestConf<<endl;
			
		}
  
	}

	MPI_Barrier(MPI_COMM_WORLD);

	L=p.GetLength(cities);	
	r.lunghezza=L;
	r.rank_i=rank;

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce( &r, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD );
	MPI_Bcast(&out.rank_i,1, MPI_INT, 0, MPI_COMM_WORLD);

	cout<<" rank " << r.rank_i <<" has length:" << r.lunghezza <<endl;
	
	//BEST RANK
	if(rank==out.rank_i) {
		cout<<"best rank is " << rank <<" with length: " << L << endl;
		for(int i = 0; i < p.GetN(); i++) BestConf<< p.GetElement(i)<<"\t";
		BestConf << endl;	

		bestLength<< L << endl;
		bestLength.close();
	}
	
	BestConf.close();
	MPI::Finalize();
	return 0;
}

	

/////////////////funzioni

void read_input(int& configurazione,int& N,int& Nmosse, double& temp, vector<double>& t_i){
	double alpha=0.99;
	ifstream dataIn("in.dat");	

	dataIn>>configurazione;		//0 cerchio 1 quadrato
	dataIn>>N;			//Numero città
	dataIn>>temp;			//Temperatura iniziale
	dataIn>>Nmosse;			//Numero di mosse per ogni Temperatura
	dataIn>>alpha;
	dataIn.close();

	t_i.push_back(temp);
	for (unsigned int i = 0; i < N_TOT; i += 1) t_i.push_back(t_i.back()*alpha);
}

void generate_input(Random* rnd, vector<posizione>& cities, vector<int>& sequence, int quadrato_cerchio){
	ofstream CitiesPosition("out/posizioni.dat");
	for(int i=0; i<cities.size(); ++i) {
		if(quadrato_cerchio==1){		
			double phi=rnd->Rannyu(0.,2*M_PI);
			cities.at(i).x=cos(phi);
			cities.at(i).y=sin(phi);
		}else{
			double x=rnd->Rannyu(-0.5,0.5);
			double y=rnd->Rannyu(-0.5,0.5);
			cities.at(i).x=x;
			cities.at(i).y=y;
		}
		sequence.at(i)=i;	
		CitiesPosition<< cities.at(i).x <<"\t" << cities.at(i).y << endl;

	}
	CitiesPosition.close();
}



