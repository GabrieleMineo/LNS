/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


#include<vector>
#include <numeric>
#include <random>
#define UNIFORM 0 
#define GAUSS 1
#define EXPONENTIAL 2
#define CAUCHY 3






#ifndef __Random__
#define __Random__



using namespace std;



class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
	Random(const Random &a){
		m1 = a.m1;
		m2 = a.m2;
		m3 = a.m3;
		m4 = a.m4;
		l1 = a.l1;
		l2 = a.l2;
		l3 = a.l3;
		l4 = a.l4;
		n1 = a.n1;
		n2 = a.n2;
		n3 = a.n3;
		n4 = a.n4;	

	}
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
	double Exponential(double lambda);
	double Cauchy(double x0, double gamma);
  double Gauss();
	double Exponential();
	double Cauchy();
	vector<double> random_vector(int tipo_distribuzione=0, int numero_elementi=0, double primo_parametro=0., double secondo_parametro=1.);
	vector<int> i_random_vector(int tipo_distribuzione=0, int numero_elementi=0, double primo_parametro=0., double secondo_parametro=1.);

};



class Statistica{

 private:
	double media, stdev;
	int numero_simulazioni, numero_blocchi, larghezza_blocco;
	vector<double> v_m, prog_media, prog_error, prog_quadr;
	void media_devstd();
	void vec_media_devstd();

 public:
	Statistica();
  ~Statistica(){};

	void get_mean_stdev(double& media_out, double& devstd_out);
	void set_vec(vector<double> v_p){v_m=v_p;}
	void set_num(int nblocks=1){
		numero_simulazioni=v_m.size();
		numero_blocchi=nblocks;
		larghezza_blocco = numero_simulazioni/numero_blocchi;

	}		

	Statistica(vector<double> v_p, int nblocks=1){
		v_m=v_p;
		set_num(nblocks);		
	}	

  vector<double> get_progressi_media(){
  	vec_media_devstd();		
		return prog_media;
	};

  vector<double> get_progressi_error(){
  	vec_media_devstd();				
		return prog_error;
	};

};



class path{

	private:
		double S0,St, T_maturity, t_attuale, vol, r, dt_eul;
		int step;
		int nstep_eulero;	
		Random g;	
		std::mt19937 gen;

	public:
		path(){};
 		~path(){};

		path(Random rnd, double S=100., double T=1., double t=0., double s=0.25, double r=0.1, int nstep=100){
			g = rnd;
			S0=S;
			St=S;	
			T_maturity=T;
			t_attuale=t;
			vol=s;
			dt_eul=(T_maturity-t_attuale)/nstep;
			std::random_device rd;
			std::mt19937 gen(rd()); // Create and seed the generator
			 
		}
		double Esatto(){
			std::normal_distribution<> d(0.,1.);			
			return S0*exp( (r-0.5*vol*vol)*T_maturity + vol*sqrt(T_maturity)*d(gen) );
		}
		double Eulero(){
			std::normal_distribution<> d(0.,1.);
			for (int j = 0; j < nstep_eulero; j++) St*=exp((r-0.5*vol*vol)*dt_eul+ vol*sqrt(dt_eul)*d(gen) );
			return St;
		}
		
		vector<double> get_prezzi_finali(int nsim, bool eulero=false){
			vector<double>prezzi;
			double p;
			for (int j = 0; j < nsim; j++){ 
				if(eulero==false) p = Esatto();
				else p = Eulero(); 
				prezzi.push_back(p);
			}
			return prezzi;
		}

		vector<double> get_gaussiana(){
			std::normal_distribution<> d(0.,1.);			
			vector<double> gaussiana;
			for(int j = 0; 100000 < j; j++) gaussiana.push_back(d(gen));
			return gaussiana;

		}

};





vector<double> pricer_vanilla(Random *rnd, bool Call=true, int Nsim = 1E4, bool eulero=false, 
															double S=100., double T=1., double t=0., double s=0.25, 
															double r=0.1, int nstep=100, double K = 100.);


#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
