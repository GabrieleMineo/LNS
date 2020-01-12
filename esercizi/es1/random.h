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



double x_up(int a);
double y_up(int a);
double z_up(int a);
double x_up(double,double);
double y_up(double,double);
double z_up(double,double);

double r(double x, double y, double z);



class walker{

	private:
		double x_,y_,z_;
		vector<int> mosse;
		void aggiorna_n_volte(int nsteps){ for (unsigned int i = 0; i < nsteps; i += 1) aggiorna();		}
		int step;

	public:
		walker(){
			x_=0.;
			y_=0.;
			z_=0.;
		};
		walker(vector<int> rand_m){
			x_=0.;
			y_=0.;
			z_=0.;
			step=0;
			mosse=rand_m;
		};
  	~walker(){};	

		void aggiorna(){
			x_+=x_up(mosse[step]);
			y_+=y_up(mosse[step]);
			z_+=z_up(mosse[step]);
			step++;
		}
		double get_distance(){
			return r(x_,y_,z_);
		}
		double get_final_distance(){
			aggiorna_n_volte(mosse.size());			
			return r(x_,y_,z_);
		}	

		vector<double> get_vec_distance(){
			vector<double> v;
			for (unsigned int i = 0; i < mosse.size(); i += 1){
				aggiorna();				
				v.push_back(get_distance());
			}			
			return v;
		}

		double get_x(){return x_;}
		double get_y(){return y_;}
		double get_z(){return z_;}
		int get_mosse(int i){return mosse[i];}
};


class walker_c{

	private:
		double x_,y_,z_;
		vector<double> theta, phi;
		int step;

	public:
		walker_c(){
			x_=0.;
			y_=0.;
			z_=0.;
		};
		walker_c(vector<double> theta_in, vector<double> phi_in){
			x_=0.;
			y_=0.;
			z_=0.;
			step=0;
			theta=theta_in;
			phi=phi_in;
		};
  	~walker_c(){};	

		void aggiorna(){
			x_+=x_up(theta[step],phi[step]);
			y_+=y_up(theta[step],phi[step]);
			z_+=z_up(theta[step],phi[step]);
			step++;
		}
		double get_distance(){
			return r(x_,y_,z_);
		}

		vector<double> get_vec_distance(){
			vector<double> v;
			for (unsigned int i = 0; i < theta.size(); i += 1){
				aggiorna();				
				v.push_back(get_distance());
			}			
			return v;
		}

		double get_x(){return x_;}
		double get_y(){return y_;}
		double get_z(){return z_;}
};



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
