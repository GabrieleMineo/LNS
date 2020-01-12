
int Fattoriale(int n);
double Bin(int j, int k);



class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  //Random(int *);
  // destructor
  ~Random();
  // methods

  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
};


//funzione d'onda
class Psi {

	private:
		int n, m, l;

	public:
		Psi(unsigned int configurazione = 0);  //conf = 0 primo esercizio
		~Psi(){};
		int get_n(){return n;};
		int get_l(){return l;};
		int get_m(){return m;};
		void set_n(int n_){n=n_;};
		void set_l(int l_){l=l_;};
		void set_m(int m_){m=m_;};
		double Norm();		
		double Laguerre(double);	
		double Harmonic(double);	
		double Probability(double r[3]); 	// |psi|^2
};


//sampler metropolis si occuperà di gestire un po' tutta la simulazione
class MT{

	private:
		double x[3];			
		double a; //questo è il mio passo		
		int count;			// # mosse accettate
		int n_step_equilibrazione;			
		int gauss_step;

		Psi orb;	
		Random rnd;			

	public:
		MT(int configurazione, int is_gauss_step);	
		~MT(){};
		
		void Equilibrate();	
		void Stepper();		
		bool AR(double, double);	//test acc rej
		double GetCoord(int);
		double GetRadius();	

		int GetCount(){return count;};
		int Getn(){	return orb.get_n();};		
		int Getl(){	return orb.get_l();};
		int Getm(){	return orb.get_m();};
};



