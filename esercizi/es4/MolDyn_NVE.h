/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//variabili mie
int nconf = 0;			//Variabile globale che tiene conto dello step di integrazione
bool ctrl = 0;			//Variabile di controllo per la fase di equilibrazione
int restart_steps = 1; 		//Variabile che contiene il numero di step di integrazione da eseguire nella fase di equilibrazione

int Nblocks;			//Numero di blocchi
int nblock = 0;			//Variabile per le medie nei blocchi

double epot_BlockMean = 0, ekin_BlockMean= 0, etot_BlockMean = 0, temp_BlockMean = 0;
double epot_BlockVar = 0, ekin_BlockVar= 0, etot_BlockVar = 0, temp_BlockVar = 0;	

double stima_pr = 0;					//Pressione istantanea
double ave_pr = 0, pr_BlockMean = 0, pr_BlockVar = 0;	//Pressione media con incertezza statistica
//////

//funzioni mie
bool Interazione1(void);
void Interazione2(void);

void Equilibrate(void);		//Funzione che gestisce l'equilibrazione del sistema con i seguenti step:
double get_temp(void);		// - mostra la temperatura dopo n passi di integrazione
double get_fs(void);		// - in caso negativo calcola il fattore di riscalamento per le velocità
void Rescale(void);		// - riscala le velocità
void OldConfFinal(void);	//Funzione che salva su file le posizioni vecchie

//parameters, observables
const int m_props=4;
int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart, nbins;
const int gr_bin = 250;		// Bin per l'istogramma della g(r)
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);


double ave_epot = 0, ave_ekin = 0, ave_etot = 0, ave_temp = 0;				//Variabili per il calcolo  
 
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
