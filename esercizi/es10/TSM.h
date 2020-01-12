#ifndef  TSM
#define TSM

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include "random.h"

using namespace std;

struct risultato{	//ottimo per un bel MIN_LOC :)
	double lunghezza;
	int rank_i;	
};

struct posizione{
	double x;
	double y;
};

class Path {
	private:
		int N;								//numero nodi
		vector <int> cammino;	//il path vero e proprio
		double L;							//sua lunghezza
		
	public:
		Path(){};
		~Path(){};	
		Path(int, vector<int>, vector<posizione>);
		Path(const Path& p){					//un copy-constructor per sicurezza
			N=p.N;
		  cammino = p.cammino;
			L=p.L;	
		}

		void Controllo();														//controllo validità percorso
		void LengthCalc(vector<posizione>);				//calcola lunghezza percorso
		double GetLength(vector<posizione>);				//calcola e restituisce

		int GetN(){return N;};							
		int GetElement(int i){return cammino[i];};						
		void SetElement(int i, int sost){cammino[i]=sost;}				
	
};

Path Scambio(Random *, Path );
Path T_Random(Random*, Path );		 

#endif
