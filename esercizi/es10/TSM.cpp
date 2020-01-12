
#include "TSM.h"
using namespace std;

Path :: Path(int N_cities, vector <int> path, vector <posizione> nodi_) {
	N=N_cities;		
	cammino=path;		
	LengthCalc(nodi_);
}

void Path :: LengthCalc(vector <posizione> nodi_) {
	double l_accumulata=0.;
	for( int i = 0; i < N-1; i++ ) {
		int j = (i+1)%N; //forzo le pbc
		l_accumulata += sqrt( (nodi_[cammino[i]].x-nodi_[cammino[j]].x) * (nodi_[cammino[i]].x-nodi_[cammino[j]].x)+
                          (nodi_[cammino[i]].y-nodi_[cammino[j]].y) * (nodi_[cammino[i]].y-nodi_[cammino[j]].y));
	}
	L=l_accumulata;
}

double Path ::	GetLength(vector <posizione> nodi_) {
	LengthCalc(nodi_);
	return L;
}


//Controlla che le città ci siano tutte e una volta sola.
void Path :: Controllo() {
	int indice=0;
	int controllo=0;
	int appoggio;
	for(indice=0; indice<N; ++indice) {
		controllo=0;	//resetta il controllo
		for(int i=0; i<N; ++i) {
			appoggio=cammino[i];	//che indice è?
			if(controllo==1 and appoggio==indice) {
				cout<<"ripetizione"<<endl;
				return ;
			}
			if(appoggio==indice) {
				controllo=1; 	//trovato
			}
		}
		if(controllo==0) {
			cout<<"nodo non trovato"<<endl;
		}
	}
}	



Path Scambio(Random * rnd, Path sequence) {
	Path seqC=sequence;					
	int quale1=(int)floor(rnd->Rannyu(0,sequence.GetN()));		//sceglie il primo indice
	int quale2=(int)floor(rnd->Rannyu(0,sequence.GetN()));		//sceglie il secondo indice
	int app=seqC.GetElement(quale1);			//permuta
	seqC.SetElement(quale1, seqC.GetElement(quale2));	//...
	seqC.SetElement(quale2, app);				//...
	seqC.Controllo();					//controlla che sia tutto ok
	return seqC;						//ritorna la sequenza permutata
}


		
                





