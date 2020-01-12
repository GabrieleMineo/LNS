#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <stdlib.h>
#include "vector"
#include "classi.h"


using namespace std;


int main(){  
 
 int N=0,M = 1E6,BLOCKS=1E2,L=M/BLOCKS;  	
 double R=0., R_ave=0., R_var=0.;		// Distanza dall'origine e valori medi
 ofstream out;	
 ofstream out_ave;

	for (unsigned int gauss = 0; gauss <= 1; gauss += 1)
	{
		for (unsigned int configurazione = 0; configurazione < 2; configurazione += 1){
			if(gauss==0){			
				if(configurazione==0) {  
				 out.open("out/1s.out");	
				 out_ave.open("out/r1s.out");
		 		}else{
				 out.open("out/2p.out");	
				 out_ave.open("out/r2p.out");
				}
			}else{
				if(configurazione==0) {  
				 out.open("out/1s_gauss.out");	
				 out_ave.open("out/r1s_gauss.out");
		 		}else{
				 out.open("out/2p_gauss.out");	
				 out_ave.open("out/r2p_gauss.out");
			}
		 }

		 R = 0.;
		 R_ave = 0.;	
		 R_var = 0.;
		 N = 0;
		 MT sampler(configurazione, gauss);		 
		 out << sampler.Getn() << " " << sampler.Getl() << " " << sampler.Getm() << " " << endl;  // intestazione del file: n l m   
		 sampler.Equilibrate();				// Fase di equilibrazione (alcuni passi del random walk)
		 for(int i = 0; i < M; i++){		
		 
			 	if(i%BLOCKS == 0){			// Stampa su file le coordinate attuali
					out << fixed << setprecision(6) << sampler.GetCoord(0) << " " << sampler.GetCoord(1)  << " " << sampler.GetCoord(2)  << endl;
			 	}			

			 	R += sampler.GetRadius()/(double)L;	// Aggiorna la media della distanza dall'origine nel blocco 
				if((i+1)%L == 0){			// medie blocco	
					N++;
					R_ave += R;
					R_var += R*R;
					out_ave << R_ave/(double)N << " " << R_var/(double)N << " " <<  pow(R_ave/(double)N, 2) << endl;
					R = 0;
				}
				sampler.Stepper();		
			}
			
			cout << " acceptance " << (double)sampler.GetCount()/(double)M << endl;

			out.close();
			out_ave.close();
		}
	}

	return 0;
}


