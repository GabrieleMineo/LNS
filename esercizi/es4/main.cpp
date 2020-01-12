/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h> 
#include <stdio.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"
#include "random.h"


using namespace std;


int main(){ 
  

  Input();             //Inizialization
  if(ctrl) Equilibrate();	//Equilibration
  for(int istep=1; istep <= nstep; ++istep){
     Measure();     //Measure properties and save to file (not every time!)
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0){ 
     cout << "Number of time-steps: " << istep << endl;
     cout << "The actual temperature is: " << stima_temp << endl;
     }
     nconf += 1;
  }
  ConfFinal();         //Write final configuration to restart
  OldConfFinal();         //Write final configuration to restart
  
 return 0;
}




void Input(){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;
  cout << endl;
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;
  cout << "Temperature = " << temp << endl;
  
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> Nblocks;
  ReadInput >> nbins;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks for blocking = " << Nblocks << endl;
  cout << "Number of bins for g(r) = " << nbins << endl << endl;
  ReadInput.close();
  
 
  // BLOCCO IMPORTANTE // qui l'utente agisce a schermo sulle modalità di simulazione
  bool load;
  load = Interazione1(); 	// Possibilità di caricare posizioni vecchie
  Interazione2();		// Possibilità di iniziare con una fase di equilibrazione 
  if(load) return ;	// Se si sono caricate posizioni vecchie, esce dalla funzione Input

  Random rnd;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities from a Maxwell-Boltzmann distribution " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Gauss(0, sqrt(temp));
     vy[i] = rnd.Gauss(0, sqrt(temp));
     vz[i] = rnd.Gauss(0, sqrt(temp));

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];     
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;      

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
  
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;            

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}




void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}


double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}




void Measure(){ //Properties measurement
 
  double N = nstep/(double)Nblocks;	// Numero di step di integrazione in ogni blocco
  double v, t, vij, p, pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  v = 0.0; //reset observables
  t = 0.0;
  p = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       pij = 16.0/pow(dr, 12) - 8.0/pow(dr, 6);
       v += vij;
       p += pij/vol;
     }
      
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;			//Potential energy
    stima_kin = t/(double)npart; 			//Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; 	//Temperature
    stima_etot = (t+v)/(double)npart; 			//Total enery
    stima_pr = rho * stima_temp + p;			//Pressure

    ave_epot += stima_pot/N; 	//Average Potential energy
    ave_ekin += stima_kin/N;	//Average Kinetic energy
    ave_etot += stima_etot/N;	//Average total energy
    ave_temp += stima_temp/N;	//Average temperature
    ave_pr += stima_pr/N;	//Average pressure


    if (nconf%(int)N == 0){	//Ogni N passi, i.e. alla fine di ogni blocco, aggiorna il valore delle medie nei blocchi
  
	nblock++;	
	Epot.open("out/ave_epot_eq(0).dat",ios::app);
  	Ekin.open("out/ave_ekin(0).dat",ios::app);
  	Temp.open("out/ave_temp_eq(0).dat",ios::app);
  	Etot.open("out/ave_etot_eq(0).dat",ios::app);
	Press.open("out/ave_press(0).dat",ios::app);
	

	epot_BlockMean += ave_epot;
	ekin_BlockMean += ave_ekin;
	etot_BlockMean += ave_etot;
	temp_BlockMean += ave_temp;
	pr_BlockMean += ave_pr;	
	
	epot_BlockVar += pow(ave_epot, 2);
	ekin_BlockVar += pow(ave_ekin, 2);
	etot_BlockVar += pow(ave_etot, 2);
	temp_BlockVar += pow(ave_temp, 2);
	pr_BlockVar += pow(ave_pr, 2);		

  Epot << epot_BlockMean/nblock <<" "<< sqrt((epot_BlockVar/nblock - pow(epot_BlockMean/nblock, 2))/nblock) << endl;
	Ekin << ekin_BlockMean/nblock <<" "<< sqrt((ekin_BlockVar/nblock - pow(ekin_BlockMean/nblock, 2))/nblock) << endl;
  Etot << etot_BlockMean/nblock <<" "<< sqrt((etot_BlockVar/nblock - pow(etot_BlockMean/nblock, 2))/nblock) << endl;
	Temp << temp_BlockMean/nblock <<" "<< sqrt((temp_BlockVar/nblock - pow(temp_BlockMean/nblock, 2))/nblock) << endl;
	Press << pr_BlockMean/nblock <<" "<< sqrt((pr_BlockVar/nblock - pow(pr_BlockMean/nblock, 2))/nblock) << endl;
  
	ave_epot = 0;	
  ave_ekin = 0;
  ave_etot = 0;
  ave_temp = 0;
	ave_pr = 0;     
	Epot.close();
  Ekin.close();
	Temp.close();
	Etot.close();
	Press.close();
  }

    if (nconf%10 == 0){		//printo ogni 10 passi
  
	Epot.open("out/output_epot(0).dat",ios::app);
  	Ekin.open("out/output_ekin_eq(0).dat",ios::app);
  	Temp.open("out/output_temp(0).dat",ios::app);
  	Etot.open("out/output_etot(0).dat",ios::app);
	Press.open("out/output_press(0).dat",ios::app);
	
   	Epot << stima_pot  << endl;
   	Ekin << stima_kin  << endl;
   	Temp << stima_temp << endl;
	Etot << stima_etot << endl;
	Press << stima_pr  << endl;
    	 
	Epot.close();
    	Ekin.close();
	Temp.close();
	Etot.close();
	Press.close();

    }

    return;
}


/*************************************************************************************/


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void OldConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print old configuration to file oldconfig.final " << endl << endl;
  WriteConf.open("oldconfig.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}




void ConfXYZ(int conf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(conf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}



double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}



void Interazione2(){

char a;
  cout << "equilibration phase? Y/N" << endl;
  cin >> a ;
  while( a != 'Y' && a != 'N'){
  	cout << endl << "Invalid input, try again! [Y/N]" << endl;
  	cin >> a;
  }
  if (a == 'Y') {
  	ctrl = 1;
 	cout << "specify the number of steps (suggestion:1000)" << endl;
	cin >> restart_steps;
  	while(restart_steps < 100 || restart_steps > 10000) {
		cout << "Invalid input, try again." << endl;
		cin >> restart_steps;	
	}
  	cout << "Starting the equilibration phase... " << endl << endl;
  }
  else cout << "Starting the simulation..." << endl << endl;

return;
}

bool Interazione1(){

  char a;
  ifstream input;
  cout << "Do you want to start with the last saved configurations?[Y/N]..." << endl;
  cin >> a ;
  while( a != 'Y' && a != 'N'){
  	cout << endl << "Invalid input, try again! [Y/N]" << endl;
  	cin >> a;
  }
  if (a == 'Y') {
  	input.open("config.final");
	if(input){
		cout << "Loading configurations..." << endl;
		for (int i=0; i<npart; ++i){
    			input >> x[i] >> y[i] >> z[i];
			x[i] *= box;
			y[i] *= box;
			z[i] *= box;
  		}
	}
	
	else {
		
		cout << "ERROR: file config.final doesn't exist! The program will use fcc configuration" << endl;
		return 0;	
	}
	
	input.close();
	input.open("oldconfig.final");
	if(input){
		cout << "Loading old configurations..."<< endl;
		for (int i=0; i<npart; ++i){
    			input >> xold[i] >> yold[i] >> zold[i];
			xold[i] *= box;
			yold[i] *= box;
			zold[i] *= box;
  		}
	}
	
	else {
		cout << "ERROR: file oldconfig.final doesn't exist! The program will use fcc configuration" << endl;
		return 0;	
	}
  
  	cout << "Do you want to rescale the position in order for the velocities to match the desired temperature?" << endl;
   	cin >> a;
   	while( a != 'Y' && a != 'N'){
  		cout << endl << "Invalid input, try again! [Y/N]" << endl;
  		cin >> a;
   	}
   	
	if(a == 'Y') Rescale();
	return 1;
  }

else return 0; 

}

double get_fs (void){			//Velocity scale factor
  double v[3] = {0.0, 0.0, 0.0};

  double sumv[3] = {0.0, 0.0, 0.0};
  //Calcola fs 
  for (int i=0; i<npart; ++i){
  sumv[0] += vx[i];
  sumv[1] += vy[i];
  sumv[2] += vz[i];
  }

  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   
  double sumv2 = 0.0;
   
  for (int i=0; i<npart; ++i){
     v[0] = vx[i] - sumv[0];
     v[1] = vy[i] - sumv[1];
     v[2] = vz[i] - sumv[2];

     sumv2 += v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
   }
   sumv2 /= (double)npart;

   return sqrt(3 * temp / sumv2); 
}

void Rescale(void){

 double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - x[i])/delta;
    vy[i] = Pbc(ynew - y[i])/delta;
    vz[i] = Pbc(znew - z[i])/delta;

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }			

//infine qui aggiungo le posizioni vecchie tramite le velocità riscalate

  double fs = get_fs();

  for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }

  return;        
}


void Equilibrate(void){
	for(int i = 1; i <= restart_steps; i++) Move();	
	ConfFinal();					// Salva le attuali configurazioni
	Rescale();					// Riscala le posizioni vecchie compatibilmente con la temperatura desiderata
  return;
}


double get_temp(void){
  double t = 0;
  
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
 
  return (2.0 / 3.0) * t/(double)npart;
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
