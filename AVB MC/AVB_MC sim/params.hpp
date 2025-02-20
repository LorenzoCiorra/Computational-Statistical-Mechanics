#ifndef _PARAMS_
#define _PARAMS_
#include "./pvector.hpp"
#include <fstream>
using namespace std;

class simpars
{
	using ntype = double;

public:
	int nx, ny, nz; /* nx*ny*nz particelle */
	double T, P;	// temperature and pressure
	int Np;			// numero di particelle
	long int maxadjstps, eqstps, adjstps, save_mgl_snapshot;
	long int savemeasure, outstps, totsteps; // Nsteps = simulations steps, outstps steps print something on current simulation status
	double rho, rc;							 // density
	int simtype;							 // simulation type (see below)
	int seed;								 // -1 means random
	pvector<double, 3> L, old_L;			 // box
	double sigma, epsilon, mass;			 // Lennard-Jones parameters
	double dt, deltra, vmax, bias_max;		 // parameter of MC moves
	ntype vtail;							 // tail correction to the potential energy
	simpars()
	{
		// qui leggo il file input.in e metto il risultato in tutti i parametri qui sotto
		ifstream ReadInput;
		string s;
		cout << "2" << endl;
		// Read input informations
		ReadInput.open("input.in");

		// cout << "Classic Lennard-Jones fluid        " << endl;

		ReadInput >> simtype >> s; // 0 NTV, 1 NPT
		if (simtype == 0)
		{
			cout << "MC(NVT) simulation       " << endl
				 << endl;
		}
		if (simtype == 1)
		{
			cout << "MC(NPT) simulation       " << endl
				 << endl;
		}
		if (simtype == 2)
		{
			cout << "MC(NVT) simulation w/ Aggregation Volume Bias       " << endl
				 << endl;
		}
		cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl
			 << endl;

		ReadInput >> nx >> s; // number of particles along each direction
		ReadInput >> ny >> s;
		ReadInput >> nz >> s;
		cout << "Number of particles along x " << nx << " along y " << ny << " along z " << nz << endl;
		ReadInput >> sigma >> s;
		ReadInput >> epsilon >> s;
		ReadInput >> rho >> s;
		//		L = {pow(nx/rho, 1./3.), pow(ny/rho, 1./3.), pow(nz/rho, 1./3.)};
		cout << "Density of particles = " << rho << endl;
		ReadInput >> rc >> s;
		cout << "Cutoff of the LJ interatomic potential = " << rc << endl
			 << endl;
		ReadInput >> seed >> s;
		ReadInput >> mass >> s;
		ReadInput >> adjstps >> s;
		ReadInput >> maxadjstps >> s;
		ReadInput >> eqstps >> s;
		ReadInput >> totsteps >> s;
		ReadInput >> save_mgl_snapshot >> s;
		ReadInput >> savemeasure >> s;
		ReadInput >> outstps;
		ReadInput >> T;
		//cout << "Temperature = " << T << endl;
		ReadInput >> P;
		ReadInput >> deltra;
		ReadInput >> vmax;
		ReadInput >> dt;
		ntype ratio = sigma / rc;
		vtail = (8. * 3.14159265358979323846 / 3.) * epsilon * rho * pow(sigma, 3) * ((1. / 3.) * pow(ratio, 12) - pow(ratio, 3));
		// cout << vtail << endl;
		ReadInput.close();
		bias_max = 1;

		/*
				simtype = 0; // 0 NTV, 1 NPT
				nx = 8;		 // number of particles along each direction
				ny = 8;
				nz = 8;
				sigma = 1.0;
				epsilon = 1.0;
				rho = 0.5; // rho=0.5 T=2.0 see johnson for potential energy
				rc = 5.0;
				seed = 0;
				mass = 1.0;
				adjstps = -200;
				maxadjstps = 2000;
				eqstps = 500;
				totsteps = 20000;
				save_mgl_snapshot = 1000;
				savemeasure = 20;
				outstps = 200;
				T = 2.0;
				P = 3.838; // se P*=\beta*P*v0, 1 < P* < 10 dove v0 è il volume di una particella
				deltra = 0.2;
				vmax = 10.0;
				dt = 0.002;
		*/
	}
};
#endif
