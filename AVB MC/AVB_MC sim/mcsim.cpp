#include "./sim.hpp"
int main(int argc, char **argv)
{
	mcsim<particleSW> mc;
    particleSW particles;
	mc.init_rng(0);			   // inizializzo il generatore di numeri casuali
	mc.prepare_initial_conf(); // creo la configurazione iniziale
	mc.run();				   // lancio la simulazione MC
	return 0;
}


