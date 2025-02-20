#ifndef _SIMCLASS_
#define _SIMCLASS_

#include <vector>
#include <string>
#include <fstream>
#include "./params.hpp"
#include "./particle.hpp"
#include "./randnumgen.hpp"
#include <iomanip> // for setprecision()

using ntype = double;

template <typename particle_type>
class sim
{
    using simp = simpars;

protected:
    simp pars; // parametri per composizione
    std::vector<particle_type> parts;

    // if opt=1 calculate energies only for i < j
    ntype calcenergyi(int i, int opt = 0)
    {
        int j;
        ntype enei = 0.0;
        for (j = 0; j < pars.Np; j++) // pars.Np è il numero totale di particelle
        {
            if (opt == 1 && i >= j)
                continue;
            if (i == j)
                continue;
            // la classe particelle deve implementare un metodo vij per il calcolo dell'energia d'interazione
            enei += parts[i].vij(parts[j], pars.L);
            // pars.L è un vettore con i lati del box
        }
        return enei;
    }

    ntype totenergy(void)
    {
        ntype ene = 0.0;
        for (auto i = 0; i < parts.size(); i++)
        {
            ene += calcenergyi(i, 1);
        }
        return ene;
    }

    void pbc(int i)
    {
        auto Dr = parts[i].r;
        Dr = pars.L.mulcw(rint(Dr.divcw(pars.L))); // L*rint(Dr/L)
        parts[i].r = parts[i].r - Dr;
    }

    void save_mgl_snapshot(long int t)
    {
        std::fstream f;
        std::string s;

        s = "mgl_frames/cnf-" + std::to_string(t) + ".mgl";
        f.open(s, std::ios::out | std::ios::trunc);
        for (int i = 0; i < pars.Np; i++)
        {
            f << parts[i].r(0) << " " << parts[i].r(1) << " " << parts[i].r(2) << " @ " << pars.sigma * 0.5 << "\n";
        }
        f.close();
    }

    void save_xyz_snapshot(int nconf)
    {
        ofstream WriteXYZ;
        WriteXYZ.open("xyz_frames/config_" + to_string(nconf) + ".xyz");
        WriteXYZ << pars.Np << endl;
        WriteXYZ << "This is only a comment!" << endl;
        for (int i = 0; i < pars.Np; ++i)
        {
            WriteXYZ << pars.sigma * 0.5 << " " << parts[i].r(0) << "   " << parts[i].r(1) << "   " << parts[i].r(2) << endl;
        }
        WriteXYZ.close();
    }

public:
    void prepare_initial_conf(void)
    {
        // SC
        int ix, iy, iz;
        int cc = 0;
        pars.Np = pars.nx * pars.ny * pars.nz;
        parts.resize(pars.Np);
        ntype vcell = pow(pars.sigma, 3.0);
        ntype rhomax = 1.0 / vcell;
        ntype sf;
        sf = cbrt(rhomax / pars.rho);
        pars.L = {ntype(pars.nx), ntype(pars.ny), ntype(pars.nz)};
        ntype clen = sf * pars.sigma;
        pars.L *= clen;
        for (ix = 0; ix < pars.nx; ix++)
            for (iy = 0; iy < pars.ny; iy++)
                for (iz = 0; iz < pars.nz; iz++)
                {
                    parts[cc].r = {ix * clen, iy * clen, iz * clen};
                    parts[cc].r -= pars.L * 0.5;
                    parts[cc].set_sigma(pars.sigma);
                    parts[cc].set_epsilon(pars.epsilon);
                    parts[cc].set_rcut(pars.rc);
                    cc++;
                }
        // ...or BCC or FCC lattice
    }

    void set_sim_type(int type)
    {
        pars.simtype = type;
    }

    void init_rng(int n)
    {
        if (n < 0)
            rng.rseed();
        else
            rng.seed(n);
    }

    void run(void) {
        // intentionally void
    };
};

template <typename particle_type>
class mcsim : public sim<particle_type>
{
    // for calc_acceptance_and_adjust: total trial moves and accepted ones
    // for calculating acceptance rates.
    using bc = sim<particle_type>;
    using bc::calcenergyi, bc::parts, bc::pars,
        bc::totenergy, bc::save_mgl_snapshot, bc::pbc, bc::save_xyz_snapshot;

    particle_type particles;
    string mu_string;
    string sim_type_string;

    // counters used to calculate acceptance rates
    long int tot_tra,
        tot_vol, tra_rej, vol_rej, bias_rej, tot_bias;

    void alpha(int i)
    {
        // muovere la particella usando come max
        // displacement pars.deltra
        // applico le pbc (tramite il metodo pbc)
        pvector<ntype, 3> delr;
        delr = {pars.deltra * 2.0 * (rng.ranf() - 0.5), pars.deltra * 2.0 * (rng.ranf() - 0.5),
                pars.deltra * 2.0 * (rng.ranf() - 0.5)};
        parts[i].tra_move_plus(delr);
        pbc(i);
    }

    void acc(int i, ntype eno)
    {
        // assume kB=1 (Boltzmann constants) so that beta=1/T
        // calcolo la nuova energia
        // accetto la mossa con probabilità min{1taU)}
        // se deltaU < 0 accetto la mossa
        // altrimenti genero xi un numero a caso tra 0 e 1
        // e se xi < exp(-beta*deltaU) accetto altrimenti rifiuto
        // accetta la mossa con criterio Metropolis MC
        ntype enn = calcenergyi(i);
        ntype delu = enn - eno;
        ntype xi = rng.ranf();
        if (delu > 0.0 && xi >= exp(-delu / pars.T))
        {
            // reject move
            tra_rej++;
            parts[i].restore();
        }
    }

    void move_NTV(int i)
    {
        // 1) calcolare l'energia della particella i-esima
        ntype eno;
        eno = calcenergyi(i);
        // 2) memorizzo la posizione della particella i
        parts[i].store();
        // 3) trial move
        alpha(i);
        // 4) acceptance per la particella i
        acc(i, eno);
    }
    void calc_acceptance_and_adjust(void)
    {
        // CALC ACCEPTANCE RATES AND
        // ADJUST pars.deltra
        ntype r_tra, r_vol, r_bias;
        if (tot_tra > 0)
        {
            // rate di accettazione delle mosse
            r_tra = ((double)(tot_tra - tra_rej)) / tot_tra;
            std::cout << "rate tra: " << r_tra << " deltra=" << pars.deltra << "\n";
            if (r_tra > 0.5)
            {
                pars.deltra *= 1.1;
            }
            else
            {
                pars.deltra /= 1.1;
            }
            tot_tra = tra_rej = 0;
        }

        // adjust maximum volume "displacement" in alpha_box
        // so that acceptance rate of box move is around 0.5
        if (tot_vol > 0) // <------------------------------------------------
        {
            r_vol = ((ntype)(tot_vol - vol_rej)) / tot_vol;
            std::cout << "rate vol: " << r_vol << " vmax=" << pars.vmax << "\n";
            if (r_vol > 0.5)
            {
                pars.vmax *= 1.1;
            }
            else
            {
                pars.vmax /= 1.1;
            }
            tot_vol = vol_rej = 0.0;
        }

        // adjust maximum biased "displacement" in the aggregation ring
        // so that acceptance rate of biased move is around 0.5
        if (tot_bias > 0) // <------------------------------------------------
        {
            r_bias = ((ntype)(tot_bias - bias_rej)) / tot_bias;
            std::cout << "rate biased move: " << r_bias << " rmax=" << pars.bias_max << "\n";
            if (r_bias > 0.5)
            {
                pars.bias_max *= 1.1;
            }
            else
            {
                pars.bias_max /= 1.1;
            }
            tot_bias = bias_rej = 0.0;
        }
    }

    void init_measures(void)
    {
        // open files in writing mode to truncate them to 0
        std::fstream f;
        f.open("mu=" + mu_string + "_energy_" + sim_type_string + ".dat", std::ios::out | std::ios::trunc);
        f.close();
        if (pars.simtype == 1) // simtype == 1 means an NPT simluations
        {
            // clear file used to store densities in NPT simulation
            f.open("density.dat", std::ios::out | std::ios::trunc);
            f.close();
        }
    }

    // ntype ratio = pars.sigma/pars.rc;
    // ntype vtail = (8.*3.14159265358979323846/3.)*pars.epsilon*pars.rho*pow(pars.sigma,3)*((1./3.)*pow(ratio,12)-pow(ratio,3)) ;

    void save_measures(long int t)
    {
        std::fstream f;
        // f.open("energy.dat", std::ios::out | std::ios::app);
        f.open("mu=" + mu_string + "_energy_" + sim_type_string + ".dat", std::ios::out | std::ios::app);
        // save potential energy per particle
        f << t << " " << totenergy() / pars.Np + pars.vtail << "\n";
        f.close();
        if (pars.simtype == 1) // 1 means NPT simulation, save density in this case
        {
            f.open("density.dat", std::ios::out | std::ios::app);
            f << t << " " << pars.Np / (pars.L(0) * pars.L(1) * pars.L(2)) << "\n";
            f.close();
        }
    }

    void restore_all_pars()
    {
        // restore all particle positions
        for (int i = 0; i < pars.Np; i++)
        {
            parts[i].restore();
        }
        pars.L = pars.old_L;
    }

    void store_all_pars()
    {
        // store all particle position
        for (int i = 0; i < pars.Np; i++)
        {
            parts[i].store();
        }
        pars.old_L = pars.L;
    }

    void alpha_box(ntype &DG, ntype &fact)
    {

        // trial box move
        //
        // 1) choose randomly a new volume
        //
        // 2) scale all particle positions and box size by factor "fact"
        //
        // 3) calculate new total interaction energy
        //
        // 4) calculate \DeltaG

        ntype old_U = totenergy();
        ntype old_V = pars.L(1) * pars.L(2) * pars.L(3);
        ntype new_V = old_V + pars.vmax * (rng.ranf() * 2.0 - 1.0);
        fact = pow(new_V / old_V, 1. / 3.); // fact = cbrt(Vn/Vo);
        pars.L *= fact;
        for (int i = 0; i < pars.Np; i++)
        {
            parts[i].r *= fact;
        }
        ntype new_U = totenergy();
        DG = pars.P * (new_V - old_V) + (new_U - old_U) - pars.Np * pars.T * log(new_V / old_V);

        /* Prof
                ntype delta_V = pars.vmax*(rng.ranf()*2.0-1.0);
                ntype old_U = totenergy();
                ntype old_V = pars.L(0)*pars.L(1)*pars.L(2);
                ntype new_V = old_V + delta_V;
                cout << "new Volume " << new_V << endl;
                fact = cbrt(new_V/old_V); // (Vn/Vo)^(1/3)
                for (int i=0; i < pars.Np; i++)
                {
                    parts[i].r *= fact;
                }
                pars.L *= fact;
                ntype new_U = totenergy();
                DG = pars.P * delta_V + (new_U-old_U) - pars.T*pars.Np*log(new_V/old_V);
        */
    }

    void acc_box(ntype DG, ntype fact)
    {
        // accept or reject box move
        // 1) calculate e^{-\Delta G} (see pdf of lectures)
        // 2) generate a random numner \xi and check wheter to accept the box trial move
        // 3) if move is rejected:
        //    i)   restore all particle positions thourgh method restore_all_pars()
        //    ii)  restore box size
        //    iii) update counter vol_rej of rejected box moves for calculating acceptance ratio
        ntype xi = rng.ranf();
        if (xi >= exp(-DG / pars.T))
        {
            // reject move
            restore_all_pars();
            vol_rej++;
            pars.L /= fact; // per versione Prof
        }
    }
    // Volume move for NPT ensemble
    void move_box(void)
    {
        ntype DG, fact;   // DG is \Delta G as discussed during lectures (see pdf of lectures)
        store_all_pars(); // store all particle positions and old L before attempting a box move
        alpha_box(DG, fact);
        acc_box(DG, fact);
    }

    void move_outin(int i)
    {
        int j;
        ntype rn;
        pvector<ntype, 3> Dr;
        int N_in = 0; // SONO DA CONTARE
        for (int l = 0; l < pars.Np; l++)
        {
            Dr = parts[i].r - parts[l].r;
            Dr = Dr - pars.L.mulcw(rint(Dr.divcw(pars.L))); // Dr - L*rint(Dr/L); Minimum image convention
            rn = Dr.norm();
            if (i != l && rn > particles.get_sigma() && rn < particles.get_sigma() + particles.get_delta())
            {
                N_in++;
            }
        }
        int N_out = pars.Np - N_in - 1;

        do
        {
            j = rng.ranf() * pars.Np;
            Dr = parts[i].r - parts[j].r;
            Dr = Dr - pars.L.mulcw(rint(Dr.divcw(pars.L))); // Dr - L*rint(Dr/L)
            rn = Dr.norm();
        } while (i == j || (rn > particles.get_sigma() && rn < particles.get_sigma() + particles.get_delta()));

        ntype eno, enn;
        eno = calcenergyi(j);
        // 2) memorizzo la posizione della particella j
        parts[j].store();
        // 3) trial move
        // MOVE J INTO BONDING VOLUME OF I:


//NUOVO https://math.stackexchange.com/questions/1111894/random-points-in-spherical-shell#:~:text=You%20could%20generate%20a%20random,origin%2C%20r%E2%80%B2%20the%20distance%20after
        pvector<ntype, 3> travel, direction;
        direction.random_orient();
        ntype radius, extracted = rng.ranf(pow(particles.get_sigma(), 3), pow(particles.get_sigma() + particles.get_delta(), 3));
        radius = cbrt(extracted);
        travel = radius*direction;
        parts[j].assign_position(travel);
        pbc(j);

        enn = calcenergyi(j);
        ntype V_in = particles.aggregation_volume();
        ntype V_tot = pars.L(1) * pars.L(2) * pars.L(3);
        ntype V_out = V_tot - V_in;

        ntype prefactor = (N_out * V_in) / (V_out * (N_in + 1));
        ntype xi = rng.ranf();
        ntype delu = enn - eno;
        if (delu > 0.0 && xi >= prefactor * exp(-delu / pars.T))
        {
            // reject move
            bias_rej++;
            parts[j].restore();
        }
    }

    void move_inout(int i)
    {
        int j;
        ntype rn;
        pvector<ntype, 3> Dr;
        int N_in = 0; // SONO DA CONTARE
        for (int l = 0; l < pars.Np; l++)
        {
            Dr = parts[i].r - parts[l].r;
            Dr = Dr - pars.L.mulcw(rint(Dr.divcw(pars.L))); // Dr - L*rint(Dr/L); Minimum image convention
            rn = Dr.norm();
            if (i != l && rn > particles.get_sigma() && rn < particles.get_sigma() + particles.get_delta())
            {
                N_in++;
            }
        }
        int N_out = pars.Np - N_in - 1;

        do
        {
            j = rng.ranf() * pars.Np;
            Dr = parts[i].r - parts[j].r;
            Dr = Dr - pars.L.mulcw(rint(Dr.divcw(pars.L))); // Dr - L*rint(Dr/L)
            rn = Dr.norm();
        } while (i == j || (rn < particles.get_sigma() && rn > particles.get_sigma() + particles.get_delta()));

        ntype eno, enn;
        eno = calcenergyi(j);
        // 2) memorizzo la posizione della particella j
        parts[j].store();
        // 3) trial move
        // MOVE J OUT THE BONDING VOLUME OF I:
        // we randomly insert j in the box with a random orientation till i and j are no longer bonded
        do
        {
            pvector<ntype, 3> delr;
            delr = {pars.deltra * 2.0 * (rng.ranf() - 0.5), pars.deltra * 2.0 * (rng.ranf() - 0.5),
                    pars.deltra * 2.0 * (rng.ranf() - 0.5)};
            parts[j].tra_move_plus(delr);
            pbc(j);

            Dr = parts[i].r - parts[j].r;
            Dr = Dr - pars.L.mulcw(rint(Dr.divcw(pars.L))); // Dr - L*rint(Dr/L)
            rn = Dr.norm();
        } while (rn > particles.get_sigma() && rn < particles.get_sigma() + particles.get_delta());

        enn = calcenergyi(j);
        ntype V_in = particles.aggregation_volume();
        ntype V_tot = pars.L(1) * pars.L(2) * pars.L(3);
        ntype V_out = V_tot - V_in;

        ntype prefactor = (N_in * V_out) / (V_in * (N_out + 1));
        ntype xi = rng.ranf();
        ntype delu = enn - eno;
        if (delu > 0.0 && xi >= prefactor * exp(-delu / pars.T))
        {
            // reject move
            bias_rej++;
            parts[j].restore();
        }
    }

public:
    void run(void)
    {
        int counter = 1;
        vector<ntype> mu = {1., 5., 10., 15., 20.};
        sim_type_string = to_string((int)pars.simtype);
        for (int w = 0; w < mu.size(); w++)
        {
            particles.set_mu(mu[w]);
            cout << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
            cout << "Value of Mu number " << counter << " with value: mu = " << particles.get_mu() << endl;
            counter++;
            mu_string = to_string((int)particles.get_mu());
            pars.T = 1. / particles.get_mu();

            // ciclo sui passi MC
            int i, t, ip;
            tot_tra = tot_vol = tra_rej = vol_rej = 0;
            tot_bias = bias_rej = 0;
            init_measures();
            for (t = 0; t < pars.totsteps; t++)
            {
                // ATTEMPT TO MOVE ALL PARTICLES:
                // ogni passo MC sono Np tentativi di mossa
                for (i = 0; i < pars.Np; i++)
                {
                    // choose between biased move and particle move
                    if (pars.simtype == 2 && rng.ranf() < 0.1) //< 0.1 per avere 10%
                    {
                        if (rng.ranf() < 0.5)
                        {
                            ip = rng.ranf() * pars.Np;
                            move_inout(ip);
                            tot_bias++;
                        }
                        else
                        {
                            ip = rng.ranf() * pars.Np;
                            move_outin(ip);
                            tot_bias++;
                        }
                    }
                    else
                    {
                        ip = rng.ranf() * pars.Np;
                        move_NTV(ip);
                        tot_tra++;
                    }
                }

                if (t > 0 && pars.savemeasure > 0 && t % pars.savemeasure == 0)
                {
                    save_measures(t);
                }

                if (t > 0 && t % pars.outstps == 0)
                {
                    std::cout << "Step #" << t << "\n";
                    // per confrontarsi con il Johnson si deve calcolare l'energia interna ovvero sommare il contributo
                    // di energia cinetica
                    std::cout << "total energy per particle is " << totenergy() / pars.Np << "\n";
                    std::cout << "box size: " << pars.L(0) << " " << pars.L(1) << " " << pars.L(2) << "\n";
                }

                if (t > 0 && pars.save_mgl_snapshot > 0 &&
                    t % pars.save_mgl_snapshot == 0)
                {
                    save_mgl_snapshot(t);
                    save_xyz_snapshot(t);
                }

                if (t > 0 && pars.adjstps > 0 && t % pars.adjstps == 0 &&
                    t < pars.maxadjstps)
                {
                    calc_acceptance_and_adjust();
                }
            }
        }
    }
};

#endif
