#ifndef _PARTICLE_
#define _PARTICLE_

#include "./pvector.hpp"
#include "./params.hpp"
using ntype = double;
class particle
{
    // 03/12/24 ADDED forces and velocities
    pvector<ntype, 3> rold; // to store particle's position
protected:
    ntype vcut;

public:
    ntype sigma, epsilon, rc, m;
    pvector<ntype, 3> r, v, f; // particle's position and velocity

    // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
    // this method should be called when initializing particles
    void set_vcut(void)
    {
        ntype rijsq, srij2, srij6, srij12;
        ntype epsilon4 = epsilon * 4.0;
        rijsq = rc * rc;
        srij2 = sigma * sigma / rijsq;
        srij6 = srij2 * srij2 * srij2;
        srij12 = srij6 * srij6;
        vcut = epsilon4 * (srij12 - srij6);
    }

    void tra_move_plus(pvector<ntype, 3> delr)
    {
        r += delr;
    }

    void tra_move_minus(pvector<ntype, 3> delr)
    {
        r += delr;
    }
    void assign_position(pvector<ntype, 3> position) 
    {
        r=position;
    }

    void store()
    {
        rold = r;
    }

    void restore()
    {
        r = rold;
    }

    void set_sigma(ntype sig)
    {
        sigma = sig;
    }

    void set_epsilon(ntype eps)
    {
        epsilon = eps;
    }

    void set_rcut(ntype rcut)
    {
        rc = rcut;
    }

    void set_mass(ntype mass)
    {
        m = mass;
    }

    particle()
    {
        sigma = 1.0;
        epsilon = 1.0;
        rc = 2.5;
        vcut = 0.0;
        m = 1.0;
    }

    // methods for MD
    void expiLp(ntype dt)
    {
        // [TODO] implementation of action of operator exp(iLp*dt)
        v += f * (dt / m);
    }

    void expiLq(ntype dt)
    {
        // [TODO] implementation of action of operator exp(iLq*dt)
        r += v * dt;
    }
};

class particleLJ : public particle
{
public:
    ntype vij(particleLJ P, pvector<ntype, 3> L)
    {
        ntype ene;
        pvector<ntype, 3> Dr;

        Dr = r - P.r;
        // MINIMUM IMAGE CONVENTION
        Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
        ntype rsq, rn = Dr.norm();
        rsq = rn * rn;
        if (rsq < rc * rc) // interaction potential cut-off
            ene = 4.0 * epsilon * (pow(sigma / rn, 12.0) - pow(sigma / rn, 6));
        else
            ene = 0.0;
        return ene;
    }

    // vcut = vLJ(rc), i.e. the value of LJ potential at r=rc
    // this method should be called when initializing particles
    void set_vcut(void)
    {
        ntype rijsq, srij2, srij6, srij12;
        ntype epsilon4 = epsilon * 4.0;
        rijsq = rc * rc;
        srij2 = sigma * sigma / rijsq;
        srij6 = srij2 * srij2 * srij2;
        srij12 = srij6 * srij6;
        vcut = epsilon4 * (srij12 - srij6);
    }

    pvector<ntype, 3> fij(particleLJ P, pvector<ntype, 3> L, ntype &vij, ntype &vijs, ntype &wij)
    {
        // L is a vector with box dimensions
        // vij will be the interaction potential between i and j
        // vijs will be the shifted interaction potential
        // wij is the virial (i.e. rij*fij) for calculating the pressure
        ntype fij, srij2, srij6, srij12;
        ntype rijsq, epsilon24 = epsilon * 24.0;
        ntype epsilon4 = epsilon * 4.0;
        pvector<ntype, 3> Dr, fijv;
        Dr = r - P.r;

        // minimum image convention
        Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
        ntype rn = Dr.norm();
        if (rn >= rc)
        {
            wij = vij = vijs = 0;
            fijv = {0, 0, 0};
            return fijv;
        }
        rijsq = rn * rn;
        srij2 = sigma * sigma / rijsq;
        srij6 = srij2 * srij2 * srij2;
        srij12 = srij6 * srij6;
        vij = srij12 - srij6;

        // virial (for calculating pressure)
        wij = vij + srij12;
        wij *= epsilon24;

        // modulus of the force divided by modulus of rij
        // fij = 24*epsilon*(2*(sigma/r)^12 - (sigma/r)^6)
        fij = wij / rijsq;

        // force between two atoms
        fijv = fij * Dr;
        vij = epsilon4 * vij;
        vijs = vij - vcut;
        // return the force acting between "calling" particles and particle P
        return fijv;
    }
};

class particleSW : public particle
{
public:
    ntype sigma_SW = 1., delta_SW = 0.05 * sigma_SW, mu = 1.;

    void set_delta_SW(ntype delta)
    {
        delta_SW = delta;
    }

    void set_sigma_SW(ntype sigma)
    {
        sigma_SW = sigma;
    }

    void set_mu(ntype mu_input)
    {
        mu = mu_input;
        // pars.T = 1./mu; //Lo faccio nel run della MCsim
    }

    ntype get_mu()
    {
        return mu;
    }

    ntype get_sigma()
    {
        return sigma_SW;
    }
    
    ntype get_delta()
    {
        return delta_SW;
    }

    ntype vij(particleSW P, pvector<ntype, 3> L)
    {
        ntype ene = 0.;
        pvector<ntype, 3> Dr;

        Dr = r - P.r;
        // MINIMUM IMAGE CONVENTION
        Dr = Dr - L.mulcw(rint(Dr.divcw(L))); // Dr - L*rint(Dr/L)
        ntype rsq, rn = Dr.norm();
        rsq = rn * rn;
        if (rn > sigma_SW && rn <= sigma_SW + delta_SW)
            ene = -mu; // Square well attraction
        else if (rn > sigma_SW + delta_SW)
            ene = 0.0;
        else if (rn <= sigma_SW)
            ene = 1e10; // Hard core repulsion
        return ene;
    }

    ntype aggregation_volume()
    {
        const double pi = 3.141592653589793;
        // Hard-core volume
        ntype innerVolume = (4.0 / 3.0) * pi * std::pow(sigma_SW, 3);
        // Interaction volume
        ntype  outerVolume = (4.0 / 3.0) * pi * std::pow(delta_SW + sigma_SW, 3);
        // Aggregation volume
        ntype aggregationVolume = outerVolume - innerVolume;

        return aggregationVolume;
    }
};

#endif