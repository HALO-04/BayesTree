#ifndef PG_Sampler
#define PG_Sampler

#include "Node.h"
#include "Queue.h"
#include "global.h"

extern NumParticles;

class Cachetemp{
public:
    //param
    double mu_prec;
    double mu_mean;

    //cache
    double half_log_2pi;
    double nn_prior_term;
    double sumy;
    double sumy2;
    double loglik;
}

void InitCachetemp(Cachetemp* cachetemp);


class Particle{
public:
    Node* thetree;

    //variable used for resample in pg sampler
    double log_weight;

    //queue for expansion
    Queue equeue;

}

void Initparticles(Particle** particle_vec, int len);

void Setparticlebytree(Particle* particle, Node* oldtree);

void Settreebyparticle(Node* dsttree, Particle* particle);

void Releaseparticle(Particle* particle);
#endif
