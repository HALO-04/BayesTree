#ifndef PG_Sampler
#define PG_Sampler

#include "Node.h"
#include "Queue.h"
#include "global.h"

class Cachetemp{
public:
    //param
    double mu_prec;
    double mu_mean;
    double lambda_bart;
    double log_lambda_bart;

    //cache
    double half_log_2pi;
    double nn_prior_term;
    double sumy;
    double sumy2;
    double loglik;
}

void InitCachetemp(Cachetemp* cachetemp, double m_bart, double k_bart, double mlambda_bart);

void UpdateCachetemp(Cachetemp* cachetemp, double m_bart);

class Particle{
public:
    Node* thetree;

    //variable used for resample in pg sampler
    double log_weight;

    //queue for expansion
    Queue equeue;

    bool growable;

}

void Runsample(Particle** particle_vec, int len, Cachetemp* cachetemp);

void Initparticles(Particle** particle_vec, int len, const Cachetemp* cachetemp);

void Setparticlebytree(Particle* particle, Node* oldtree);

void Settreebyparticle(Node* dsttree, Particle* particle);

void Releaseparticle(Particle* particle);
#endif
