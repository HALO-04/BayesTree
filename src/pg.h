#ifndef PG_Sampler
#define PG_Sampler



#include "Node.h"
#include "Queue.h"
#include "Likelihood.h"
#include "global.h"
#include "Lib.h"

#include <float.h>
#include <Random>
#include <R.h>
#include <Rmath.h>

/*
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
*/

class Particle{
public:
    Node* thetree;

    //queue for expansion
    Queue equeue;

    bool growable;

}

void RunSample(Particle** particle_vec, int len, Cachetemp* cachetemp);

void InitParticles(Particle** particle_vec, int len, const Cachetemp* cachetemp);

void SetParticlebyTree(Particle* particle, Node* oldtree);

void SetTreebyParticle(Node* dsttree, Particle* particle);

void ReleaseParticle(Particle* particle);

bool DrValidSplit(Node* gnode);

int PGLowerBound(int *vec, int len);


#endif
