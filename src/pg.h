#ifndef PG_Sampler
#define PG_Sampler



#include "Node.h"
#include "Queue.h"
#include "Likelihood.h"
#include "global.h"
#include "Lib.h"
#include "Rlob.h"
#include "Prior.h"

#include <cfloat>
//#include <Random>
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

    void CopyFrom(Particle* src);

    void SetFlag();

    void retrieve();

};

void RunSample(Node* thetree);

void InitParticles(Particle** particle_vec, double* log_weights, int len);

bool GrowParticle(Particle* p, Node** pgrow_node);

double UpdateWeight(Node* gnode);

void Resample(Particle** particle_vec, double* log_weight_vec, int size);

void ReleaseParticle(Particle* particle);

bool DrValidSplit(Node* gnode);

int SelectParticle(Particle** particle_vec, double* weight_vec, int size);

int PGLowerBound(double* vec, int len, double key);

bool CheckGrow(Particle** particle_vec, int size);

#endif
