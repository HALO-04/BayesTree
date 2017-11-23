#include "pg.h"



void InitCachetemp(Cachetemp* cachetemp, double m_bart, double k_bart, double mlambda_bart){
	this->mu_mean = 0;
	this->mu_prec = m_bart * (2 * k_bart) * (2 * k_bart);
	this->half_log_2pi = 0.5 * std::log(2 * PI);
	this->nn_prior_term = 0.5 * std::log(mu_prec);
	this->lambda_bart = mlambda_bart;
	this->log_lambda_bart = std::log(mlambda_bart);

	UpdateCachetemp(cachetemp, m_bart);
}

void UpdateCachetemp(Cachetemp* cachetemp, double m_bart){
	this->sumy = 0;
	this->sumy2 = 0;
	for(int i = 1; i <= NumObs; i++){
		this->sumy += YDat1[i];
		this->sumy2 += (YDat1[i] / m_bart) * (YDat1[i] / m_bart);
	}
	this->sumy /= NumObs;

	double mu_prec_post = this->lambda_bart * NumObs + this->mu_prec;
	double mu_mean_post = (this->lambda_bart * this->sum_y) / mu_prec_post;
	//double mu_mean_post = (this->mu_prec * this->mu_mean + this->lambda_bart * this->sum_y) / mu_prec_post;

	this->loglik = this->nn_prior_term - NumObs * this->half_log_2pi
					+ 0.5 * (NumObs * this->log_lambda_bart - std::log(mu_prec_post))
					+ mu_prec_post * mu_mean_post * mu_mean_post - this->lambda_bart * this->sumy2;
}


void Runsample(Node* thetree, Cachetemp* cachetemp){

    Particle**  particle_vec = new Particle*[NumParticle + 1];
    Initparticles(particle_vec, NumParticle, cachetemp);
    Particle* first_particle = *(particle_vec + 1);
    Setparticlebytree(first_particle, thetree);

    while(true){
        double sum_weight = first_particle->weight;
        for(int i = 2; i <= len; i++){
            Particle* cur_particle = *(particle_vec + i);
            bool done = Growparticle(cur_particle);
            if(done){
                Updateweight(cur_particle);
            }
            sum_weight += cur_particle->weight;
        }
        for(int i = 1; i < NumParticle; i++) *(particle_vec + i)->weight /= sum_weight;
        resample(particle_vec, NumParticle);
        if(!Checkgrow(particle_vec, NumParticle)) break;
    }

    int select_idx = Selectparticle(particle_vec, NumParticle);
    Settreebyparticle(thetree, *(particle_vec + idx));

    for(int i = 1; i <= NumParticle; i++) delete *(particle_vec + i);
    delete[] particle_vec;
}




void Initparticles(Particle** particle_vec, int len, const Cachetemp* cachetemp){
    for(int i = 1; i <= len; i++){
        Particle* cur_particle = *(particle_vec + i);
        cur_particle->thetree = new Node;
        cur_particle->equeue.append(cur_particle->thetree);
        cur_particle->growable = true;

        cur_particle->weight = cachetemp->loglik;
    }
}

void Growparticle(Particle* p){

}


void Setparticlebytree(Particle* particle, Node* oldtree);

void Settreebyparticle(Node* dsttree, Particle* particle);

void Releaseparticle(Particle* particle){
    if(!(particle->thetree))
        delete particle->thetree;
}
