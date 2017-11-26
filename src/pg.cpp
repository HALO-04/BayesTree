#include "pg.h"


/*
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
*/

void Runsample(Node* thetree, Cachetemp* cachetemp){

    Particle**  particle_vec = new Particle*[NumParticle + 1];
    double* weight_vec = new double[NumParticle + 1];
    Initparticles(particle_vec, weight_vec, NumParticle, cachetemp);

    Particle* first_particle = *(particle_vec + 1);
    thetree->CopyTree(first_particle->thetree);

    while(true){
        double sum_weight = weight_vec[1];
        for(int i = 2; i <= len; i++){
            Particle* cur_particle = *(particle_vec + i);
            bool done = Growparticle(cur_particle);
            if(done){
                Updateweight(cur_particle, weight_vec, i);
            }
            sum_weight += weight_vec[i];
        }
        for(int i = 1; i < NumParticle; i++) weight_vec[i] /= sum_weight;
        resample(particle_vec, weight_vec, NumParticle);
        if(!Checkgrow(particle_vec, NumParticle)) break;
    }

    int select_idx = Selectparticle(particle_vec, weight_vec, NumParticle);
    Particle* selected = *(particle_vec + select_idx);
    thetree->deall();
    selected->CopyTree(thetree);

    for(int i = 1; i <= NumParticle; i++) delete *(particle_vec + i);
    delete[] particle_vec;
    delete[] weight_vec;
}




void Initparticles(Particle** particle_vec, double* weight_vec, int len, const Cachetemp* cachetemp){
    Particle* cur_particle;
    for(int i = 1; i <= len; i++){
        cur_particle = *(particle_vec + i);
        cur_particle->thetree = new Node;
        cur_particle->thetree->SetData();

        cur_particle->equeue.append(cur_particle->thetree);
        cur_particle->growable = true;

        weight_vec[i] = cachetemp->loglik;
    }

    double loglik = LogLNode(cur_particle->thetree);
    for(int i = 1; i <= len; i++) weight_vec[i] = loglik;

}

bool Growparticle(Particle* p){
    Queue* q = &(p->equeue);

    if(q->empty()) return false;

    Node* grow_node = (Node*)q->pop();
    double psplit = PriParams.base / pow(1.0+Depth(grow_node), PriParams.power);

    bool status;
    int split_var;
    int split_idx;

    if(!Bern(psplit))
        return false;
    status = DrValidSplit(grow_node, &split_var, &split_idx);
    if(!status)
        return false;



}


bool DrValidSplit(Node* gnode, int* var, int* split_idx){
    int Ngood = SumGoodVar(gnode);

    int* n_dim = new int[Ngood + 1];

    int i;
    int k = 0;
    for(i = 1; i <= NumObs; i++) {
        if(gnode->VarAvail[i]) {
            k++;
            n_dim[k] = i;
        }
    }
    std::random_shuffle(n_dim + 1, n_dim + Ngood + 1);

    int tvar, length;
    int first_bound, second_bound;
    double x_min, x_max, tmp_value;

    Cell* cur_cell;
    int* pIvec;

    length = grow_node->DataList.length;

    for(i = 1; i <= Ngood; i++){
        tvar = n_dim[i];
        cur_cell = grow_node->DataList.first;
        // risk of cur_cell = NULL

        if(VarType[tvar] == ORD){
            x_min = DBL_MAX;
            x_max = -1.0 * DBL_MAX;
            while(true){
                pIvec = (int*)cur_cell->contents;
                tmp_value = XDat[*pIvec][tvar];
                if(tmp_value > x_max)
                    x_max = tmp_value;
                if(tmp_value < x_min)
                    x_min = tmp_value;
                if(cur_cell->End)
                    break;
                cur_cell = cur_cell->after;
            }
            first_bound = PGLowerBound(RuleMat[tvar] + 1, RuleNum[tvar], x_min);
            second_bound = PGLowerBound(RuleMat[tvar] + 1, RuleNum[tvar], x_max);
            if(first_bound < second_bound){
                (*var) = tvar;

                //need a uniform sample

                (*split_idx) = first_bound + 1;
                delete[] n_dim;
                return true;
            }
        }else{
            // CAT cases
        }
    }

    delete[] n_dim;

}

int PGLowerBound(int *vec, int len, double value){

}


void Releaseparticle(Particle* particle){
    if(!(particle->thetree))
        delete particle->thetree;
}

bool Checkgrow(particle_vec, NumObs){
    Particle* cur;
    for(int i = 1; i <= NumObs; i++){
        cur = *(particle_vec + i);
        if(cur->growable)
            return true;
    }
    return false;
}



