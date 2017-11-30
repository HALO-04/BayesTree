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

void Particle::CopyFrom(Particle* src){
    this->growable = src->growable;
    src->thetree->CopyTree(this->thetree);
    this->equeue.CopyFrom(&(src->equeue));
}

void RunSample(Node* thetree){

    Particle**  particle_vec = new Particle*[NumParticle + 1];
    double* log_weight_vec = new double[NumParticle + 1];
    InitParticles(particle_vec, log_weight_vec, NumParticle);

    Particle* first_particle = *(particle_vec + 1);
    first_particle->thetree->deall();
    thetree->CopyTree(first_particle->thetree);

    //only for compile test
    first_particle->growable = false;

    Node* gnode;

    int itr = 0;

    while(true){
        itr++;
        for(int i = 2; i <= NumParticle; i++){
            Particle* cur_particle = *(particle_vec + i);

            bool done = GrowParticle(cur_particle, &gnode);
            if(done){
                log_weight_vec[i] += UpdateWeight(gnode);
            }
        }
        Resample(particle_vec, log_weight_vec, NumParticle);
        if(!CheckGrow(particle_vec, NumParticle)) break;
    }

    int select_idx = SelectParticle(particle_vec, log_weight_vec, NumParticle);
    Particle* selected = *(particle_vec + select_idx);
    thetree->deall();
    selected->thetree->CopyTree(thetree);

    for(int i = 1; i <= NumParticle; i++) {
        ReleaseParticle(particle_vec[i]);
        particle_vec[i] = NULL;
    }
    delete[] particle_vec;
    delete[] log_weight_vec;
}


void InitParticles(Particle** particle_vec, double* weight_vec, int len){
    Particle* cur_particle;
    for(int i = 1; i <= len; i++){
        cur_particle = new Particle;
        cur_particle->thetree = new Node;
        cur_particle->thetree->SetData();
        cur_particle->equeue.append(cur_particle->thetree);
        cur_particle->growable = true;
        particle_vec[i] = cur_particle;
    }

    double loglik = LogLNode(cur_particle->thetree);
    for(int i = 1; i <= len; i++) weight_vec[i] = loglik;

}

bool GrowParticle(Particle* p, Node** pgrow_node){
    Queue* q = &(p->equeue);

    if(q->empty()){
        p->growable = false;
        return false;
    }

    Node* grow_node = (Node*)q->pop();
    double psplit = PriParams.base / pow(1.0+Depth(grow_node), PriParams.power);

    bool status;

    if(!Bern(psplit))
        return false;
    status = DrValidSplit(grow_node);
    if(!status)
        return false;

    int LeftEx, RightEx;
    SpawnChildren(grow_node, LeftEx, RightEx);
    (*pgrow_node) = grow_node;
    return true;
}


bool DrValidSplit(Node* gnode){
    int Ngood = SumGoodVar(gnode);
    int* n_dim = new int[Ngood + 1];
    int i;
    int k = 0;
    for(i = 1; i <= NumX; i++) {
        if(gnode->VarAvail[i]) {
            k++;
            n_dim[k] = i;
        }
    }
    Lib::shuffle(n_dim + 1, Ngood);

    int tvar, length;
    int first_bound, second_bound;
    double x_min, x_max, tmp_value;

    Cell* cur_cell;
    int* pIvec;

    length = gnode->DataList.length;

    for(i = 1; i <= Ngood; i++){
        tvar = n_dim[i];
        cur_cell = gnode->DataList.first;
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
                (gnode->rule).Var = tvar;
                double u = unif_rand();
                //need a uniform sample
                int offset = (int)floor(u * (second_bound - first_bound));
                (gnode->rule).OrdRule = first_bound + offset;
                delete[] n_dim;
                return true;
            }
        }else{
            // CAT cases
            Rprintf("\n\nPG Sampler only supports ORD x variable\n\n");
            exit(1);
            /*
            int* n_split = new int[ReluNum[tvar] + 1];
            while(true){
                pIvec = (int*)cur_cell->contents;
                tmp_value = XDat[*pIvec][tvar];

                if(cur_cell->End)
                    break;
                cur_cell = cur_cell->after;
            }*/
        }
    }
    delete[] n_dim;
    return false;
}

double UpdateWeight(Node* gnode){
    double old_log = LogLNode(gnode);
    double new_log = LogLNode(gnode->LeftC) + LogLNode(gnode->RightC);
    return (new_log - old_log);
}


int PGLowerBound(double *array, int size, double key){
    int first = 0, middle;
    int half, len;
    len = size;

    while(len > 0) {
        half = len >> 1;
        middle = first + half;
        if(array[middle] < key) {
            first = middle + 1;
            len = len-half-1;       //sreach in the rigth part
        }
        else
            len = half;            //sreach in the left part
    }
    return first + 1;
}

void Resample(Particle** particle_vec, double* log_weight_vec, int size){
    double* weight_norm = new double[size + 1];
    int* sample_index = new int[size + 1];
    int i;
    double tmax, sum_log, log_numP, log_pd;

    sum_log = Lib::softmax(log_weight_vec, weight_norm, size, tmax);
    log_pd = log(sum_log) + tmax;
    log_numP = log(size);

    for(i = 1; i <= size; i++)
        log_weight_vec[i] = log_pd - log_numP;

    Lib::SampleMultinomial(weight_norm, size, sample_index, size);
    Particle** new_particle_vec = new Particle*[size + 1];

    for(i = 2; i <= size; i++){
        Particle* new_particle = new Particle;
        new_particle->thetree = new Node;
        new_particle->CopyFrom(particle_vec[sample_index[i]]);
        new_particle_vec[i] = new_particle;
    }
    for(i = 2; i <= size; i++){
        ReleaseParticle(particle_vec[i]);
        particle_vec[i] = new_particle_vec[i];
    }

    delete[] new_particle_vec;
    delete[] sample_index;
    delete[] weight_norm;
}

int SelectParticle(Particle** particle_vec, double* log_weight_vec, int size){
    double* weight = new double[size + 1];
    int i;
    double tmax;
    Lib::softmax(log_weight_vec, weight, size, tmax);
    double u = unif_rand();
    int result = size;
    double cumsum = 0;
    for(i = 1; i <= size; i++){
        cumsum += weight[i];
        if(u <= cumsum){
            result = i;
            break;
        }
    }
    return result;
}



void ReleaseParticle(Particle* particle){
    if(particle->thetree){
        particle->thetree->deall();
        delete particle->thetree;
    }
    delete particle;
}

bool CheckGrow(Particle** particle_vec, int size){
    Particle* cur;
    for(int i = 2; i <= size; i++){
        cur = *(particle_vec + i);
        if(cur->growable)
            return true;
    }
    return false;
}



