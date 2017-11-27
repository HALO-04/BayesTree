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

void Runsample(Node* thetree){

    Particle**  particle_vec = new Particle*[NumParticle + 1];
    double* log_weight_vec = new double[NumParticle + 1];
    Initparticles(particle_vec, log_weight_vec, NumParticle);

    Particle* first_particle = *(particle_vec + 1);
    first_particle->thetree->deall();
    thetree->CopyTree(first_particle->thetree);

    Node* gnode;

    while(true){
        for(int i = 2; i <= len; i++){
            Particle* cur_particle = *(particle_vec + i);

            bool done = Growparticle(cur_particle, &gnode);
            if(done){
                log_weight_vec[i] += UpdateWeight(gnode);
            }
        }
        Resample(particle_vec, log_weight_vec, NumParticle);
        if(!Checkgrow(particle_vec, NumParticle)) break;
    }

    int select_idx = SelectParticle(particle_vec, weight_vec, NumParticle);
    Particle* selected = *(particle_vec + select_idx);
    thetree->deall();
    selected->CopyTree(thetree);

    for(int i = 1; i <= NumParticle; i++) Releaseparticle(particle_vec + i);
    delete[] particle_vec;
    delete[] weight_vec;
}




void Initparticles(Particle** particle_vec, double* weight_vec, int len){
    Particle* cur_particle;
    for(int i = 1; i <= len; i++){
        cur_particle = *(particle_vec + i);
        cur_particle->thetree = new Node;
        cur_particle->thetree->SetData();

        cur_particle->equeue.append(cur_particle->thetree);
        cur_particle->growable = true;
    }

    double loglik = LogLNode(cur_particle->thetree);
    for(int i = 1; i <= len; i++) weight_vec[i] = loglik;

}

bool Growparticle(Particle* p, Node** pgrow_node){
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

    int LeftEx, RightEx;
    SpawnChildren(grow_node, &LeftEx, &RightEx);
    (*pgrow_node) = grow_node;
    return true;
}


bool DrValidSplit(Node* gnode){
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
            int* n_split = new int[ReluNum[tvar] + 1];
            Rprintf("unused now\n");
            exit(0);
            while(true){
                pIvec = (int*)cur_cell->contents;
                tmp_value = XDat[*pIvec][tvar];

                if(cur_cell->End)
                    break;
                cur_cell = cur_cell->after;
            }
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


int PGLowerBound(int *array, int size, double key){
    int first = 0, middle;
    int half, len;
    len = size;

    while(len > 0) {
        half = len >> 1;
        middle = first + half;
        if(array[middle] < key) {
            first = middle + 1;
            len = len-half-1;       //在右边子序列中查找
        }
        else
            len = half;            //在左边子序列（包含middle）中查找
    }
    return first + 1;
}

void Resample(Particle** particle_vec, double* log_weight_vec, int size){
    double* weight_norm = new double[size + 1];
    int* sample_index = new int[size + 1];
    double tmax = -1.0 * DBL_MAX;
    int i;
    double sum, log_pd;
    for(i = 1; i <= size; i++){
        weight_norm[i] = log_weight_vec[i];
        if(weight_norm[i] > tmax)
            tmax = weight_norm[i];
    }
    sum = 0;
    for(i = 1; i <= size; i++){
        weight_norm[i] -= tmax;
        weight_norm[i] = std::exp(weight_norm[i]);
        sum += weight_norm[i];
    }
    for(i = 1; i <= size; i++)
        weight_norm[i] /= sum;
    log_pd = std::log(sum) + tmax;
    double log_numP = std::log(size);
    for(i = 1; i <= size; i++)
        log_weight_vec[i] = log_pd - log_numP;

    Lib::SampleMultinomial(weight_norm, size, sample_index, size);
    Particle** new_particle_vec = new Particel*[size + 1];

    for(i = 2; i <= size; i++){
        Particle* new_particle = new Particle;
        new_particle->thetree = new Node;
        new_particle->CopyFrom(particle_vec[sample_index[i]]);
        new_particle_vec[i] = new_particle;
    }
    for(i = 2; i <= size; i++){
        delete particle_vec[i];
        particle_vec[i] = new_particle_vec[i];
    }

    delete[] new_particle_vec;
    delete[] sample_index;
    delete[] weight_norm;
}

int SelectParticle(Particle** particle_vec, double* log_weight_vec, int size){
    double* weight = new double[size + 1];
    int i;
    double max = -1.0 * DBL_MAX;
    for(i = 1; i <= size; i++){
        if(max < log_weight_vec[i])
            max = log_weight_vec[i];
    }
    double sum = 0;
    for(i = 1; i <= size; i++){
        weight[i] = log_weight_vec[i] - max;
        weight[i] = std::exp(weight[i]);
        sum += weight[i];
    }
    double u = unif_rand();
    int result = size;
    double cumsum = 0;
    for(i = 1; i <= size; i++){
        weight[i] /= sum;
        cumsum += weight[i];
        if(u <= cumsum){
            result = i;
            break;
        }
    }
    return result;
}



void Releaseparticle(Particle* particle){
    if(!(particle->thetree))
        delete particle->thetree;
    delete particle;
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



