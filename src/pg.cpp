#include "pg.h"



void InitCachetemp(Cachetemp* cachetemp){

}







void Initparticles(Particle** particle_vec, int len){
    for(int i = 0; i < len; i++){
        Particle* cur_particle = *(particle_vec + i);
        cur_particle->thetree = new Node;
        cur_particle->equeue.append(cur_particle->thetree);

        //todo add the weight
    }
}

void Setparticlebytree(Particle* particle, Node* oldtree);

void Settreebyparticle(Node* dsttree, Particle* particle);

void Releaseparticle(Particle* particle){
    if(!(particle->thetree))
        delete particle->thetree;
}
