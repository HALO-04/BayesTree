#include "track.h"

Tracker::Tracker(){
	this->cur = this->front = this->rear = new State;
	this->rear->next = NULL;
}

Tracker::~Tracker(){
	this->cur = this->front;
	State* temp;
	while(this->cur){
		temp = this->cur->next;
		delete this->cur;
		this->cur = temp;
	}
}

void Tracker::CopyFrom(Tracker* src){
	//clear dst tracker
	this->cur = this->front->next;
	State* temp;
	while(this->cur){
		temp = this->cur->next;
		delete this->cur;
		this->cur = temp;
	}
	this->rear = this->front;
	this->rear->next = NULL;

	//copy tracker
	temp = src->front->next;
	while(temp){
		this->append(temp->isgrow, temp->var, temp->split_idx);
		temp = temp->next;
	}
	this->reset();
}

void Tracker::append(bool misgrow, int mvar, int msplit_idx){
	State* new_state = new State;
	new_state->isgrow = misgrow;
	new_state->mvar = mvar;
	new_state->mvar = msplit_idx;
	new_state->next = NULL;

	this->rear->next = new_state;
	this->rear = new_state;
}

void Tracker::reset(){
	this->cur = this->front;
}

bool Tracker::gonext(){
	if(this->cur->next){
		this->cur = this->cur->next;
		return true;
	}else{
		this->reset();
		return false;
	}
}




