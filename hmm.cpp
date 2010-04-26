#include <stdlib.h>
#include "hmm.h"

using namespace std;

// HMM

HMM::HMM(const char* fn){

    TiXmlDocument doc(fn);
    doc.LoadFile();

    TiXmlElement* root = doc.RootElement();

    // Instantiate the start state
    int start;
    root->Attribute("start", &start);

    for( TiXmlElement* e = root->FirstChildElement(); e; e = e->NextSiblingElement() ){

        State *s = 0;
        const char* type = e->Attribute("type"); 

        if( type ){
            if( 0 == strcmp(type, "silent") ){
                s = new SilentState(e);
            } else if( 0 == strcmp(type, "accepting") ){
                s = new AcceptingState(e);
            } 
        } else {
            s = new State(e);
        }

        if( (int) _states.size() <= s->getId() ){
            _states.resize(s->getId()+1);
        }

        _states[s->getId()] = s;
    }

    setTransitions();
    _startState = _states[start];
}

void HMM::setTransitions(){

    vector< VState* >::iterator s_itr;

    for( s_itr = _states.begin(); s_itr != _states.end(); s_itr++ ){
            (*s_itr)->relabelTransition(_states);
    }


}

HMM::~HMM(){ }

char* HMM::generate(int request_length){

    char* res = (char*) calloc(request_length + 1, sizeof(char));
    //I don't need to set the null terminator on res, calloc does that for me.

    VState* s = _startState;
    double ep = 0.0;
    double tp = 0.0;
    int ii = 0;
    while( ii < request_length ){
        ep = _rng.rand();
        tp = _rng.rand();
        
        if( s->hasEmission() ){
            res[ii] = s->emit(ep);
            ii++;
        }

        if( s->hasTransition() ){
            s = s->transition(tp);
        } else {
            //If we don't have a valid transition, die
            return res;
        }
    }
    return res;
}

// HMM State
State::State(){};

State::State(TiXmlElement* stateElem){

    stateElem->Attribute("id", &_id);

    int ii = 0;
    for( TiXmlElement* e = stateElem->FirstChildElement(); e; e = e->NextSiblingElement() ){

        if( 0 == strcmp("transitions", e->Value()) ){
            int trans;
            if( e->Attribute("monomorphic", &trans) ){
                _transition = new MonoBehavior<VState*>((VState*) trans);
            } else {
                _transition = new PolyBehavior<VState*>(e->FirstChildElement());
            }
        } else if( 0 == strcmp("emissions", e->Value()) ){
            char emit;
            if( e->Attribute("monomorphic") ){
                emit = e->Attribute("monomorphic")[0];
                _emission = new MonoBehavior<char>(emit);
            } else {
                _emission = new PolyBehavior<char>(e->FirstChildElement());
            }
        } else {

            fprintf(stderr, "Got unexpected node type '%s'\n", e->Value());

        }

        ++ii;
    }
}

State::State(int id, char emission, int transition){
    _id = id;
    _emission = new MonoBehavior<char>(emission);
    _transition = new MonoBehavior<VState*>((VState*) transition);
}

State::~State() { 
    // delete _emission;
    // delete _transition;
}

VState* State::transition(double p){
    return _transition->emit(p); 
}

double State::transitionProbability(int state){
    return _transition->loglikelihood(state);
}

char State::emit(double p){
    return _emission->emit(p);
}

double State::emissionProbability(char em, int qual=INT_MIN){
    return _emission->loglikelihood(em, qual);
}

void State::relabelTransition(vector<VState*> &s){
   _transition->relabelTransition(s); 
}

// IndexedState
IndexedState::IndexedState(TiXmlElement *elem){

}

char IndexedState::emit(double p, int index){
    return _emissions->emit(p, index);
}

VState* IndexedState::transition(double p, int index){
    if( index < _emissions->size() ){
       return _internalTransition->emit(p, index);
    } else {
        return _terminalTransition->emit(p, 0);
    }
}

void IndexedState::relabelTransition(vector< VState* >& states){

   _internalTransition->relabelTransition(states);
   _terminalTransition->relabelTransition(states);

}
// HMM Behavior
template <class T>
double Behavior<T>::loglikelihood(bool match, double mutr, int qual){
    // err is a measure of sequencing error
    // _mutr measures mutational rate
    double seq_err = 0.0;

    //TODO What is the correct behavior here?
    // We haven't been given quality scores.
      if( match ){
        if( INT_MIN != qual ){
            seq_err = QUAL2LL(qual);
        } else {
            seq_err = 0.0;
        }
        return (log(1.0 - mutr) + seq_err);

    } else {

        if( INT_MIN != qual ){
            seq_err = QUAL2ERROR(qual);
        } else {
            seq_err = DBL_EPSILON;
        }

        double e = mutr * (1 - seq_err)
                 + (1 - mutr) * seq_err
                 + mutr * seq_err;

        return log(e);
    }
}

//      MonoBehavior
template <class T>
MonoBehavior<T>::MonoBehavior(T emission, double errorRate){
    _emission = emission;
    _mutr = errorRate;
}

template <class T>
MonoBehavior<T>::~MonoBehavior(){ }

template <class T>
T MonoBehavior<T>::emit(double n, int position){ return _emission; }

template <class T>
double MonoBehavior<T>::loglikelihood(T emit, int qual){
    return Behavior<T>::loglikelihood( emit == _emission, _mutr, qual);
}

template <class T>
void MonoBehavior<T>::relabelTransition(vector<T> &s){
    _emission = (T) s.front();
}

//      PolyBehavior
template <class T>
PolyBehavior<T>::PolyBehavior(TiXmlElement* e){
    double tally = 0.0;
    double p;
    T val;

    while( e ){

        e->QueryDoubleAttribute("prob", &p);
        
        // TODO This is kludgey
        if( e->Attribute("nval") ){
            val = (T) atoi(e->Attribute("nval"));
        } else {
            val = (T) e->Attribute("val")[0];
        }
        
        tally += p;

        e = e->NextSiblingElement();
        if( !e ){
            tally = 1.0;
        }
        _emissions.insert(pair<double, T>(tally, val));

        //TODO Is this actually what I want?
        _likelihoods.insert(pair<T, double>(val, log(p)));
    }
}


template <class T>
T PolyBehavior<T>::emit(double p, int position){
    return (*_emissions.lower_bound(p)).second;
}

template <class T>
double PolyBehavior<T>::loglikelihood(T emit, int qual){
    typename map<T, double>::iterator itr;

    itr = _likelihoods.find(emit);

    if( itr != _likelihoods.end() ){
        double p = (*itr).second;
        return log(p) + LOGERROR(qual);
    } else {
        // The probability of something else...
        double p = 1.0 - _density; // innate error probability
        double err = QUAL2ERROR(qual); // sequencing error rate
        return log((1 - p) * err + p * (1 - err) + (1 - p) * (1 - err));

    }

}

template <class T>
void PolyBehavior<T>::relabelTransition(vector<T> &s){

    typename map<double, T>::iterator e_itr;

    for( e_itr = _emissions.begin(); e_itr != _emissions.end(); e_itr++ ){
        unsigned int index = (intptr_t) (*e_itr).second;
        (*e_itr).second = (T) s[index];
    }
    
    _likelihoods.clear();



}



// IndexedState
template <class T>
IndexedBehavior<T>::IndexedBehavior(vector<T> emissions, double mutr){
    _emissions = emissions;
    _mutr = mutr;
}

template <class T>
IndexedBehavior<T>::~IndexedBehavior() { }

template <class T>
T IndexedBehavior<T>::emit(double p, int position){
    return _emissions[position];
}

template <class T>
double IndexedBehavior<T>::loglikelihood(T emission, int position, int qual){
    return Behavior<T>::loglikelihood( emission == _emissions[position], _mutr, qual);
}

// AcceptingState
AcceptingState::AcceptingState(TiXmlElement *e){

}


// SilentState
SilentState::SilentState(TiXmlElement *e){

}

SilentState::SilentState(int id, int successor){
    _id = id;
    _transition = new MonoBehavior<VState*>((VState*) successor);  
}

int main(){


    HMM h("hmm.xml");
    
    char* n = h.generate(100);
    printf("%s\n", n);
    free(n);
    return 0;
}
