#include <stdlib.h>
#include "hmm.h"

using namespace std;

// HMM
HMM::HMM(char* sequence){

    int len = strlen(sequence);
    _startState = 0;
    _states.reserve(len+2);
    _states.push_back(new SilentState(0, 1));

    for(int ii = 0; ii < len; ++ii){
       State *s = new State(ii + 1, sequence[ii],ii+2);         
       _states.push_back(s);
    }

    _states.push_back(new AcceptingState(len+1));


    vector< State >::iterator s_itr;
}

HMM::~HMM(){ }

char* HMM::generate(int request_length){
    
    char* res = (char*) calloc(request_length + 1, sizeof(char));
    //I don't need to set the null terminator on res, calloc does that for me.

    State* s = getState(_startState);
    double p = 0.0;
    int ii = 0;
    while( ii < request_length ){
        if( s->hasEmission() ){
           res[ii] = s->emit(p);
           ii++;
        }
        if( s->hasTransition() ){
            s = getState(s->transition(p));
        } else {
            //If we don't have a valid transition, die
            return res;
        }
    }
    return res;
}

State* HMM::getState(int id){
    if( id >= 0 && id < _states.size() ){
        return _states[id];
    }
}

list< int > HMM::viterbi( char* emissions, char* qualities ){
    const int len = strlen(emissions);
    const int noStates = _states.size();
    vcell vmatrix[noStates][len];

    for( int ii = 0; ii < noStates; ii++ ){
        for( int jj = 0; jj < len; jj++ ){
            vmatrix[ii][jj].predecessor = -1;
            vmatrix[ii][jj].best_score = -DBL_MAX;
        }
    }

    queue< edge > searchQueue;

    //TODO Add support for non-terminal silent states
    edge startEdge;
        startEdge.out_state = _startState;
        startEdge.out_emit = 0;
        startEdge.in_state = -1;
        startEdge.in_emit = -1;
        startEdge.cost = 0.0;
    // We always start in the start state
    vmatrix[_startState][0].predecessor = -1;
    vmatrix[_startState][0].emission = -1;
    vmatrix[_startState][0].best_score = 0.0;

    State* start = getState(_startState);
    if( start->hasEmission() ){
       vmatrix[_startState][0].best_score = start->emissionProbability(emissions[0],qualities[0]); 
    }

    searchQueue.push(startEdge);

    while( !searchQueue.empty() ){
        // This version of viterbi is "foward looking" instead of "backward looking"
        // I can preserve my runtime guarantees because I traverse each edge only once
        edge e = searchQueue.front();
        searchQueue.pop();
        State* s = getState(e.in_state);

        s->enqueueTransitions(e.out_state, searchQueue);

        double prb = 0.0;

        // Get the probability of reaching the predecessor
        if( e.out_state >= 0 && e.out_emit >= 0 ){
            vmatrix[e.out_state][e.out_emit].best_score; 
        }

        // Add the transition cost to this node
        prb += e.cost;

        // Add the emission cost
        if( s->hasEmission() ){
            prb += s->emissionProbability(emissions[e.in_emit],qualities[e.in_emit]);
        }

        if( vmatrix[e.in_state][e.in_emit].best_score < prb ){
            vmatrix[e.in_state][e.in_emit].best_score = prb;
            vmatrix[e.in_state][e.in_emit].predecessor = e.out_state;
            vmatrix[e.in_state][e.in_emit].emission = e.out_emit;
        }

    }

    /*
     * Backtrace the Viterbi matrix
     */

    int pred = -1;
    int em = -1;
    double score = -DBL_MAX;
    
    int emid = len - 1;
    list<int> res;

    // Find the best accepting state
    for(int ii = 0; ii < noStates; ii++ ){
        if( vmatrix[ii][emid].best_score > score ){
            score = vmatrix[ii][emid].best_score;
            pred = vmatrix[ii][emid].predecessor;
            em = vmatrix[ii][emid].emission;
        }
    }
   
    // And backtrace!
    while( em >= 0){
        res.push_back(pred);
        int p = vmatrix[pred][em].predecessor;
        int e = vmatrix[pred][em].emission;
        pred = p;
        em = e;
    }

    return res;
}

// HMM State
State::State(){};

State::State(int id, char emission, int transition){
    _id = id;
    _emission = new MonoBehavior<char>(emission);
    _transition = new MonoBehavior<int>(transition);
}

State::State(list<pair<double, char> > emission, list<pair<double, int> > transitions){
    
    if( emission.size() == 1 ){
        _emission = new MonoBehavior<char>(emission.front().second);
    } else {
        _emission = new PolyBehavior<char>(emission);
    }

    if( transitions.size() == 1 ){
        _transition = new MonoBehavior<int>(transitions.front().second);
    } else {
        _transition = new PolyBehavior<int>(transitions);
    }
}

State::~State() { 
    // delete _emission;
    // delete _transition;
}

int State::transition(double n){
   _transition->emit(n); 
}

double State::transitionProbability(int state){
    _transition->loglikelihood(state, 20);
}

char State::emit(double n){
    _emission->emit(n);
}

double State::emissionProbability(char em, int qual=0){
    _emission->loglikelihood(em, qual);
}

void State::enqueueTransitions(int em, queue< edge > &sq){
    int out_emit = em;
    
    if( hasEmission() ){
        out_emit += 1;
    }

    _transition->enqueueEmissions(_id, em, out_emit, sq);
    return;
}

// HMM Behavior
//      Base class
/*
template <class T>
T Behavior<T>::emit(double p){ return (T) -1; };

template <class T>
void Behavior<T>::enqueueEmissions(T a, int b, int c, std::queue< edge > &d){ return; }

template <class T>
double Behavior<T>::loglikelihood(T a, int b){ return 0.0; };
*/

//      MonoBehavior
template <class T>
MonoBehavior<T>::MonoBehavior(T emission){
    _emission = emission;
}

template <class T>
MonoBehavior<T>::~MonoBehavior(){ }

template <class T>
T MonoBehavior<T>::emit(double n){ return _emission; }

template <class T>
double MonoBehavior<T>::loglikelihood(T emit, int qual=20){
    if( emit == _emission ){
        return LOGQUALITY(qual);
    } else {
        return 0.0; // Return error probability here
    }
}

template <class T>
void MonoBehavior<T>::enqueueEmissions(T out_id, int out_em, int in_em, queue< edge > &sq){
    edge e;
        e.out_state = (int) out_id;
        e.out_emit = out_em;
        e.in_state = (int) _emission;
        e.in_emit = in_em;
        e.cost = 0.0; // deterministic edge, no traversal cost 
    sq.push(e);
    return;
}

//      PolyBehavior
template <class T>
PolyBehavior<T>::PolyBehavior(list<pair<double, T> > emissions){

    double p = 0.0;
    typename list< pair<double, T > >::iterator e_itr;
    e_itr = emissions.begin(); 
    
    while( e_itr != emissions.end() ){
        p = (*e_itr).first + p;
        T val = (*e_itr).second;
        if( ++e_itr == emissions.end() ){
            p = 1.0;
        }
     
        pair<double, T> lk(log((*e_itr).first), val);
        pair<double, T> pr(p, val);
        _likelihoods.push_back(lk);
        _emissions.insert(pr);
    }
}

template <class T>
T PolyBehavior<T>::emit(double p){
    return (*_emissions.upper_bound(p)).second;
}

template <class T>
double PolyBehavior<T>::loglikelihood(T emit, int qual=20){
    return 0.0;
}

template <class T>
void PolyBehavior<T>::enqueueEmissions(T out_id, int out_em, int in_em, queue< edge > &sq){
    edge e;
        e.out_state = (int) out_id;
        e.out_emit = out_em;
        e.in_emit = in_em;
        
   typename list< pair<double, T> >::iterator lk_itr;

        for( lk_itr = _likelihoods.begin(); lk_itr != _likelihoods.end(); lk_itr++ ){
            e.cost = (*lk_itr).first; // deterministic edge, no traversal cost 
            e.in_state = (int) (*lk_itr).second;
            sq.push(e);
        }
    return;
}

SilentState::SilentState(int id, int successor){
    _id = id;
    _transition = new MonoBehavior<int>(successor);  
}

int main(){


    HMM h("ATCGATCGATCG");
    printf("%s\n", h.generate(1000));
    return 0;
}
