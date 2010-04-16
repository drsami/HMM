#include <stdlib.h>
#include "hmm.h"

using namespace std;

// HMM

HMM::HMM(const char* fn){

    TiXmlDocument doc(fn);
    doc.LoadFile();

    TiXmlElement* root = doc.RootElement();

    // Instantiate the start state
    root->Attribute("start", &_startState);

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

    printf("%s\n", generate(50));
}

/*
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
 */
HMM::~HMM(){ }

char* HMM::generate(int request_length){

    char* res = (char*) calloc(request_length + 1, sizeof(char));
    //I don't need to set the null terminator on res, calloc does that for me.

    State* s = getState(_startState);
    double p = 0.0;
    int ii = 0;
    while( ii < request_length ){
        if( s->hasEmission() ){
            res[ii] = s->emit(p, ii);
            ii++;
        }
        if( s->hasTransition() ){
            s = getState(s->transition(p, ii));
        } else {
            //If we don't have a valid transition, die
            return res;
        }
    }
    return res;
}

State* HMM::getState(int id){
    if( id >= 0 && id < (int) _states.size() ){
        return _states[id];
    }

    raise(SIGSEGV);
    return (State*) -1;
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
        vmatrix[_startState][0].best_score = start->emissionProbability(emissions[0],qualities[0], 0); 
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
            prb += vmatrix[e.out_state][e.out_emit].best_score; 
        }

        // Add the transition cost to this node
        prb += e.cost;

        // Add the emission cost
        if( s->hasEmission() ){
            prb += s->emissionProbability(emissions[e.in_emit],qualities[e.in_emit], e.in_emit);
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

State::State(TiXmlElement* stateElem){

    stateElem->Attribute("id", &_id);

    int ii = 0;
    for( TiXmlElement* e = stateElem->FirstChildElement(); e; e = e->NextSiblingElement() ){

        if( 0 == strcmp("transitions", e->Value()) ){
            int trans;
            if( e->Attribute("monomorphic", &trans) ){
                _transition = new MonoBehavior<int>(trans);
            } else {
                _transition = new PolyBehavior<int>(e->FirstChildElement());
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

int State::transition(double p, int n){
    return _transition->emit(p, n); 
}

double State::transitionProbability(int state, int n){
    return _transition->loglikelihood(state, 20, n);
}

char State::emit(double p, int n){
    return _emission->emit(p, n);
}

double State::emissionProbability(char em, int n, int qual=0){
    return _emission->loglikelihood(em, n, qual);
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
T MonoBehavior<T>::emit(double n, int position){ return _emission; }

template <class T>
double MonoBehavior<T>::loglikelihood(T emit, int n, int qual=20){
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
PolyBehavior<T>::PolyBehavior(TiXmlElement* e){
    double tally = 0.0;
    double p;
    T val;

    while( e ){

        e->QueryDoubleAttribute("prob", &p);
        val = (T) (e->Attribute("val"))[0];

        tally += p;

        e = e->NextSiblingElement();
        if( !e ){
            tally = 1.0;
        }
        _emissions.insert(pair<double, T>(tally, (T) val));

        //TODO Is this actually what I want?
        _likelihoods.push_back(pair<double, T>(log(p), (T) val));
    }
}


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
T PolyBehavior<T>::emit(double p, int position){
    return (*_emissions.lower_bound(p)).second;
}

template <class T>
double PolyBehavior<T>::loglikelihood(T emit, int position, int qual=20){
    raise(SIGSEGV);
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

// AcceptingState
AcceptingState::AcceptingState(TiXmlElement *e){

}


// SilentState
SilentState::SilentState(TiXmlElement *e){

}

SilentState::SilentState(int id, int successor){
    _id = id;
    _transition = new MonoBehavior<int>(successor);  
}

int main(){


    HMM h("hmm.xml");
    //    printf("%s\n", h.generate(1000));
    //    return 0;
}
