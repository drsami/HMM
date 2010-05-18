#include <stdio.h>
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

        VState *s = 0;
        const char* type = e->Attribute("type"); 

        if( type ){
            if( 0 == strcmp(type, "indexed") ){
                s = new IndexedState(e);
            } else if( 0 == strcmp(type, "silent") ){
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
        //printf("Relabel state %d...\n", (*s_itr)->getId());
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
    int position = 0;
    while( ii < request_length ){
        ep = _rng.rand();
        tp = _rng.rand();

        if( s->hasEmission() ){
            res[ii] = s->emit(ep, position);
            ii++;
        }

        if( s->hasTransition() ){
            s = s->transition(tp, position);
        } else {
            //If we don't have a valid transition, die
            return res;
        }
    }
    return res;
}

void HMM::viterbi(char *seq, char *qual){

    int len = strlen(seq); 
    int numSearched = 0;
    SearchQueue dijkstraQueue;
    list<vsearch_entry<VState*>* > expandedNodes;
    set<pair<int, pair<int, int> > > searchedNodes;


    vsearch_entry<VState*>* lastNode = NULL;
    vsearch_entry<VState*> *head = new vsearch_entry<VState*>();
    head->state = _startState;
    head->incoming = NULL;
    head->position = 0;
    head->emission = 0;
    head->loglikelihood.v = 0.0;

    dijkstraQueue.push(head);

    while( !dijkstraQueue.empty() ){

        numSearched++;
        vsearch<VState*> node_wrapper = dijkstraQueue.top();
        vsearch_entry<VState*> *node = node_wrapper.getEntry();
        expandedNodes.push_back(node);
        dijkstraQueue.pop();

        pair<int, pair<int, int> > ndid(node->state->getId(), pair<int,int>(node->emission, node->position));

        if( searchedNodes.count(ndid) ){
            continue;
        } else {
            searchedNodes.insert(ndid);
        }

        printf("<%d, %d, %d>: %e [%c]\n", node->state->getId(), node->emission, node->position, node->loglikelihood.v, seq[node->emission]);

        list<string> labels;

        if( node->emission >= len ){
            printf("Last node: %d, %d, %d\n", node->state->getId(), node->emission, node->position);
            printf("Searched: %d\tQueue size: %d\n", numSearched, dijkstraQueue.size());
        
            do{
                VState *st = node->state;
                string lb = st->getLabel();
                if( lb != "" && (labels.empty() || labels.front() != lb)){
                    labels.push_front(lb);
                }
            } while( node->incoming && (node = node->incoming) );

            list<string>::iterator lb_itr;
            for( lb_itr = labels.begin() ; lb_itr != labels.end() ; lb_itr++ ){
                printf("%s\n", (*lb_itr).c_str());
            }


            list<vsearch_entry<VState*>* >::iterator sl_itr = expandedNodes.begin();

            return;
        }

        logdouble ep;
        ep.v = 0.0;

        if( qual ){
            ep = node->state->emissionProbability(seq[node->emission], node->position, qual[node->emission]);
        } else {
            ep = node->state->emissionProbability(seq[node->emission], node->position);
        }

        if( node->state->incrementing() ){
            node->emission += 1;
        }

        if( node->state->resetting() ){
            node->position = 0;
        }

        node->loglikelihood = node->loglikelihood + ep;
        node->state->enqueueTransitions(dijkstraQueue, node);
        lastNode = node;
    }
}

// HMM State
State::State(){};

State::State(TiXmlElement* stateElem){

    stateElem->Attribute("id", &_id);
    const char* label = stateElem->Attribute("label");
    if( label ){
        _label = string(label);
    } else {
        _label = "";
    }

    const char* reset = stateElem->Attribute("NoReset");
    const char* increment = stateElem->Attribute("NoIncrement");

    if( reset ){
        _positionReset = false;
    } else {
        _positionReset = true;
    }

    if( increment ){
        _positionIncrement = false;
    } else {
        _positionIncrement = true;
    }

    int ii = 0;
    for( TiXmlElement* e = stateElem->FirstChildElement(); e; e = e->NextSiblingElement() ){

        double callProb;
        if( 0 == strcmp("transitions", e->Value()) ){
            int trans;
            if( e->Attribute("monomorphic", &trans) ){
                e->Attribute("prob", &callProb);
                _transition = new MonoBehavior<VState*>((VState*) trans, callProb);
            } else {
                _transition = new PolyBehavior<VState*>(e->FirstChildElement());
            }
        } else if( 0 == strcmp("emissions", e->Value()) ){
            char emit;
            if( e->Attribute("monomorphic") ){
                emit = e->Attribute("monomorphic")[0];
                e->Attribute("prob", &callProb);
                _emission = new MonoBehavior<char>(emit, callProb);
            } else {
                _emission = new PolyBehavior<char>(e->FirstChildElement());
            }
        } else {
            fprintf(stderr, "Got unexpected node type '%s'\n", e->Value());
        }

        ++ii;
    }
}

State::~State() { 
    // delete _emission;
    // delete _transition;
}

VState* State::transition(double p, int &n){
    if( _positionReset ){ n = 0; }
    if( _positionIncrement ){ n = n + 1; }
    return _transition->emit(p); 
}

logdouble State::transitionProbability(VState* state, int position){
    return _transition->loglikelihood(state);
}

char State::emit(double p, int n){
    return _emission->emit(p);
}

logdouble State::emissionProbability(char em, int position, int qual){
    if( hasEmission() ){
        return _emission->loglikelihood(em, qual);
    } else {
        logdouble r;
        r.v = 0.0;
        return r;
    }
}

void State::relabelTransition(vector<VState*> &s){
    if( hasTransition() ){
        _transition->relabelTransition(s); 
    }
}

void State::enqueueTransitions(SearchQueue &searchQueue, vsearch_entry<VState*> *incomingNode){
    _transition->enqueueBehavior(searchQueue, incomingNode, _positionReset, _positionIncrement, hasEmission() ); 
}

// IndexedState
IndexedState::IndexedState(TiXmlElement *elem){

    elem->Attribute("id", &_id);
    const char* label = elem->Attribute("label");
    if( label ){
        _label = string(label);
    } else {
        _label = "";
    }
    int ii = 0;
    for( TiXmlElement* e = elem->FirstChildElement(); e; e = e->NextSiblingElement() ){
        double callProb;
        if( 0 == strcmp("internalTransition", e->Value()) ){
            int trans;
            if( e->Attribute("monomorphic", &trans) ){
                e->Attribute("prob", &callProb);
                _internalTransition = new MonoBehavior<VState*>((VState*) trans, callProb);
            } else {
                _internalTransition = new PolyBehavior<VState*>(e->FirstChildElement());
            }
        } else if( 0 == strcmp("terminalTransition", e->Value()) ){
            int trans;
            if( e->Attribute("monomorphic", &trans) ){
                e->Attribute("prob", &callProb);
                _terminalTransition = new MonoBehavior<VState*>((VState*) trans, callProb);
            } else {
                _terminalTransition = new PolyBehavior<VState*>(e->FirstChildElement());
            }
        } else if( 0 == strcmp("emissions", e->Value()) ){
            _emissions = new IndexedBehavior<char>(e);
        } else {
            fprintf(stderr, "Got unexpected node type '%s'\n", e->Value());
        }

        ++ii;
    }

}

char IndexedState::emit(double p, int index){
    return _emissions->emit(p, index);
}

VState* IndexedState::transition(double p, int &index){
    if( index < _emissions->size() - 1 ){
        VState* rtrn = _internalTransition->emit(p, index);
        ++index;
        return rtrn;
    } else {
        index = 0;
        return _terminalTransition->emit(p, 0);
    }
}


void IndexedState::enqueueTransitions(SearchQueue &searchQueue, vsearch_entry<VState*> *incomingNode){ 
    if( incomingNode->position < _emissions->size() - 1 ){
        _internalTransition->enqueueBehavior(searchQueue, incomingNode, false, true, hasEmission());
    } else {
        _terminalTransition->enqueueBehavior(searchQueue, incomingNode, true, false, hasEmission());
    }
}

logdouble IndexedState::emissionProbability(char emission, int position, int quality){

    if( position < _emissions->size() - 1 ){
        return _emissions->loglikelihood(emission, position, quality);
    } else {
        logdouble r;
        r.v = 0.0;
        return r;
    }
}

logdouble IndexedState::transitionProbability(VState* state, int position){
    if( position < _emissions->size() - 1 ){
        return _internalTransition->loglikelihood(state);
    } else {
        return _terminalTransition->loglikelihood(state);
    }

}

void IndexedState::relabelTransition(vector< VState* >& states){

    _internalTransition->relabelTransition(states);
    _terminalTransition->relabelTransition(states);

}

// HMM Behavior
template <class T>
logdouble Behavior<T>::loglikelihood(bool match, double prob, int qual){
    logdouble res;
    // err is a measure of sequencing error
    // _mutr measures mutational rate

    //TODO What is the correct behavior here?
    // We haven't been given quality scores.
    if( match ){
        assert(prob > 0.0);
        res.v = log(prob);
        return res;
    } else {
        assert(prob <= 1.0);
        res.v = log(1.0 - prob);
        return res;
    }

}

//      MonoBehavior
template <class T>
MonoBehavior<T>::MonoBehavior(T emission, double errorRate){
    _emission = emission;
    _prob = errorRate;
    _likelihood.v = log(_prob);
    _notlikelihood.v = log(1.0 - _prob);
    assert(_prob > 0.0 && _prob <= 1.0);
}

template <class T>
MonoBehavior<T>::~MonoBehavior(){ }

template <class T>
T MonoBehavior<T>::emit(double n, int position){ return _emission; }


template <class T>
logdouble MonoBehavior<T>::loglikelihood(T emit, int qual){
    return Behavior<T>::loglikelihood(emit == _emission, _prob, qual);
}

template <class T>
void MonoBehavior<T>::enqueueBehavior(priority_queue<vsearch<T> > &searchQueue, vsearch<T> entryWrapper, bool reset, bool increment, bool silent){

    vsearch_entry<T> *entry = entryWrapper.getEntry();
    vsearch_entry<T> *newBehavior = new vsearch_entry<T>();
    newBehavior->state = _emission;
    newBehavior->incoming = entry;
    if( reset ){
        newBehavior->position = 0;
    } else if( increment ){
        newBehavior->position = entry->position + 1;
    }

    logdouble r = entry->loglikelihood;
    r.v += log(_prob);
    newBehavior->loglikelihood = r;

    if( silent ){
        newBehavior->emission = entry->emission;
    } else {
        newBehavior->emission = entry->emission + 1;
    }

    printf("Adding a transition to node %d with position %d and emission %d\n", ((VState*) _emission)->getId(), newBehavior->position, newBehavior->emission);
    vsearch<T> e(newBehavior);
    searchQueue.push(e);
}

template <class T>
void MonoBehavior<T>::relabelTransition(vector<T> &s){
    _emission = (T) s[(intptr_t)_emission];
}

//      PolyBehavior
template <class T>
PolyBehavior<T>::PolyBehavior(TiXmlElement* e){
    double tally = 0.0;
    double p = 0.0;
    T val;

    while( e ){

        assert( TIXML_SUCCESS == e->QueryDoubleAttribute("prob", &p) );

        // TODO This is kludgey
        if( e->Attribute("nval") ){
            val = (T) atoi(e->Attribute("nval"));
        } else {
            val = (T) e->Attribute("val")[0];
        }

        tally += p;
        _emissions.insert(pair<double, T>(tally, val));
        //TODO Is this actually what I want?
        logdouble r;
        r.v = log(p);
        _likelihoods.insert(pair<T, logdouble>(val, r));
        e = e->NextSiblingElement();
    }

    _density = tally;
    // ensure that density is within proper bounds,
    // account for accumulated float precision loss
    assert( _density > 0.0 && _density <= 1.0000001  );
}

template <class T>
T PolyBehavior<T>::emit(double p, int position){
    return (*_emissions.lower_bound(p)).second;
}

template <class T>
logdouble PolyBehavior<T>::loglikelihood(T emit, int qual){
    typename map<T, logdouble>::iterator itr;

    itr = _likelihoods.find(emit);

    if( itr != _likelihoods.end() ){
        logdouble p = (*itr).second;
        return p; //+ LOGERROR(qual);
    } else {
        // The probability of something else...
        double p = 1.0 - _density; // innate error probability
        assert(p > 0.0);
        logdouble r;
        r.v = log(p);
        return r;
    }

}

template <class T>
void PolyBehavior<T>::enqueueBehavior(priority_queue<vsearch<T> > &searchQueue, vsearch<T> entryWrapper, bool reset, bool increment, bool silent){

    vsearch_entry<T>* entry = entryWrapper.getEntry();
    typename map<T, logdouble>::iterator l_itr;
    for( l_itr = _likelihoods.begin(); l_itr != _likelihoods.end(); l_itr++ ){
        vsearch_entry<T> *newBehavior = new vsearch_entry<T>();
        newBehavior->incoming = entry;
        if( reset ){
            newBehavior->position = 0;
        } else if( increment ){
            newBehavior->position = entry->position + 1;
        }

        if( silent ){
            newBehavior->emission = entry->emission;
        } else {
            newBehavior->emission = entry->emission + 1;
        }

        newBehavior->state = l_itr->first;
        logdouble r;
        r = entry->loglikelihood;
        r = r + l_itr->second;
        newBehavior->loglikelihood = r;
        vsearch<T> e(newBehavior);
        searchQueue.push(e);
    }
}

template <class T>
void PolyBehavior<T>::relabelTransition(vector<T> &s){

    typename map<double, T>::iterator e_itr;

    map<T, logdouble> newLikelihoods;

    for( e_itr = _emissions.begin(); e_itr != _emissions.end(); e_itr++ ){
        logdouble ll = _likelihoods.find(e_itr->second)->second;
        int ptr = (intptr_t) e_itr->second;
        e_itr->second = (T) s[ptr];

        newLikelihoods.insert(pair<T, logdouble>((T) s[ptr], ll));
    }

    _likelihoods.clear();
    _likelihoods = newLikelihoods;
}

// IndexedState
template <class T>
IndexedBehavior<T>::IndexedBehavior(TiXmlElement *elem){

    const char* c = elem->Attribute("str");
    string s(c);
    //This is dangerous for non string types.
    _emissions = vector<T>(s.begin(), s.end());
    assert( TIXML_SUCCESS == elem->QueryDoubleAttribute("prob", &_prob) );
    _likelihood.v = log(_prob);
    _notlikelihood.v = log(1.0 - _prob);
    assert(_prob > 0.0 && _prob <= 1.000001);
}

template <class T>
IndexedBehavior<T>::~IndexedBehavior() { }

template <class T>
T IndexedBehavior<T>::emit(double p, int position){
    return _emissions[position];
}

template <class T>
logdouble IndexedBehavior<T>::loglikelihood(T emission, int position, int qual){
    bool match = (emission == _emissions[position]);
    return Behavior<T>::loglikelihood(match, _prob, qual);
}

template <class T>
void IndexedBehavior<T>::enqueueBehavior(priority_queue<vsearch<T> > &s, vsearch<T> n, bool reset, bool increment, bool silent){
    return;
}

// AcceptingState
AcceptingState::AcceptingState(TiXmlElement *e) : SilentState(e) {
}


// SilentState
SilentState::SilentState(TiXmlElement *e) : State(e) {
}

int main(int argc, char* argv[]){

    if( argc == 1 ){
        fprintf(stderr, "Usage: %s <model.xml>\n", argv[0]);
        return 1;
    }

    const char* fp = argv[1];
    HMM h(fp);


    char* n = h.generate(320);
    printf("%s\n", n);
    h.viterbi(n);
    //printf("%s\n", n);
    free(n);
    return 0;
}
