#include <stdlib.h>
#include "hmm.h"

using namespace std;

// HMM

HMM::HMM(const char* fn){
    _rng = MTRand(); 
    TiXmlDocument doc(fn);
    doc.LoadFile();

    TiXmlElement* root = doc.RootElement();

    // Instantiate the start state
    int start = INT_MIN;
    root->Attribute("start", &start);
    assert( start != INT_MIN );

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

HMM::~HMM(){
    vector<VState*>::iterator sv_itr;
    for( sv_itr = _states.begin(); sv_itr != _states.end(); ){
        VState* s = (*sv_itr);
        sv_itr++;
        delete s;
    }
}

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
    priority_queue<vsearch<VState*> > dijkstraQueue;
    list<vsearch_entry<VState*>* > allNodes;

    vsearch_entry<VState*> *head = new vsearch_entry<VState*>();
    head->state = _startState;
    head->incoming = NULL;
    head->position = 0;
    head->emission = 0;
    head->loglikelihood = 0.0;

    dijkstraQueue.push(vsearch<VState*>(head));

    while( !dijkstraQueue.empty() ){

        numSearched++;
        vsearch<VState*> search = dijkstraQueue.top();
        vsearch_entry<VState*> *node = search.get();
        dijkstraQueue.pop();

        //printf("<%d, %d, %d>: %e\n", node.state->getId(), node.emission, node.position, node.loglikelihood);

        if( node->emission >= len ){
            printf("Last node: %d, %d, %d\n", node->state->getId(), node->emission, node->position);
            printf("Searched: %d\tQueue size: %d\n", numSearched, dijkstraQueue.size());

            list<string> labels;

            do{
            
                string str = node->state->getLabel();
                if( !str.empty() && str != labels.front() ){
                    labels.push_front(str);
                }

            } while( node->incoming && (node = node->incoming) );


            list<string>::iterator label_itr;

            for( label_itr = labels.begin(); label_itr != labels.end(); label_itr++ ){
                printf("%s\n", (*label_itr).c_str());
            }

            /* Throw everything out */
            list<vsearch_entry<VState*>* >::iterator all_itr = allNodes.begin();
            for( ; all_itr != allNodes.end(); ){
                vsearch_entry<VState*> *delnode = (*all_itr);
                all_itr++;
                delete delnode;
            }

            return;
        }

        double ep = 0.0;
        if( qual ){
            ep = node->state->emissionProbability(seq[node->emission], node->position, qual[node->emission]);
        } else {
            //printf("Position: %d\n", node->position);
            ep = node->state->emissionProbability(seq[node->emission], node->position);
        }

        if( node->state->incrementing() ){
            node->emission += 1;
        }

        if( node->state->resetting() ){
            node->position = 0;
        }

        node->loglikelihood = node->loglikelihood + ep;
        //node->state->enqueueTransitions(allNodes, dijkstraQueue, node);
    }
}

// HMM State
State::State(){};

State::State(TiXmlElement* stateElem){

    stateElem->Attribute("id", &_id);

    _emission = NULL;
    _transition = NULL;

    const char* reset = stateElem->Attribute("NoReset");
    const char* increment = stateElem->Attribute("Increment");

    const char* label = stateElem->Attribute("label");

    if( label ){
        _label = string(label);
    }

    if( reset ){
        _positionReset = false;
    } else {
        _positionReset = true;
    }

    if( increment ){
        _positionIncrement = true;
    } else {
        _positionIncrement = false;
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
    if( _emission ){ delete _emission;_emission = NULL; }
    if( _transition ){ delete _transition; _transition = NULL; }
}

VState* State::transition(double p, int &n){
    if( _positionReset ){ n = 0; }
    if( _positionIncrement ){ n = n + 1; }
    return _transition->emit(p); 
}

double State::transitionProbability(VState* state, int position){
    return _transition->loglikelihood(state);
}

char State::emit(double p, int n){
    return _emission->emit(p);
}

double State::emissionProbability(char em, int position, int qual){
    if( hasEmission() ){
        return _emission->loglikelihood(em, qual);
    } else {
        return 0.0;
    }
}

void State::relabelTransition(vector<VState*> &s){
    if( hasTransition() ){
        _transition->relabelTransition(s); 
    }
}

void State::enqueueTransitions(
        list<vsearch_entry<VState*>*> &allNodes,
        priority_queue<vsearch_entry<VState*>* > &searchQueue, 
        vsearch_entry<VState*> *incomingNode){
    _transition->enqueueBehavior(allNodes, searchQueue, incomingNode, _positionReset, _positionIncrement, hasEmission() ); 
}

// IndexedState
IndexedState::IndexedState(TiXmlElement *elem){

    elem->Attribute("id", &_id);

    const char* label = elem->Attribute("label");

    if( label ){
        _label = string(label);
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
                //e->Attribute("prob", &callProb);
                callProb = 1.0;
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

IndexedState::~IndexedState(){
    delete _emissions;
    delete _internalTransition;
    delete _terminalTransition;
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


void IndexedState::enqueueTransitions(
        list<vsearch_entry<VState*>*> &allNodes,
        priority_queue<vsearch_entry<VState*>* > &searchQueue, 
        vsearch_entry<VState*> *incomingNode){ 
    if( incomingNode->position < _emissions->size() - 1 ){
        _internalTransition->enqueueBehavior(allNodes, searchQueue, incomingNode, false, true, hasEmission());
    } else {
        _terminalTransition->enqueueBehavior(allNodes, searchQueue, incomingNode, true, false, hasEmission());
    }
}

double IndexedState::transitionProbability(VState* state, int position){
    if( position < _emissions->size() - 1){
        return _internalTransition->loglikelihood(state);
    } else {
        return _terminalTransition->loglikelihood(state);
    }
}

double IndexedState::emissionProbability(char em, int position, int qual){
    if( hasEmission() ){
        return _emissions->loglikelihood(em, position, qual);
    } else {
        return 0.0;
    }
}


void IndexedState::relabelTransition(vector< VState* >& states){

    _internalTransition->relabelTransition(states);
    _terminalTransition->relabelTransition(states);

}

// HMM Behavior
template <class T>
double Behavior<T>::LogLikelihood(bool match, double prob, int qual){
    // err is a measure of sequencing error
    // _mutr measures mutational rate

    //TODO What is the correct behavior here?
    // We haven't been given quality scores.
    if( match ){
        return log(prob);
    } else {
        return log(1.0 - prob);
    }

}

//      MonoBehavior
template <class T>
MonoBehavior<T>::MonoBehavior(T emission, double errorRate){
    _emission = emission;
    _prob = errorRate;
    _likelihood = log(_prob);
    _notlikelihood = log(1.0 - _prob);
    assert(_prob > 0 && _prob <= 1.0);
}

template <class T>
MonoBehavior<T>::~MonoBehavior(){ }

template <class T>
T MonoBehavior<T>::emit(double n, int position){ return _emission; }


template <class T>
double MonoBehavior<T>::loglikelihood(T emit, int qual){
    return Behavior<T>::LogLikelihood(emit == _emission, _prob, qual);
}

template <class T>
void MonoBehavior<T>::enqueueBehavior(
        list<vsearch_entry<T>*> &all_nodes,
        priority_queue<vsearch_entry<T>* > &searchQueue, 
        vsearch_entry<T> *entry, bool reset, bool increment, bool silent){

    vsearch_entry<T> *newBehavior = new vsearch_entry<T>();
    newBehavior->state = _emission;
    newBehavior->incoming = entry;
    if( reset ){
        newBehavior->position = 0;
    } else if( increment ){
        newBehavior->position = entry->position + 1;
    }

    newBehavior->loglikelihood = entry->loglikelihood + _prob;

    if( silent ){
        newBehavior->emission = entry->emission;
    } else {
        newBehavior->emission = entry->emission + 1;
    }
    all_nodes.push_back(newBehavior);
    searchQueue.push(newBehavior);
}

template <class T>
void MonoBehavior<T>::relabelTransition(vector<T> &s){
    _emission = (T) s[(intptr_t)_emission];
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
        _likelihoods.insert(pair<T, double>(val, p));
    }
}

template <class T>
PolyBehavior<T>::~PolyBehavior(){
    _emissions.clear();
    _likelihoods.clear();
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
        return log(p); //+ LOGERROR(qual);
    } else {
        // The probability of something else...
        double p = 1.0 - _density; // innate error probability
        return log(p);
    }

}

template <class T>
void PolyBehavior<T>::enqueueBehavior(
        list<vsearch_entry<T>*> &all_nodes,
        priority_queue<vsearch_entry<T>* > &searchQueue, 
        vsearch_entry<T> *entry, bool reset, bool increment, bool silent){

    typename map<T, double>::iterator l_itr;
    for( l_itr = _likelihoods.begin(); l_itr != _likelihoods.end(); l_itr++ ){
        vsearch_entry<T> *newBehavior = new vsearch_entry<T>();
        newBehavior->incoming = entry;
        if( reset ){
            newBehavior->position = 0;
        } else if( increment ){
            newBehavior->position = entry->position + 1;
        } else {
            newBehavior->position = entry->position;
        }

        if( silent ){
            newBehavior->emission = entry->emission;
        } else {
            newBehavior->emission = entry->emission + 1;
        }


        newBehavior->state = l_itr->first;
        newBehavior->loglikelihood = entry->loglikelihood + log(l_itr->second);
        all_nodes.push_back(newBehavior);
        searchQueue.push(newBehavior);
    }
}

template <class T>
void PolyBehavior<T>::relabelTransition(vector<T> &s){

    typename map<double, T>::iterator e_itr;

    map<T, double> newLikelihoods;

    for( e_itr = _emissions.begin(); e_itr != _emissions.end(); e_itr++ ){
        double ll = _likelihoods.find(e_itr->second)->second;
        int ptr = (intptr_t) e_itr->second;
        e_itr->second = (T) s[ptr];

        newLikelihoods.insert(pair<T, double>((T) s[ptr], ll));
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
    elem->Attribute("prob", &_prob);
    _likelihood = log(_prob);
    _notlikelihood = log(1.0 - _prob);
    assert(_prob > 0.0 && _prob <= 1.0);
}

template <class T>
IndexedBehavior<T>::~IndexedBehavior() { }

template <class T>
T IndexedBehavior<T>::emit(double p, int position){
    return _emissions[position];
}

template <class T>
double IndexedBehavior<T>::loglikelihood(T emission, int position, int qual){
    bool match = emission == _emissions[position];
    return Behavior<T>::LogLikelihood(match, _prob, qual);
}

template <class T>
void IndexedBehavior<T>::enqueueBehavior(
        list<vsearch_entry<T>*> &all_nodes,
        priority_queue<vsearch_entry<T>* > &s, 
        vsearch_entry<T> *n, bool reset, bool increment, bool silent){
    return;
}

// AcceptingState
AcceptingState::AcceptingState(TiXmlElement *e) : SilentState(e) {
}


// SilentState
SilentState::SilentState(TiXmlElement *e) : State(e) {
}

/*
int main(){

    HMM *h = new HMM("hmm.xml");

    char* n = h->generate(300);
    h->viterbi(n);
    printf("%s\n", n);
    free(n);
    delete h;
    return 0;
}
*/
