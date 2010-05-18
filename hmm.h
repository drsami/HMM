#include "kseq.h"
#include "qual.h"

#include "tinyxml.h"
#include "MersenneTwister.h"

#include <float.h>
#include <math.h>
#include <signal.h>

#include <list>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

class HMM;
class VState;
class State;

typedef struct{
    double v;
} logdouble;

logdouble operator+(const logdouble lhs, const logdouble rhs){
    logdouble r;
    r.v = lhs.v + rhs.v;
    return r;
}

logdouble operator+(const logdouble lhs, const double rhs){
    logdouble r;
    r.v = lhs.v + log(rhs);
    return r;
}

bool operator<(const logdouble lhs, const logdouble rhs){
    return lhs.v < rhs.v;
}

template <class T>
class vsearch_entry {
    //TODO add getters and setters to vsearch_entry
    public:
        vsearch_entry<T>(){}
        T state;
        vsearch_entry<T>* incoming;
        int position;
        int emission;
        logdouble loglikelihood;
         bool operator<(const vsearch_entry<T> n) const {
            return loglikelihood < n.loglikelihood;
        }
};

template <class T>
class vsearch {
    public:
        vsearch(vsearch_entry<T> *v) { _v = v; }
        vsearch_entry<T>* getEntry(){ return _v; }
        bool operator<(const vsearch<T> n) const {
            return *_v < *(n._v);
        }

    private:
        vsearch_entry<T> *_v;
};

// No native templated typedefs.
#ifndef SearchQueue 
#define SearchQueue std::priority_queue<vsearch<VState*> >
#endif

class HMM {
    public:
        HMM();
        ~HMM();
        HMM(const char*);
        HMM(char*);
        char* generate(int);
        void viterbi(char*, char *qual =NULL);
    private:
        void setTransitions();
        VState* _startState;
        std::vector< VState* > _states;
        MTRand _rng;
};

template <class T>
class Behavior {
    public:
        Behavior(){};
        Behavior(const Behavior<T>&); //intentionally undefined to avoid object-slicing.
        ~Behavior(){};
        virtual T emit(double, int = 0){ return (T) -1; }
        virtual logdouble loglikelihood(T, int=INT_MIN){ 
          logdouble r; 
          r.v = 1.0;
          return r;
        }
        static logdouble loglikelihood(bool, double=DBL_EPSILON, int=INT_MIN);
        virtual void relabelTransition(std::vector<T>&){ return; };
        virtual void enqueueBehavior(std::priority_queue<vsearch<T> >&, vsearch<T>, bool=true, bool=false, bool=false) = 0;
};

// MonoBehavior
//   MonoBehaviors output a single emission.
//   Deterministic emission behavior.
template <class T>
class MonoBehavior : public Behavior<T> {
    public:
        MonoBehavior(T, double);
        ~MonoBehavior();
        virtual T emit(double, int = 0);
        virtual logdouble loglikelihood(T, int=INT_MIN);
        void relabelTransition(std::vector<T>&);
        void enqueueBehavior(std::priority_queue<vsearch<T> >&, vsearch<T>, bool=true,bool=false,bool=false);
    private:
        T _emission;
        double _prob;
        logdouble _likelihood;
        logdouble _notlikelihood;
};

// PolyBehavior
//  
template <class T>
class PolyBehavior : public Behavior<T> {
    public:
        PolyBehavior(TiXmlElement*);
        virtual T emit(double, int = 0);
        virtual logdouble loglikelihood(T, int=INT_MIN);
        void relabelTransition(std::vector<T>&);
        void enqueueBehavior(std::priority_queue<vsearch<T> >&, vsearch<T>, bool=true, bool=false, bool=false);
   private:
        std::map<double, T> _emissions; 
        std::map<T, logdouble> _likelihoods;
        double _density;
};

// IndexedBehavior
template <class T>
class IndexedBehavior : public Behavior<T> {
    public:
        IndexedBehavior(TiXmlElement*);
        IndexedBehavior(std::vector<T>, double = 0.0);
        ~IndexedBehavior();
        virtual T emit(double, int = 0);
        virtual logdouble loglikelihood(T, int, int=0);
        int size(){ return _emissions.size(); }
        void enqueueBehavior(std::priority_queue<vsearch<T> >&, vsearch<T>, bool = false, bool = true, bool = false);
    private:
        std::vector<T> _emissions;
        double _prob;
        logdouble _likelihood;
        logdouble _notlikelihood;
};

class VState {
    friend class HMM;
    public:
        virtual char emit(double, int=0) = 0;
        virtual VState* transition(double, int&) = 0;
        virtual bool hasTransition() = 0;
        virtual bool hasEmission() = 0;
        int getId(){ return _id; };
        std::string getLabel(){ return _label; }
        virtual void enqueueTransitions(SearchQueue&, vsearch_entry<VState*>*) = 0;
        virtual logdouble emissionProbability(char, int = 0, int=INT_MIN) = 0;
        virtual logdouble transitionProbability(VState*,int = 0) = 0;
        virtual bool incrementing() = 0;
        virtual bool resetting() = 0;
    protected:
        virtual void relabelTransition(std::vector< VState* >&) = 0;
        int _id;
        std::string _label;
};

class State : public VState {
    friend class VState;
    public:
        State();
        State(TiXmlElement*);
        State(int, char, int);
        State(std::list<std::pair<double, char> >, std::list<std::pair<double, int> >);
        ~State();
        virtual VState* transition(double, int&);
        virtual bool hasEmission(){ return true; }
        virtual bool hasTransition(){ return true; }
        char emit(double, int=0);
        logdouble emissionProbability(char, int=0, int=INT_MIN);
        logdouble transitionProbability(VState*, int = 0);
        virtual void enqueueTransitions(SearchQueue&, vsearch_entry<VState*>* );
        bool incrementing(){ return _positionIncrement; }
        bool resetting(){ return _positionReset; }
    protected:
        Behavior<VState*> *_transition;
        Behavior<char> *_emission;
        bool _positionReset;
        bool _positionIncrement;
        void relabelTransition(std::vector< VState* >&);
};

class IndexedState : public VState {
    friend class VState;
    public:
        IndexedState();
        IndexedState(TiXmlElement*);
        ~IndexedState();
        VState* transition(double, int&);
        virtual bool hasEmission(){ return true; }
        virtual bool hasTransition(){ return true; }
        char emit(double, int);
        virtual void enqueueTransitions(SearchQueue&, vsearch_entry<VState*>*);
        logdouble emissionProbability(char, int=0, int=INT_MIN);
        logdouble transitionProbability(VState*, int = 0);
        bool incrementing(){ return true; }
        bool resetting(){ return false; }

    private:
        IndexedBehavior<char> *_emissions;
        Behavior<VState*> *_internalTransition;
        Behavior<VState*> *_terminalTransition;
        void relabelTransition(std::vector< VState* >&);
        bool resetPosition(int pos){ return pos >= _emissions->size(); }
};

class SilentState : public State {
    public:
        SilentState(){};
        SilentState(TiXmlElement*);
        SilentState(int, int);
        SilentState(std::list< std::pair<double, int> >);
        bool hasEmission(){ return false; }
};

class AcceptingState : public SilentState {
    public:
        AcceptingState(int id){ _id = id; };
        AcceptingState(TiXmlElement*);
        ~AcceptingState(){};
        bool hasTransition(){ return false; }
        void enqueueTransitions(SearchQueue&, vsearch_entry<VState*>* ){ return; }
};
