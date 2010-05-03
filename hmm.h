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
#include <utility>
#include <vector>

class HMM;
class VState;
class State;

template <class T>
class vsearch_entry {
    //TODO add getters and setters to vsearch_entry
    public:
        vsearch_entry<T>(){}
        T state;
        vsearch_entry<T>* incoming;
        int position;
        int emission;
        double loglikelihood;
        bool operator<(const vsearch_entry<T> n) const {
            return loglikelihood < n.loglikelihood;
        }
};

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
        virtual double loglikelihood(T, int=INT_MIN){ return 1.0; }
        static double loglikelihood(bool, double=DBL_EPSILON, int=INT_MIN);
        virtual void relabelTransition(std::vector<T>&){ return; };
        virtual void enqueueBehavior(std::priority_queue<vsearch_entry<T>* >&, vsearch_entry<T>*, bool=true, bool=false, bool=false) = 0;
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
        virtual double loglikelihood(T, int=INT_MIN);
        void relabelTransition(std::vector<T>&);
        void enqueueBehavior(std::priority_queue<vsearch_entry<T>* >&, vsearch_entry<T>*, bool=true,bool=false,bool=false);
    private:
        T _emission;
        double _prob;
        double _likelihood;
        double _notlikelihood;
};

// PolyBehavior
//  
template <class T>
class PolyBehavior : public Behavior<T> {
    public:
        PolyBehavior(TiXmlElement*);
        virtual T emit(double, int = 0);
        virtual double loglikelihood(T, int=INT_MIN);
        void relabelTransition(std::vector<T>&);
        void enqueueBehavior(std::priority_queue<vsearch_entry<T>* >&, vsearch_entry<T>*, bool=true, bool=false, bool=false);
   private:
        std::map<double, T> _emissions; 
        std::map<T, double> _likelihoods;
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
        virtual double loglikelihood(T, int, int=0);
        int size(){ return _emissions.size(); }
        void enqueueBehavior(std::priority_queue<vsearch_entry<T>* >&, vsearch_entry<T>*, bool = false, bool = true, bool = false);
    private:
        std::vector<T> _emissions;
        double _prob;
        double _likelihood;
        double _notlikelihood;
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
        virtual void enqueueTransitions(std::priority_queue<vsearch_entry<VState*>* >&, vsearch_entry<VState*>*) = 0;
        virtual double emissionProbability(char, int = 0, int=INT_MIN) = 0;
        virtual double transitionProbability(VState*,int = 0) = 0;
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
        double emissionProbability(char, int=0, int=INT_MIN);
        double transitionProbability(VState*, int = 0);
        virtual void enqueueTransitions(std::priority_queue<vsearch_entry<VState*>* >&, vsearch_entry<VState*>* );
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
        virtual void enqueueTransitions(std::priority_queue<vsearch_entry<VState*>* >&, vsearch_entry<VState*>*);
        double emissionProbability(char, int=0, int=INT_MIN){ return DBL_EPSILON; }
        double transitionProbability(VState*, int = 0) { return DBL_EPSILON; }
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
        void enqueueTransitions(std::priority_queue<vsearch_entry<VState*>* >&, vsearch_entry<VState*>* ){ return; }
};
