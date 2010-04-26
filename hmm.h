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

class HMM {
    public:
        HMM();
        ~HMM();
        HMM(const char*);
        HMM(char*);
        char* generate(int);
        std::list< int > viterbi(char*, char*);
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
        ~Behavior(){};
        virtual T emit(double, int = 0){ return (T) -1; }
        virtual double loglikelihood(T, int=INT_MIN){ return 1.0; }
        static double loglikelihood(bool, double=DBL_EPSILON, int=INT_MIN);
        virtual void relabelTransition(std::vector<T>&){ return; };
};

// MonoBehavior
//   MonoBehaviors output a single emission.
//   Deterministic emission behavior.
template <class T>
class MonoBehavior : public Behavior<T> {
    public:
        MonoBehavior(T, double = 0.0);
        ~MonoBehavior();
        virtual T emit(double, int = 0);
        virtual double loglikelihood(T, int=INT_MIN);
        void relabelTransition(std::vector<T>&);
    private:
        T _emission;
        double _mutr;
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
    private:
        std::vector<T> _emissions;
        double _mutr;
};

class VState {
    friend class HMM;
    public:
        virtual char emit(double, int=0) =0;
        virtual VState* transition(double, int&) = 0;
        virtual bool hasTransition() = 0;
        virtual bool hasEmission() = 0;
        int getId(){ return _id; };
        std::string getLabel(){ return _label; }
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
        double emissionProbability(char, int);
        double transitionProbability(int);
    protected:
        Behavior<VState*> *_transition;
        Behavior<char> *_emission;
        bool _positionReset;
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
};
