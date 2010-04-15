#include "kseq.h"
#include "qual.h"

#include <float.h>
#include <math.h>

#include <list>
#include <map>
#include <queue>
#include <utility>
#include <vector>

class HMM;
class State;
class Transitions;
class Emissions;

typedef struct {

    int predecessor;
    int emission;
    double best_score;

} vcell;

typedef struct {
    int out_state;
    int out_emit;
    int in_state;
    int in_emit;
    double cost;
} edge;

class HMM {
    public:
        HMM();
        ~HMM();
        HMM(char*);
        char* generate(int);
        std::list< int > viterbi(char*, char*);
    private:
        int _startState;
        State* getState(int x);
        std::vector< State* > _states;

};

template <class T>
class Behavior {
    public:
        Behavior(){};
        ~Behavior(){};
        virtual T emit(double){ return (T) -1; }
        virtual void enqueueEmissions(T, int, int, std::queue< edge >&){ return; }
        virtual double loglikelihood(T, int){ return 1.0; }
};

// MonoBehavior
//   MonoBehaviors output a single emission.
//   Deterministic emission behavior.
template <class T>
class MonoBehavior : public Behavior<T> {
    public:
        MonoBehavior(T);
        ~MonoBehavior();
        virtual T emit(double);
        virtual void enqueueEmissions(T, int, int, std::queue< edge >&);
        virtual double loglikelihood(T, int);
    private:
        T _emission;
};

// PolyBehavior
//  
template <class T>
class PolyBehavior : public Behavior<T> {
    public:
        PolyBehavior(std::list<std::pair<double, T> >);
        virtual T emit(double);
        virtual void enqueueEmissions(T, int, int, std::queue< edge >&);
        virtual double loglikelihood(T, int);
   private:
        std::map<double, T> _emissions;    
        std::list< std::pair<double, T> > _likelihoods;
};

class State {
    public:
        State();
        State(int, char, int);
        State(std::list<std::pair<double, char> >, std::list<std::pair<double, int> >);
        ~State();
        int getId(){ return _id; };
        int transition(double);
        virtual bool hasEmission(){ return true; }
        virtual bool hasTransition(){ return true; }
        char emit(double);
        double emissionProbability(char, int);
        double transitionProbability(int);
        void enqueueTransitions(int, std::queue< edge >&);
    protected:
        int _id;
        Behavior<int> *_transition;
        Behavior<char> *_emission;
};

class SilentState : public State {
    public:
        SilentState(){};
        SilentState(int, int);
        SilentState(std::list< std::pair<double, int> >);
        bool hasEmission(){ return false; }
};

class AcceptingState : public SilentState {
    public:
        AcceptingState(int id){ _id = id; };
        ~AcceptingState(){};
        bool hasTransition(){ return false; }
};
