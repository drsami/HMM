#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <zlib.h>
#include "kseq.h"

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

KSEQ_INIT(gzFile, gzread);

enum nucleotide {T,G,C,A};

class Hash;
class HashSet;
class Sequence;

unsigned int HID = 0;
class Hash {
    public:
        Hash();
        Hash(boost::dynamic_bitset<>, boost::dynamic_bitset<>);
        ~Hash(){};
        bool match(boost::dynamic_bitset<>);
        bool match(boost::dynamic_bitset<>, boost::dynamic_bitset<>);
        bool operator<(const Hash) const;
        bool operator>(const Hash) const;
        unsigned int id(){ return _id; }
        unsigned int size(){ return _size; }
        boost::dynamic_bitset<> getMask(){ return _mask; }
    private:
        unsigned int _id;
        unsigned int _size;
        boost::dynamic_bitset<> _mask;
        boost::dynamic_bitset<> _value;
};

class HashSet {
    public:
        HashSet();
        bool insert(Hash);
        void match(boost::dynamic_bitset<>, boost::dynamic_bitset<>, std::multimap<int, int>&, int);
        std::multimap<int, int> matchAll(Sequence);
        int size(){ return _size; }
    private:
        std::set<Hash> _hashes;
        unsigned int _size;
};

class Sequence {
    public:
        Sequence(kseq_t*);
        Sequence ungap();
        boost::dynamic_bitset<> getEncoded(){ return _seqbits; }
        boost::dynamic_bitset<> getMask(){ return _gapMask; }
        std::string getName(){ return _name; }
    private:

        inline bool isNucleotide(char n){
            switch(n){
                case 'A':
                case 'a':
                case 'T':
                case 't':
                case 'C':
                case 'c':
                case 'G':
                case 'g': return true;
                default: return false;
            }
        }

        std::string _name;
        std::string _comment;
        std::string _seq;
        std::string _qual;
        boost::dynamic_bitset<> _seqbits;
        boost::dynamic_bitset<> _qualityMask;
        boost::dynamic_bitset<> _gapMask;
};
/*
class Reference : public Sequence {
    public:
        Reference();
    private:
        std::map<int, int> _hashMap;
};
*/

class MultipleSequenceAlign { 
    public:
        MultipleSequenceAlign(){};
        ~MultipleSequenceAlign(){};
        int insert(kseq_t* s);
        bool remove(std::string);
        std::vector<int> getMinimalEncoding();
        boost::dynamic_bitset<> getMask();
        boost::dynamic_bitset<> getGaps();
        boost::dynamic_bitset<> getConsensus();
        HashSet makeHashes();
    private:
        //Member functions
        void reset();
        //Member variables
        std::vector<Sequence> _sequences;
        boost::dynamic_bitset<> _consensus;
        boost::dynamic_bitset<> _gaps;
        boost::dynamic_bitset<> _mask;
        int _length;
};


