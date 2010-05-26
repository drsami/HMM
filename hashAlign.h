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

class Sequence;
class Hash;
class Reference;

class MultipleSequenceAlgn { 
  public:
    MultipleSequenceAlgn(){};
    ~MultipleSequenceAlgn(){};
    int insert(kseq_t* s);
    bool remove(std::string);
    std::vector<int> getMinimalEncoding();
    boost::dynamic_bitset<> getGaps();
    std::list<Hash> makeHashes(std::vector<int>);
    boost::dynamic_bitset<> consensusWithEncoding(std::vector<int> encoding);
    std::set<Reference> getReferences();
  private:

    //Member functions
    std::vector<Sequence> _sequences;

    //Member variables
    boost::dynamic_bitset<> _gaps;
    int _length;
};

class Hash {
  public:
    Hash(){ _size = 0;};
    Hash(boost::dynamic_bitset<>, boost::dynamic_bitset<>, std::vector<int>);
    ~Hash(){};
    std::set<int> matches(Sequence);
    std::set<int> matches(boost::dynamic_bitset<>);
    bool match(boost::dynamic_bitset<>);
    bool operator<(const Hash) const;
    bool operator>(const Hash) const;
    unsigned int size(){ return _size; }
    std::vector<int> getEncoding() const { return _encoding; }
    std::vector<int> getEncoding() { return _encoding; }
    boost::dynamic_bitset<> getMask(){ return _mask; }
  private:
    std::vector<int> _encoding;
    boost::dynamic_bitset<> _mask;
    boost::dynamic_bitset<> _value;
    int _size;
};

class Sequence {
  public:
    Sequence(kseq_t*);
    Sequence(){};
    ~Sequence(){};
    Sequence ungapped();
    boost::dynamic_bitset<> encode(std::vector<int> encoding);
    boost::dynamic_bitset<> getMask();
    std::string getName(){ return _name; }
    bool isGapped();
    bool operator<(const Sequence) const;
  private:
    bool isNucleotide(char);
    std::string _name;
    std::string _comment;
    std::string _seq;
    std::string _qual;
    boost::dynamic_bitset<> _mask;
    int _gapped;
};

typedef struct {
  Hash hash;
  int expected_position;
  int observed_position;
  int offset;
} hash_info;

bool operator<(const hash_info lhs, const hash_info rhs) {
  return lhs.observed_position < rhs.observed_position;
}


class Reference {
  public:
    Reference(Sequence, std::list<Hash>);
    void match(Sequence);
    bool operator<(const Reference) const;
  private:
    void vote(boost::dynamic_bitset<> &voters, std::map<int, int> &votes, hash_info&);
    Sequence _ref;
    std::map<Hash, int> _positions;
};
