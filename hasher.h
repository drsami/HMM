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
    boost::dynamic_bitset<> getMask(){ return _qualityMask & _gapMask; }
  private:
    std::string _name;
    std::string _comment;
    std::string _seq;
    std::string _qual;
    boost::dynamic_bitset<> _seqbits;
    boost::dynamic_bitset<> _qualityMask;
    boost::dynamic_bitset<> _gapMask;
};
