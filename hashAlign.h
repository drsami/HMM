#include <zlib.h>
#include <stdio.h>
#include "kseq.h"

#include <algorithm>
#include <string>
#include <vector>

#include <boost/dynamic_bitset.hpp>

KSEQ_INIT(gzFile, gzread)

class Sequence;

class MultipleSequenceAlgn { 
    public:
      MultipleSequenceAlgn(){};
      ~MultipleSequenceAlgn(){};
      int addSequence(kseq_t* s);
      std::vector<int> getMinimalEncoding();
    private:
      
      //Member functions
      std::vector<Sequence> _sequences;
      boost::dynamic_bitset<> consensusWithEncoding(std::vector<int> encoding);
      
      //Member variables
      int _length;
};

class Sequence {
    public:
      Sequence(kseq_t*);
      ~Sequence(){};
      boost::dynamic_bitset<> encode(std::vector<int> encoding);
      boost::dynamic_bitset<> getMask();
    private:
      std::string _name;
      std::string _comment;
      std::string _seq;
      std::string _qual;
      boost::dynamic_bitset<> _mask;
};
