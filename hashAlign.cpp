#include "hashAlign.h"
#include <iostream>

using namespace boost;
using namespace std;

// Insert a sequence into the multiple sequence alignment
int MultipleSequenceAlgn::insert(kseq_t *kseq){
  Sequence seq(kseq);
  _sequences.push_back(seq);
  _gaps.clear(); // reset the gaps;
  return _sequences.size();
}

bool MultipleSequenceAlgn::remove(string seq){

  vector<Sequence>::iterator s_itr;

  for(s_itr = _sequences.begin(); s_itr != _sequences.end(); s_itr++){
    if( s_itr->getName() == seq ){
      _sequences.erase(s_itr);
      _gaps.clear();
      return true;
    }
  }

  return false;
}

dynamic_bitset<> MultipleSequenceAlgn::getGaps(){
  if( _gaps.size() > 0){
    return _gaps;
  }

  dynamic_bitset<> mask;
  if( _sequences.size() ){
    mask = _sequences[0].getMask();

    vector<Sequence>::iterator s_itr = _sequences.begin() + 1;
    for( ; s_itr != _sequences.end(); s_itr++ ){
      mask &= s_itr->getMask();
    }
  }

  _gaps = mask;
  return mask;

}

list<Hash> MultipleSequenceAlgn::makeHashes(vector<int> encoding){

  unsigned int hash_size = 24;


  list<Hash> hashes;

  dynamic_bitset<> sequence = _sequences[0].encode(encoding);
  dynamic_bitset<> gaps = getGaps();
  dynamic_bitset<> mask = consensusWithEncoding(encoding);
  unsigned int len = gaps.size();

  cout << mask << endl;

  assert(gaps.size() == sequence.size());
  assert(gaps.size() == mask.size());

  dynamic_bitset<> gapFrame(hash_size);
  dynamic_bitset<> maskFrame(hash_size);
  dynamic_bitset<> seqFrame(hash_size);
  if( gaps.size() < hash_size ){
    //TODO throw an exception?
    return hashes;
  }

  for(unsigned int ii = 0; ii < len; ii++){
    gapFrame     <<= 1;
    gapFrame[0]  =   gaps[len - 1 - ii];
    maskFrame    <<= 1;
    maskFrame[0] =   mask[len - 1 - ii];
    seqFrame     <<= 1;
    seqFrame[0]  =   sequence[len - 1 - ii];

    // Ensure that we've buffered, that we aren't spanning a gap
    // and that the high order bit is set
    if( ii >= hash_size && (~gapFrame).none() && maskFrame[hash_size - 1] ){
      if( maskFrame.count() >= 6 ){
        Hash h(maskFrame, maskFrame & seqFrame, encoding);
        cout << ii - hash_size << ":\t" << maskFrame.count() << "\t" << maskFrame << endl;
        hashes.push_back(h);
      }
    }
  }

  return hashes;
}

vector<int> MultipleSequenceAlgn::getMinimalEncoding(){
  int baseEncoding[] = {0,1,2,3};
  vector<int> encoding(baseEncoding, baseEncoding + 4);

  // How many concordant bits do we have?
  vector< pair< int, vector<int> > > bitCounts;

  // Under each possible encoding, how many concordant bits do we have?
  do{
    dynamic_bitset<> mask = consensusWithEncoding(encoding);
    int concordance = mask.count();
    bitCounts.push_back( pair< int, vector<int> >(concordance, encoding));
  } while( next_permutation(encoding.begin(), encoding.end() ) );

  sort(bitCounts.begin(), bitCounts.end());

  //fprintf(stderr, "%d bits of entropy\n", bitCounts.rbegin()->first);
  // Return the optimal encoding
  return (bitCounts.rbegin())->second;
}

// Find the consensus bits given a particular encoding
dynamic_bitset<> MultipleSequenceAlgn::consensusWithEncoding(vector<int> encoding){

  dynamic_bitset<> consensus;
  dynamic_bitset<> mask;
  if( _sequences.size() ){
    consensus = _sequences[0].encode(encoding);
    mask      = _sequences[0].getMask();

    vector<Sequence>::iterator seq_itr;
    seq_itr = _sequences.begin() + 1;

    // we'll begin by turning mask bits to 0
    for(; seq_itr != _sequences.end(); seq_itr++ ){
      dynamic_bitset<> s = seq_itr->encode(encoding);
      mask &= ~(consensus ^ s); // rm discordant bits from the mask
      mask &= seq_itr->getMask();
      consensus &= s;
    }
  }

  return mask;
}

ReferenceSet MultipleSequenceAlgn::getReferences(){

  vector<int> encoding = getMinimalEncoding();
  list<Hash> hashes = makeHashes(encoding);
  set<Hash> hset = set<Hash>(hashes.begin(), hashes.end());
  ReferenceSet refs(encoding, hset);
  vector<Sequence>::iterator s_itr;

  for( s_itr = _sequences.begin(); s_itr != _sequences.end(); s_itr++){
    refs.insert(Reference(*s_itr, hashes, encoding));
  }

  return refs;
}

Sequence::Sequence(kseq_t* kseq){
  _gapped = -1;
  int len;
  if((len = kseq->name.l)){
    _name = string(kseq->name.s, len);
  } else {
    _name = "";
  }

  if((len = kseq->comment.l)){
    _comment = string(kseq->comment.s, len);
  } else {
    _comment = "";
  }

  if((len = kseq->seq.l)){
    _seq = string(kseq->seq.s, len);
  } else {
    _seq = "";
  }

  if((len = kseq->qual.l)){
    _qual = string(kseq->qual.s, len); 
  } else {
    _qual = ""; 
  };

  // construct the mask
  string::iterator s_itr;
  for( s_itr = _seq.begin(); s_itr != _seq.end(); s_itr++ ){
    if( isNucleotide(*s_itr) ){
      _mask.push_back( true );  _mask.push_back( true );
    } else {
      _mask.push_back( false ); _mask.push_back( false );
    }
  }
}

Sequence::Sequence(string seq){
    _seq = seq;
    _name = "";
    _comment = "";
    _qual = "";
    _mask = dynamic_bitset<>(seq.size(), true);
}

bool Sequence::operator<(const Sequence rhs) const {
  return _name < rhs._name;
}

bool Sequence::isNucleotide(char n){
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

Sequence Sequence::ungapped(){

  Sequence ugs;
  ugs._name = _name;
  ugs._comment = _comment;
  ugs._mask = dynamic_bitset<>();
  ugs._gapped = -1;
  for( unsigned int ii = 0 ; ii < _mask.size() ; ++ii ){
    if( _mask[ii] ){
      ugs._seq.push_back(_seq[ii]);
      ugs._qual.push_back(_qual[ii]);
    }
  }

  ugs._mask.resize(ugs._seq.size(), true);

  return ugs;

}

dynamic_bitset<> Sequence::getMask(){
  return _mask;
}

// Return a sequence given a base encoding scheme
dynamic_bitset<> Sequence::encode(vector<int> encoding){

  dynamic_bitset<> enc;

  string::iterator s_itr;
  for( s_itr = _seq.begin(); s_itr != _seq.end(); s_itr++ ){
    int bits;
    switch( *s_itr ){
      case 'A': 
      case 'a': bits = encoding[0]; break;
      case 'T': 
      case 't': bits = encoding[1]; break;
      case 'C': 
      case 'c': bits = encoding[2]; break;
      case 'G': 
      case 'g': bits = encoding[3]; break;
                // if we encounter something else, we need to have
                // a mask available to store it.
      default: bits = encoding[0];
    }

    enc.push_back( bits % 2 );
    enc.push_back( (bits / 2) % 2 );
  }

  return enc;
}

bool Sequence::isGapped(){ 
  if( _gapped < 0 ){
    _gapped = (~_mask.any() ? 1 : 0);
  }

  return _gapped;
}


// REFERENCE SET

void ReferenceSet::insert(Reference r){
    _references.insert(r);
    return;
}

void ReferenceSet::match(Sequence s){

    dynamic_bitset<> sbits = s.encode(_encoding); 

    for(unsigned int ii = 0; ii < _hashes.size(); ii++){
    }

}

Reference::Reference(Sequence seq, list<Hash> hashes, vector<int> enc){
  _ref = seq;

  list<Hash>::iterator h_itr;

  for( h_itr = hashes.begin(); h_itr != hashes.end(); h_itr++ ){
    set<int> m = h_itr->matches(_ref);
    if(m.size() == 1){
      _positions.insert(pair<Hash, int>(*h_itr, *(m.begin())));
    } else {
        // printf("Hash matched %d times\n", m.size());
    }
  }
}

bool Reference::operator<(const Reference rhs) const {
  return _ref < rhs._ref;
}


Hash::Hash(dynamic_bitset<> mask, dynamic_bitset<> val, vector<int> encoding){
  _value = val;
  _mask = mask;
  _size = mask.size();
  _encoding = encoding;
  printf("got hash %d\n", encoding.size());
}

bool Hash::operator<(const Hash rhs) const {
  return _value < rhs._value;
}

bool Hash::operator>(const Hash rhs) const {
  return _value > rhs._value;
}

set<int> Hash::matches(dynamic_bitset<> encoded){
  //TODO It might make more sense to iterate over the sequence once
  // trying all the hashes as we go. First get it working, then get
  // it right.

  set<int> res;
  dynamic_bitset<> buffer;

  //fill the buffer
  for(unsigned int ii = 0; ii < _mask.size(); ii++){
    buffer.push_back(encoded[ii]);
  }

  // march over the sequence and find hits
  unsigned int index = 0;
  do{
    if( match(buffer) ){
      res.insert(index);
    }
  }while( index++ < encoded.size() - _mask.size() );

  // and return the set
  return res;
}


set<int> Hash::matches(Sequence seq){

  dynamic_bitset<> encoded;
  if( seq.isGapped() ){
    encoded = seq.ungapped().encode(_encoding);   
  } else {
    encoded = seq.encode(_encoding);
  }

  return matches(encoded);

}

bool Hash::match(dynamic_bitset<> t){
  //cout << (_value) << endl << (t & _mask) << endl;
  return (_value == (t & _mask));
}


int main(int argc, char *argv[])
{

  MultipleSequenceAlgn m;

  vector<kseq_t*> sequences;
  gzFile fp;
  kseq_t *seq;
  int l;
  if (argc == 1) {
    fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
    return 1;
  }
  fp = gzopen(argv[1], "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    sequences.push_back(seq);
    m.insert(seq);
  }

  vector<int> encoding = m.getMinimalEncoding();
 
  printf("A -> %d\n", encoding[0]);
  printf("T -> %d\n", encoding[1]);
  printf("C -> %d\n", encoding[2]);
  printf("G -> %d\n", encoding[3]);

  dynamic_bitset<> consensus = m.consensusWithEncoding(encoding);

  dynamic_bitset<> gaps = m.getGaps();

  printf("Gaps:\n");
  for(int ii = 0; ii < (int) gaps.size(); ii++){
    printf("%d", (gaps[ii] ? 1 : 0));
  }
  printf("\n");

  consensus &= gaps;
  vector<int>::iterator e_itr;
  for(e_itr = encoding.begin(); e_itr != encoding.end(); e_itr++){
    printf("%d\t", *e_itr);
  }
  printf("\n");

  for( int ii = 0; ii < (int) consensus.size(); ii++ ){
    printf("%d", (consensus[ii] ? 1 : 0));
  }

  printf("\n");

  printf("%d\t%d\n", consensus.count(), consensus.size());

  list<dynamic_bitset<> >::iterator itvl_itr;

  list<Hash> hashes = m.makeHashes(encoding);

  printf("Got %d hashes\n", hashes.size());

  ReferenceSet refset = m.getReferences();


  Sequence testSeq("CAGGTCACCTTGAAGGAGTCTGGTCCTGTGCTGGTGAAACCCACAGAGACCCTCACGCTGACCTGCACCGTCTCTGGGTTCTCACTCAGCAATGCTAGAATGGGTGTGAGCTGGATCCGTCAGCCCCCAGGGAAGGCCCTGGAGTGGCTTGCACACATTTTTTCGAATGACGAAAAATCCTACAGCACATCTCTGAAGAGCAGGCTCACCATCTCCAAGGACACCTCCAAAAGCCAGGTGGTCCTTACCATGACCAACATGGACCCTGTGGACACAGCCACATATTACTG");

  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

