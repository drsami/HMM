#include "hasher.h"

using namespace boost;
using namespace std;

Hash::Hash(dynamic_bitset<> mask, dynamic_bitset<> val){
  _id = HID++;
  _value = val & mask;
  _mask = mask;
  _size = mask.size();
}

bool Hash::operator<(const Hash rhs) const {
  return (_value & _mask) < (rhs._value & rhs._mask);
}

bool Hash::operator>(const Hash rhs) const {
  return (_value & _mask) > (rhs._value & rhs._mask);
}

bool Hash::match(dynamic_bitset<> t){
  return (_value == (t & _mask));
}

bool Hash::match(dynamic_bitset<> val, dynamic_bitset<> valmask){
  if( _mask.is_subset_of(valmask) ){
    return (_value == (val & _mask));
  } else {
    return false;
  }
}

// HashSet
HashSet::HashSet(){

}

bool HashSet::insert(Hash h){

    if( _hashes.size() ){
        assert(h.size() == _size );
    } else {
        _size = h.size();
    }
    
    return (_hashes.insert(h)).second;
}

void HashSet::match(dynamic_bitset<> val, dynamic_bitset<> mask, multimap<int,int> &matches, int pos){
    set<Hash>::iterator h_itr;

    for(h_itr = _hashes.begin(); h_itr != _hashes.end(); h_itr++){
        Hash h = *h_itr;
        if(h.match(val, mask)){
            int ids = h.id();
            matches.insert(pair<int, int>(pos, ids));
        }
    }
    return;
}

multimap<int, int> HashSet::matchAll(Sequence seq){

    multimap<int, int> res;
    dynamic_bitset<> eseq = seq.getEncoded();
    unsigned int seqlen = eseq.size();

    dynamic_bitset<> val_set;
    dynamic_bitset<> mask_set;

    // Buffer
    for( unsigned int ii = 0; ii < _size; ii++ ){
      val_set.push_back(eseq[_size-ii-1]);
      mask_set.push_back(1);
    }

    for( unsigned int ii = 0; ii < seqlen - _size; ii++ ){
        match(val_set, mask_set, res, ii);
        val_set >>= 1;
        val_set[_size-1] = eseq[_size+ii];
    }

    return res;
}

Sequence::Sequence(kseq_t* kseq){
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

  _gapMask = dynamic_bitset<>(len*2);

  // construct the mask
  string::iterator s_itr;
  for( s_itr = _seq.begin(); s_itr != _seq.end(); s_itr++ ){
    _gapMask <<= 2;
    if( isNucleotide(*s_itr) ){
      _gapMask |= dynamic_bitset<>(2,3);
    }
  }
}

// MultipleSequenceAlign

void MultipleSequenceAlign::reset(){
    _consensus.clear();
    _gaps.clear();
    _mask.clear();
}

// Insert a sequence into the multiple sequence alignment
int MultipleSequenceAlign::insert(kseq_t *kseq){
  Sequence seq(kseq);
  _sequences.push_back(seq);
  reset();
  return _sequences.size();
}

bool MultipleSequenceAlign::remove(string seq){
  vector<Sequence>::iterator s_itr;
  for(s_itr = _sequences.begin(); s_itr != _sequences.end(); s_itr++){
    if( s_itr->getName() == seq ){
      _sequences.erase(s_itr);
      reset();
      return true;
    }
  }
  return false;
}

dynamic_bitset<> MultipleSequenceAlign::getGaps(){
  if( _gaps.size() > 0){
    return _gaps;
  }

  dynamic_bitset<> mask(0);
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

dynamic_bitset<> MultipleSequenceAlign::getConsensus(){
    if( _consensus.size() > 0 ){
        return _consensus;
    }

    dynamic_bitset<> consensus(0);
  if( _sequences.size() ){
    consensus = _sequences[0].getEncoded();

    vector<Sequence>::iterator s_itr = _sequences.begin() + 1;
    for( ; s_itr != _sequences.end(); s_itr++ ){
      consensus &= s_itr->getEncoded();
    }
  }

  _consensus = consensus;
  return consensus;
}

dynamic_bitset<> MultipleSequenceAlign::getMask(){
    if( _mask.size() > 0 ){
        return _mask;
    }

    dynamic_bitset<> prev = _sequences[0].getEncoded();
    dynamic_bitset<> accumulator = dynamic_bitset<>(prev.size());
    accumulator.flip();


    vector<Sequence>::iterator s_itr = _sequences.begin() + 1;
    for( ; s_itr != _sequences.end(); s_itr++ ){
        dynamic_bitset<> curr = s_itr->getEncoded();
        accumulator -= (curr ^ prev);
        prev = curr;
    }

    dynamic_bitset<> gaps = getGaps();
    accumulator &= gaps;
    _mask = accumulator;
    return accumulator;
}

HashSet MultipleSequenceAlign::makeHashes(){
  
  HashSet hs;
  multiset<Hash> h;
  
  int hash_length = 24;

  dynamic_bitset<> consensus = getConsensus();
  dynamic_bitset<> gaps      = getGaps();
  dynamic_bitset<> mask      = getMask();

  int len = mask.size();

  dynamic_bitset<> hashMask(hash_length);     // tell us /where/ the invariant bits are in the hash
  dynamic_bitset<> hashSequence(hash_length); // tell us /what/ the invariant bits are in the hash

  dynamic_bitset<> gapWindow(hash_length); // Is there a gap here?
  // buffer
  for( int ii = 0; ii < hash_length; ii++ ){
    hashMask[ii]     = mask[ii];
    hashSequence[ii] = consensus[ii];
    gapWindow[ii]    = gaps[ii];
  }

  for( int ii = hash_length; ii < len; ii++ ){
      // The least significant bit is unmasked and we aren't in a gap
      if( hashMask[0] && (~gapWindow).none() ){
         Hash hash(hashMask, hashSequence);
         h.insert(hash);
      }

      // right shift everything dumping the least significant bit
      hashMask     >>= 1;
      hashSequence >>= 1;
      gapWindow    >>= 1;
      
      // reset the most significant bit
      hashMask[hash_length - 1] = mask[ii];
      hashSequence[hash_length - 1] = consensus[ii];
      gapWindow[hash_length - 1] = gaps[ii];
  }

  multiset<Hash>::iterator h_itr;

  // Remove duplicates
  for( h_itr = h.begin() ; h_itr != h.end(); h_itr++ ){
    if(h.count(*h_itr) == 1){
        hs.insert(*h_itr);
    }
  }

  return hs;
}

int main(int argc, char *argv[])
{

  dynamic_bitset<> hash(5);
  hash[0] = 1;
  hash[1] = 1;
  hash[4] = 1;
  dynamic_bitset<> mask(5);
  mask[0] = 1;
  mask[1] = 1;
  mask[4] = 1;
  dynamic_bitset<> val(5);
  val[0] = 1;
  val[1] = 1;
  val[2] = 1;
  val[3] = 0;
  val[4] = 1;
  Hash h(mask, hash);

  return 0;
}

