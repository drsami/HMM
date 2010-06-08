#include "hasher.h"

using namespace boost;
using namespace std;

Hash::Hash(dynamic_bitset<> mask, dynamic_bitset<> val){
  _id = HID++;
  _value = val;
  _mask = mask;
  _size = mask.size();
  assert(_mask.size() == _value.size());
}

bool Hash::operator<(const Hash rhs) const {
  return _value < rhs._value;
}

bool Hash::operator>(const Hash rhs) const {
  return _value > rhs._value;
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

  // construct the mask
/*  string::iterator s_itr;
  for( s_itr = _seq.begin(); s_itr != _seq.end(); s_itr++ ){
    if( isNucleotide(*s_itr) ){
      _mask.push_back( true );  _mask.push_back( true );
    } else {
      _mask.push_back( false ); _mask.push_back( false );
    }
  }
  */
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

