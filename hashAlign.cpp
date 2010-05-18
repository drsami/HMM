#include "hashAlign.h"

using namespace boost;
using namespace std;

// Insert a sequence into the multiple sequence alignment
int MultipleSequenceAlgn::addSequence(kseq_t *kseq){
    Sequence seq(kseq);
    _sequences.push_back(seq);
    return _sequences.size();
}

dynamic_bitset<> MultipleSequenceAlgn::getGaps(){
    dynamic_bitset<> mask;
    if( _sequences.size() ){
        mask = _sequences[0].getMask();

        vector<Sequence>::iterator s_itr = _sequences.begin() + 1;
        for( ; s_itr != _sequences.end(); s_itr++ ){
            mask &= s_itr->getMask();
        }
    }

    for(int ii = 0; ii < (int) mask.size(); ii++){
        printf("%d", (mask[ii] ? 1 : 0));
    }
    printf("\n");
    return mask;

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
  string::iterator s_itr;
  for( s_itr = _seq.begin(); s_itr != _seq.end(); s_itr++ ){
    switch( *s_itr ){
      case 'A':
      case 'a':
      case 'T':
      case 't':
      case 'C':
      case 'c':
      case 'G':
      case 'g': _mask.push_back( true );  _mask.push_back( true ); break;
      default:  _mask.push_back( false ); _mask.push_back( false );
    }
  }
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
/*
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
    m.addSequence(seq);
    printf("return value: %d\n", l);
  }

  vector<int> encoding = m.getMinimalEncoding();
  dynamic_bitset<> consensus = m.consensusWithEncoding(encoding);
  consensus &= m.getGaps();
  vector<int>::iterator e_itr;
  for(e_itr = encoding.begin(); e_itr != encoding.end(); e_itr++){
    printf("%d\t", *e_itr);
  }
  printf("\n");

  for( int ii = 0; ii < (int) consensus.size(); ii++ ){
    printf("%d", (consensus[ii] ? 0 : 1));
  }

  printf("\n");

  printf("%d\t%d\n", consensus.count(), consensus.size());

  //m.getGaps();

  vector<kseq_t*>::iterator s_itr;

  for( s_itr = sequences.begin(); s_itr != sequences.end(); s_itr++ ){
    seq = *s_itr;
    printf("name: %s\n", seq->name.s);
    if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    printf("seq: %s\n", seq->seq.s);
    if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
  }

  kseq_destroy(seq);
  gzclose(fp);
  return 0;
}

*/
