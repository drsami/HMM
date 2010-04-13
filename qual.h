#include "kseq.h"

#ifndef _QUALITY_HMM_
#define _QUALITY_HMM_

// Log-likeihood error probabilities from phred scores
double PHRED_TABLE[] = {0.5};

//default: assume that we're using Sanger reads

#define LOGQUALITY(X) (PHRED_TABLE[X-33])
#define LOGQUALITY_SANGER(X) (PHRED_TABLE[X-33])
#define LOGQUALITY_ILLUMINA(X) (PHRED_TABLE[X-64])
#define LOGQUALITY_SOLEXA(X) (PHRED_TABLE[X-59])

#endif
