#ifndef __FILE_OUTPUT_H__
#define __FILE_OUTPUT_H__

using namespace std;

extern int chain_lenth_distribution[LEN_OF_CHAIN];
extern Init init_values;

void file_distribution_output(ofstream& f2out);
void file_chain_output(bool * chain_bool, ofstream& f2out); 



#endif
