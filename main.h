#ifndef __MAIN_H__
#define __MAIN_H__

using namespace std;

#define LEN_OF_CHAIN 1164
#define NUM_OF_ADATOMS 364


struct Init {
	int iterations_of_exp;
	int tay;
	int temperature1;
	int temperature2;
	int time_of_nap;
};


extern unsigned int avg_atoms;
extern unsigned int avg_ot;
extern unsigned int avg_nap;
extern double time1; //время

#endif
