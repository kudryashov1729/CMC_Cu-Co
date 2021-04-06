#ifndef __CMC_H__
#define __CMC_H__

struct Rate_Catalog
{
	double move_left_barieer;   
	double move_right_barieer;   
	double rate_left;
	double rate_right;
};
struct Possible_Events 
{
	int from;
	int to;
	double rate;
}; 

extern Rate_Catalog * mas_rate;
extern Possible_Events * events;


void find_rang(bool * chain_bool, int pos, double t);
void find_rate(bool * chain_bool, int pos, double temperature);
void create_rate_catalog(bool * chain_bool, double t);
void change_rate_catalog(bool * chain_bool, double t);
void choose_event(int * chain_int, bool* chain_bool, bool new_atoms);

extern int mem_pos_to; 

#endif
