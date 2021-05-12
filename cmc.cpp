#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include "main.h"
#include "cmc.h"

#define K0_const 1e12 //секунда^(-1)
#define k_bol 8.625E-5 //константа Больцмана еВ/К

int mem_pos_to = 0; 


Rate_Catalog * mas_rate = new Rate_Catalog[LEN_OF_CHAIN];
Possible_Events * events = new Possible_Events[NUM_OF_ADATOMS * 2 + 10];


/*
    // инициализация двумерного массива:
    int a[5][3] = { {4, 7, 8}, {9, 66, -1}, {5, -5, 0}, {3, -3, 30}, {1, 1, 1} };
    
    a[лев - 1][прав - 1]
    
    a[4][4] =   
    a[0][0]     a[0][1]     a[0][2]     a[0][3]
    a[1][0]     a[1][1]     a[1][2]     a[1][3]    
    a[2][0]     a[2][1]     a[2][2]     a[2][3]
    a[3][0]     a[3][1]     a[3][2]     a[3][3]
    
 */

//все прыжки вправо
double energy_barrier[5][5] ={
    {   1000,   0.7,    0.9,    0.8,    0.8},
    {   1000,   0.2,    0.3,    0.2,    0.2},
    {   1000,   0.2,    0.2,    0.3,    0.3},
    {   1000,   0.4,    0.2,    0.3,    0.3},
    {   1000,   0.4,    0.2,    0.3,    0.3}
};

/*
 
O_____O_____O__ Left:   0.3     Ridht:  0.3
O_____O____O___ Left:   0.3     Ridht:  0.3
O_____O___O____ Left:   0.3     Ridht:  0.3
O_____O__O_____ Left:   0.3     Ridht:  0.2
O_____O_O______ Left:   0.2     Ridht:  0.4
O_____OO_______ Left:   0.8     Ridht:  1000
_O____O_____O__ Left:   0.3     Ridht:  0.3
_O____O____O___ Left:   0.3     Ridht:  0.3
_O____O___O____ Left:   0.3     Ridht:  0.3
_O____O__O_____ Left:   0.3     Ridht:  0.2
_O____O_O______ Left:   0.2     Ridht:  0.4
_O____OO_______ Left:   0.8     Ridht:  1000
__O___O_____O__ Left:   0.3     Ridht:  0.3
__O___O____O___ Left:   0.3     Ridht:  0.3
__O___O___O____ Left:   0.3     Ridht:  0.3
__O___O__O_____ Left:   0.3     Ridht:  0.2
__O___O_O______ Left:   0.2     Ridht:  0.4
__O___OO_______ Left:   0.8     Ridht:  1000
___O__O_____O__ Left:   0.2     Ridht:  0.3
___O__O____O___ Left:   0.2     Ridht:  0.3
___O__O___O____ Left:   0.2     Ridht:  0.3
___O__O__O_____ Left:   0.2     Ridht:  0.2
___O__O_O______ Left:   0.3     Ridht:  0.2
___O__OO_______ Left:   0.9     Ridht:  1000
____O_O_____O__ Left:   0.4     Ridht:  0.2
____O_O____O___ Left:   0.4     Ridht:  0.2
____O_O___O____ Left:   0.4     Ridht:  0.2
____O_O__O_____ Left:   0.2     Ridht:  0.3
____O_O_O______ Left:   0.2     Ridht:  0.2
____O_OO_______ Left:   0.7     Ridht:  1000
_____OO_____O__ Left:   1000    Ridht:  0.8
_____OO____O___ Left:   1000    Ridht:  0.8
_____OO___O____ Left:   1000    Ridht:  0.8
_____OO__O_____ Left:   1000    Ridht:  0.9
_____OO_O______ Left:   1000    Ridht:  0.7
_____OOO_______ Left:   1000    Ridht:  1000

 */

void find_rang(bool * chain_bool, int pos, double t) 
{ // содержит вызов find_rate()
	int l = 1;
	while (chain_bool[((int)LEN_OF_CHAIN + (pos - l)) % LEN_OF_CHAIN] != true && l < 5) {// ищем разность позициий данного атома и ближайшего слева
		l++;
	}
	int r = 1;
	while (chain_bool[(pos + r) % LEN_OF_CHAIN] != true && r < 5) {// ищем разность позициий данного атома и ближайшего справа
		r++;
	}
	l--;
    r--;
	l = l > 4 ? 4 : l;
    r = r > 4 ? 4 : r;
	mas_rate[pos].move_right_barieer = energy_barrier[l][r];
    
    int m = l;
    l = r;
    r = m;
    mas_rate[pos].move_left_barieer = energy_barrier[l][r];
	
	find_rate(chain_bool, pos, t);
}; 

void find_rate(bool * chain_bool, int pos, double temperature)
{
	if(mas_rate[pos].move_left_barieer > 0)
	mas_rate[pos].rate_left = (double)( K0_const )* exp((double)((-1)* mas_rate[pos].move_left_barieer / (k_bol) / temperature));
	if (mas_rate[pos].move_right_barieer > 0)
	mas_rate[pos].rate_right = (double)( K0_const )* exp((double)((-1)* mas_rate[pos].move_right_barieer / (k_bol) / temperature));
};

void create_rate_catalog(bool * chain_bool, double t) 
{
	for (int i = 0; i < LEN_OF_CHAIN; i++) {
		if (chain_bool[i]) {
			find_rang(chain_bool, i, t);
		}
	}
};

void change_rate_catalog(bool * chain_bool, double t)
{
	find_rang(chain_bool, mem_pos_to, t);
	for (int i = 1; i < 10; i++) {
		if (chain_bool[(LEN_OF_CHAIN + mem_pos_to - i) % LEN_OF_CHAIN]) {
			find_rang(chain_bool, (LEN_OF_CHAIN + mem_pos_to - i) % LEN_OF_CHAIN, t);
			break;
		}
	}
	for (int i = 1; i < 10; i++) {
		if (chain_bool[(mem_pos_to + i) % LEN_OF_CHAIN]) {
			find_rang(chain_bool, (mem_pos_to + i) % LEN_OF_CHAIN, t);
			break;
		}
	}
};
void choose_event(int * chain_int, bool* chain_bool, bool new_atoms)
{
	//
	int k = 1; //индекс свободного места в масиве возможных событий events[]
	int m = 0; //индекс атома в масиве chain_int[], ещё недобавленного в возможное событие 
	while (chain_int[m] >= 0) {
		int pos = (chain_int[m] - 1 + LEN_OF_CHAIN) % LEN_OF_CHAIN;
		if (!chain_bool[pos]) { // если слева от атома есть место
			events[k].index_from = m; //сохряняем номер ячейки массива chain_int[], который хранит позицию откуда прыгает атом
			events[k].to = pos;
			events[k].rate = mas_rate[chain_int[m]].rate_left;
			k++;
		}
		pos = (pos + 2) % LEN_OF_CHAIN;
		if (!chain_bool[pos]) { // если справа от атома есть место
			events[k].index_from = m; //сохряняем номер ячейки массива chain_int[], который хранит позицию откуда прыгает атом
			events[k].to = pos;
			events[k].rate = mas_rate[chain_int[m]].rate_right;
			k++;
		}
		m++;
	}
	double sum[NUM_OF_ADATOMS * 2 + 10]; //массив сумм всех вероятностей


	if (new_atoms) { // если напыление включено
		sum[0] = 2.7;
	}
	else { //напыление выключено
		sum[0] = 0;
	}

	int i = 1;
	while (i < k) {
		sum[i] = sum[i - 1] + events[i].rate;
		i++;
	}
	int r1 = rand();
	double rand1 = ((double)(r1) / RAND_MAX) * sum[i - 1];
	int j = new_atoms ? 0 : 1;
	for (j; sum[j] < rand1; j++) {}
	if (j >= k) {
		j = 1;
		std::cout << "Out of array events" << endl;
	}

	if (j == 0) { //событие напыления
		if (!new_atoms) {
			std::cout << "Try to add atom" << endl;
		}
		int rand2 = 0;
		if (chain_int[NUM_OF_ADATOMS] >= 0) {
			std::cout << "Too many atoms" << endl;
		}
		do {
			rand2 = (int)(((double)(rand()) / RAND_MAX) * (LEN_OF_CHAIN - 1)); //случайная позиция
		} while (chain_bool[rand2]); //до тех пор пока не попадем в пустую позицию
		chain_int[m] = rand2;
		chain_bool[rand2] = true;
		mem_pos_to = rand2;
	}
	else { //событие движения
		chain_bool[chain_int[events[j].index_from]] = false;
        chain_int[events[j].index_from] = events[j].to;
		chain_bool[events[j].to] = true;
		mem_pos_to = events[j].to;
	}
	//Счет времени
	int r2 = rand() + 1;
    double random_max_value = RAND_MAX;
	double l = log((random_max_value + 1) / double(r2));
	time1 = time1 + (1 / sum[i - 1]) * l;
}

