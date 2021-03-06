#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include "main.h"
#include "file_output.h"

int chain_lenth_distribution[LEN_OF_CHAIN];

void file_chain_output(bool * chain_bool, ofstream& f2out)
{
	int k = 0;
	int cnt = 0;
	for (k; chain_bool[k]; k++) {}
	for (int n = (k + 1) % LEN_OF_CHAIN; n != k; n = (n + 1) % LEN_OF_CHAIN) {
		if (chain_bool[n]) {
			cnt++;
		}
		else {
			if (cnt != 0) {
				f2out << cnt << " | ";
				chain_lenth_distribution[cnt]++;
				cnt = 0;
			}
		}
	}
	if (cnt != 0) {
		f2out << cnt << endl;
		chain_lenth_distribution[cnt]++;
		cnt = 0;
	}
	else {
		f2out << endl;
	}
}

void file_distribution_output(ofstream& f2out) 
{
	f2out << "Среднее число напыленных атомов: " << avg_atoms / (init_values.iterations_of_exp) << endl;
	f2out << endl << "(Длина, число отсчетов)" << endl;
	for (int i = 0; i < NUM_OF_ADATOMS + 1; i++) {
		f2out << i << " " << chain_lenth_distribution[i] << endl;
	}
	for (int i = NUM_OF_ADATOMS + 1; i < LEN_OF_CHAIN; i++) {
		if (chain_lenth_distribution[i] != 0) f2out << i << " " << chain_lenth_distribution[i] << endl;
	}
}  

void file_output_description(ofstream& f2out)
{

	f2out << -1 << endl;
    std::cout << "Время выполнения программы " << search_time << " мс";
	f2out << "Время выполнения программы " << search_time << " мс" << endl;
    f2out << -1 << endl;
	f2out << "T1=" << init_values.temperature1 << "K ";
	f2out << "T2=" << init_values.temperature2 << "K ";
	f2out << "Exp=" << init_values.iterations_of_exp << endl;
	f2out << "Адатомов: " << avg_atoms / (init_values.iterations_of_exp) << "/" << LEN_OF_CHAIN;
	f2out << " Тау=" << init_values.tay << " секунд" << endl;

	f2out << "Ср. вр. нап: " << avg1_time / (init_values.iterations_of_exp) << "сек" << endl;
	f2out << "Ср. ч-ло ит. нап: " << avg_nap / (init_values.iterations_of_exp) << endl;
	f2out << "Ср. т-ра после нап: " << avg1_temp / (init_values.iterations_of_exp) << "K" << endl;

	f2out << "Ср. вр. отжига: " << avg2_time / (init_values.iterations_of_exp) << "cек" << endl;
	f2out << "Ср. ч-ло ит. отжига: " << avg_ot / (init_values.iterations_of_exp) << endl;
	f2out << "Ср. т-ра после отжига: " << avg2_temp / (init_values.iterations_of_exp) << endl;
}
