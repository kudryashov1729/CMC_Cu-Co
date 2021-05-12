#include <iostream>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <math.h>

#include "main.h"
#include "file_output.h"
#include "cmc.h"


Init init_values;

unsigned int avg_atoms;
unsigned int avg_ot;
unsigned int avg_nap;

double avg1_time = 0;
double avg2_time = 0;
double avg1_temp = 0;
double avg2_temp = 0;
unsigned int search_time;
/*
string to_string(int n)
{
    char buf[40];
    sprintf(buf,"%d",n);
    return buf;
}*/

double time1; //время
double T1; //температура
string s = "exp="; //имя выводимого файла 

bool init(int argc, char* argv[], ofstream& f2out) {
	if (argc == 6) {
		init_values.iterations_of_exp =   atoi(argv[1]);
		init_values.tay =                 atoi(argv[2]);
		init_values.temperature1 =        atoi(argv[3]);
		init_values.temperature2 =        atoi(argv[4]);
		init_values.time_of_nap =         atoi(argv[5]);
        
        s = "exp=" + to_string(init_values.iterations_of_exp) + "_tay=" + to_string(init_values.tay) + "_T1=" + to_string(init_values.temperature1) + "_T2=" + to_string(init_values.temperature2) + "_time=" + to_string(init_values.time_of_nap);
        
		f2out.open(s, ios::out);
		return false;
	}
	std::cout << endl << "Input:\n iterations_of_exp\n tay\n temperature1\n temperature2\n time_of_nap" << endl;
	return true;
}

int main(int argc, char* argv[])
{
	
    setlocale(LC_ALL, "Russian");
    
    ofstream f2out;
	
    if (init(argc, argv, f2out)) {
		delete[] mas_rate;
		delete[] events;
		return 0;
	}

	unsigned int start_time = clock();

	int * chain_int = new int[LEN_OF_CHAIN];
	bool * chain_bool = new bool[LEN_OF_CHAIN];

	avg_atoms = 0;
	avg_ot = 0;
	avg_nap = 0;
    
    for (int i = 0; i < LEN_OF_CHAIN; i++) {
		chain_lenth_distribution[i] = 0;
	}
	

	for (int index = 0; index < init_values.iterations_of_exp; index++) {
		if (int(init_values.iterations_of_exp / 100) != 0) {
			if (int(index % int(init_values.iterations_of_exp / 100)) == 0) std::cout << "\r" << int(100 * double(index) / init_values.iterations_of_exp) << "%";
			if ((index + 1) % init_values.iterations_of_exp == 0) std::cout << "\r100%" << endl;
		}

		T1 = init_values.temperature1;
		time1 = 0;

		f2out << "\\\\\\\\\\\\\\ \t Серия " << index + 1 << "\t \\\\\\\\\\\\\\" << endl;
		srand(index);
		for (int i = 0; i < LEN_OF_CHAIN; i++) {
			chain_int[i] = -1;
			chain_bool[i] = false;
		}
		
		//1 ЭТАП (Напыление)
		int i_nap = 1;

		for (i_nap = 1; chain_int[NUM_OF_ADATOMS - 1] < 0; i_nap++) {
			choose_event(chain_int, chain_bool, true);
			change_rate_catalog(chain_bool, T1);
		}

		//вывод в файл
		f2out << "НАПЫЛЕНИЕ" << endl;
		f2out << "Температура после напыления: " << T1 << endl;
		int count_number_of_adatoms = 0;
		for (int i = 0; i < LEN_OF_CHAIN; i++) { if (chain_bool[i]) count_number_of_adatoms++; }
		avg_atoms = avg_atoms + count_number_of_adatoms;
		f2out << "Время напыления: " << time1 << endl;
		f2out << "Число итераций напыления: " << i_nap << endl;

		avg1_time += time1;
		avg1_temp += T1;
		avg_nap += i_nap;

		//2 ЭТАП (Отжиг)
		time1 = 0;
		create_rate_catalog(chain_bool, T1);
		int i_ot = 0;
		for ( i_ot = 1;  time1 < init_values.time_of_nap; i_ot++) {
			choose_event(chain_int, chain_bool, false);
			T1 = (init_values.temperature1 - init_values.temperature2) * exp((-1) * time1 / init_values.tay) + init_values.temperature2;
			create_rate_catalog(chain_bool, T1);
		}

		avg_ot += i_ot;

		f2out << "ОТЖИГ:" << endl;
		f2out << "Длины цепочек:" << endl;
		file_chain_output(chain_bool, f2out); //подсчет распределения
		f2out << "Температура после отжига: " << T1 << endl;
		f2out << "Время отжига: " << time1 << endl;
		f2out << "Число итераций отжига: " << i_ot << endl;
		f2out << endl;

		avg2_time += time1;
		avg2_temp += T1;

	}

	//Вывод результатов в файл
	file_distribution_output(f2out); 
    file_output_description(ofstream& f2out)

	unsigned int end_time = clock(); // конечное время
	search_time = end_time - start_time; // искомое время
    
	
	ofstream file_time;
	file_time.open("run_time.txt", ios_base::app);
	file_time << s << "\t\t" << search_time / double(init_values.iterations_of_exp) << "ms" << endl;
	file_time.close();
	

	delete[] chain_bool;
	delete[] chain_int;
	delete[] mas_rate;
	delete[] events;

	

	f2out.close();
	return 0;
}


