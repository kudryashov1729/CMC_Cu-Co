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

double time1; //время
double T1; //температура
char s[200]; //имя выводимого файла 


bool init(int argc, char* argv[], ofstream& f2out) {
	if (argc == 6) {
		init_values.iterations_of_exp = atoi(argv[1]);
		init_values.tay = atoi(argv[2]);
		init_values.temperature1 = atoi(argv[3]);
		init_values.temperature2 = atoi(argv[4]);
		init_values.time_of_nap = atoi(argv[5]);
		strcpy(s, "exp=");
		strcat(s, argv[1]);
		strcat(s, "_tay=");
		strcat(s, argv[2]);
		strcat(s, "_T1=");
		strcat(s, argv[3]);
		strcat(s, "_T2=");
		strcat(s, argv[4]);
		strcat(s, "_time=");
		strcat(s, argv[5]);
		strcat(s, ".txt");
		f2out.open(s, ios::out);
		return false;
	}
	std::cout << endl << "Input:\n iterations_of_exp\n tay\n temperature1\n temperature2\n time_of_nap" << endl;
	return true;
}

int main(int argc, char* argv[])
{
	ofstream f2out;
	if (init(argc, argv, f2out)) {
		delete[] mas_rate;
		delete[] events;
		return 0;
	}
	setlocale(LC_ALL, "Russian");

	unsigned int start_time = clock();

	int * chain_int = new int[LEN_OF_CHAIN];
	bool * chain_bool = new bool[LEN_OF_CHAIN];

	for (int i = 0; i < LEN_OF_CHAIN; i++) {
		chain_lenth_distribution[i] = 0;
	}


	avg_atoms = 0;
	avg_ot = 0;
	avg_nap = 0;
	double avg1_time = 0;
	double avg2_time = 0;
	double avg1_temp = 0;
	double avg2_temp = 0;

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

	unsigned int end_time = clock(); // конечное время
	unsigned int search_time = end_time - start_time; // искомое время
	std::cout << "Время выполнения программы " << search_time << " мс";
	f2out << "Время выполнения программы " << search_time << " мс" << endl;

	f2out << -1 << endl;
	f2out << "T1=" << init_values.temperature1 << "K ";
	f2out << "T2=" << init_values.temperature2 << "K ";
	f2out << "Exp=" << init_values.iterations_of_exp << endl;
	f2out << "Адатомов: " << avg_atoms / (init_values.iterations_of_exp) << "/100";
	f2out << " Тау=" << init_values.tay << " секунд" << endl;

	f2out << "Ср. вр. нап: " << avg1_time / (init_values.iterations_of_exp) << "сек" << endl;
	f2out << "Ср. ч-ло ит. нап: " << avg_nap / (init_values.iterations_of_exp) << endl;
	f2out << "Ср. т-ра после нап: " << avg1_temp / (init_values.iterations_of_exp) << "K" << endl;

	f2out << "Ср. вр. отжига: " << avg2_time / (init_values.iterations_of_exp) << "cек" << endl;
	f2out << "Ср. ч-ло ит. отжига: " << avg_ot / (init_values.iterations_of_exp) << endl;
	f2out << "Ср. т-ра после отжига: " << avg2_temp / (init_values.iterations_of_exp) << endl;
	
	ofstream file_time;
	file_time.open("run_time.txt", ios_base::app);
	file_time << s << "\t\t" << search_time / double(init_values.iterations_of_exp) << "ms" << endl;
	file_time.close();
	

	delete[] chain_bool;
	delete[] chain_int;
	delete[] mas_rate;
	delete[] events;

	

	f2out.close();
	//std::system("python py_vis.py");
	return 0;
}


