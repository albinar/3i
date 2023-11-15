#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include "mystring.h"
#include <Rinternals.h>
#include <R.h>
//#include <cassert>

using namespace std;
///// string manupulation part

void removeAllCharsFromString(string &str, char c)
{
	for ( ; str.find(c) != string::npos; str = str.replace(str.find(c), 1, "")) ;
}

void replaceAll(string &str, char c1, char c2)
{
	char tmp[2] = {c2, 0};
	for ( ; str.find(c1) != string::npos; str = str.replace(str.find(c1), 1, tmp)) ;
}

//// replaces 0 with 1 and 1 with 0 in input string.
//// that's needed to make input key numbers more sensible.
void inversZeroOne(string &str)
{
	replaceAll(str, '1', '_');
	replaceAll(str, '0', '1');
	replaceAll(str, '_', '0');
}

///// file reading part

////////////// not needed yet //////////////
void readAllProportions(const char *name)
{
}

void readPvals(const char *signs_file_name, 
		const char *pval_file_name, 
		double *pvals, 
		int *s_values, 
		char **signs, 
		int &signs_row_num,
		char **proteins)
{
	FILE *f_p = fopen(pval_file_name, "r");
	FILE *f_s = fopen(signs_file_name, "r");
	char line[2000];
	char *ltmp;
	char *save_ptr;
	ltmp = fgets(line, 2000, f_p);
	if (!ltmp)
		Rprintf("dataRead.cpp:readPvals: Code should have beed asserted here 1!!!\n");
	ltmp = fgets(line, 2000, f_s);
	if (!ltmp)
		Rprintf("dataRead.cpp:readPvals: Code should have beed asserted here 2!!!\n");

	char *tmp;
	my_strtok_r(line, ",", &save_ptr);
	int protein_it = 0;
	while (1)
	{
		tmp = my_strtok_r(NULL, ",", &save_ptr);
		if (!tmp)
			break;
		string pr = tmp;
		removeAllCharsFromString(pr, '\n');
		removeAllCharsFromString(pr, '"');
		my_strcpy(proteins[protein_it], pr.c_str());
		protein_it++;
	}

	int count = 0;

	signs_row_num = 0;
	int s_values_it = 0;
	while (!feof(f_p) && !feof(f_s))
	{
		line[0] = 0;
		ltmp = fgets(line, 2000, f_p);
		if (!line[0])
			break;

		char *tmp;
		tmp = my_strtok_r(line, ",\n\r", &save_ptr);
		tmp = my_strtok_r(NULL, ",\n\r", &save_ptr);
		double pvalue = strtod(tmp, NULL);
		pvals[signs_row_num] = pvalue;

		line[0] = 0;
		ltmp = fgets(line, 2000, f_s);
		if (!line[0])
			break;

		tmp = my_strtok_r(line, ",", &save_ptr);
		string name = tmp;
		removeAllCharsFromString(name, '"');
		my_strcpy(signs[signs_row_num], name.c_str());

		while (1)
		{
			tmp = my_strtok_r(NULL, ",", &save_ptr);
			if (!tmp)
				break;
			string vl = tmp;
			removeAllCharsFromString(vl, '\n');
			removeAllCharsFromString(vl, '"');
			s_values[s_values_it] = atoi(vl.c_str());
			s_values_it++;
		}


		signs_row_num++;
		count++;
		/*
		fprintf(stderr, "%s\t%lf\n", tmp, pvalue);
		if (count == 10)
			return;
		*/
	}
	fclose(f_p);
	fclose(f_s);
}
