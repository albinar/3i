#ifndef _DATA_READ_
#define _DATA_READ_

void removeAllCharsFromString(string &str, char c);
void replaceAll(string &str, char c1, char c2);
void inversZeroOne(string &str);
void readAllProportions(const char *name);
void readPvals(const char *signs_file_name, 
		const char *pval_file_name, 
		double *pvals, 
		int *s_values, 
		char **signs, 
		int &signs_row_num,
		char **proteins);
#endif
