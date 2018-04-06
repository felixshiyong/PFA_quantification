#pragma once
#include <cstring>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <map>
#include <fstream>
#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <random>
#include <time.h>
#include <cmath>
#include <sstream>
using namespace std;


class protein_uniprot
{
public:
	protein_uniprot(void);
	~protein_uniprot(void);

	string IDLine;
	vector <string> AC_accessions;
	
	vector <int> FT_modifications;
	vector <int> FT_Variants;

	string protein_SEQ; 

	void initial_class();
	//// give null values of all the properties.
};

