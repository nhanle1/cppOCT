// stdafx.cpp : source file that includes just the standard includes
// octProc.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information

#include "stdafx.h"

// TODO: reference any additional headers you need in STDAFX.H
// and not in this file


//output image file for checkings
/*
ofstream outputfile("img.txt");
if (outputfile.is_open()) {
for (int i = 0; i < (2* nX); i++) {
for (int j = 0; j < nZ; j++) {
outputfile << abs(specC[i*nZ + j]) << '\t';
if (j == nZ - 1) {
outputfile << '\n';
*/



//ofstream outputfile("ED.bin", ios::binary);
//if (outputfile.is_open()) {
//	for (size_t i = 0; i < nZ*nXR; i++) {
//		outputfile.write((char*)(EDdata.data() + i), sizeof(complex<double>));
//	}

//}
//outputfile.close();


/*
//Linear interpolation function

for (int i = 0; i < nXR; i++) {
dinterp1(buffer.data()+i*nZ, nZ, KES, nZ, tbuf.data() + i * nZ);
}
copy(tbuf.begin(), tbuf.end(), specC.begin());

double interpolate(double x, vector<pair<double, double> > table) {
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > table.back().first) return 0; //or INF depending on usage
	if (x < table[0].first) return 0;
	vector<pair<double, double> >::iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(table.begin(), table.end(), make_pair(x, INF));
	// Corner case
	if (it == table.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first) / (it->first - it2->first);
}
*/