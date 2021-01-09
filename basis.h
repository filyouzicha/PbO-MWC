#ifndef BASIS_H
#define BASIS_H

#include "mersenne.h"
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <sys/times.h>
#include <unistd.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <string.h>
#include <math.h>
#include <assert.h>
using namespace std;

//#define FORMAT_WCLQ
#define FORMAT_CLQ
//#define GRAPH_REDUCTION
#define NDEBUG
struct Remaining_vertex {
	vector<long> vertex;
	vector<vector<long>::size_type> index;

	vector<long>::iterator begin() {
		return vertex.begin();
	}
	vector<long>::iterator end() {
		return vertex.end();
	}
	void init(vector<long>::size_type vertex_size) {
		vertex.reserve(vertex_size);
		index.resize(vertex_size);
		for (vector<long>::size_type i = 0; i < vertex_size; i++) {
			vertex.push_back(i);
			index[i] = i;
		}
	}
	void remove(long v) {
		index[*vertex.rbegin()] = index[v];
		vertex[index[v]] = *vertex.rbegin();
		vertex.pop_back();
	}

	vector<long>::size_type size() {
		return vertex.size();
	}

	bool empty() {
		return vertex.empty();
	}
};
extern Mersenne rng;
extern Remaining_vertex remaining_vertex;
extern long long* We;
extern double* We_penalty;
extern vector< vector<int> > neighbor;
extern int* neighbor_len;
extern long long* vertex_neighbor_weight;// = NULL;

extern vector<long> is_pending;
extern vector<long> working_vertex;
extern vector<long> next_working_vertex;

extern bool label_reduction;
extern long long* tabuin; 
extern long long Iter;
extern int seed;
void graph_reduction(int threshold);
void graph_reduction_iterative(int threshold);

const int MY_RAND_MAX_INT = 10000000; 
#endif
