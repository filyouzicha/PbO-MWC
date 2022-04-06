#ifndef _HEURISTIC_H
#define _HEURISTIC_H

#include "basis.h"
#include "cliqueHash.h"
//#define clique_hash_mode
#define first_improved_revisit_restart_mode
//#define rand_drop_mode

extern long long* time_stamp;
extern long long Num_Iter;
extern int* conf_change;

int breaking_tie_by_random(int Index[], int length, int Node[]);
int breaking_tie_by_timestamp(int Index[], int length, int Node[]);
int breaking_tie_by_degree(int Index[], int length, int Node[]);

extern int len1;
extern int* C1;
extern int* BC;
extern int* FC1;
extern int* TC1;
extern int BMS;
void select_var_from_C1_BMS(int& l1, int& l2, long long& w1, long long& w2);
void select_var_from_C1_traverse(int& l1, int& l2, long long& w1, long long& w2);
void select_var_from_C1_crafted_weight_BMS(int& l1, double& w1);
void select_var_from_C1_crafted_weight_traverse(int& l1, double& w1);

//restart
extern CliqueHash *ptr_to_hashed_clique;
extern int last_step_improved;

extern int* funch;
extern int* C0;
extern int* cruset;
extern int len;
extern int* vectex;
extern long long Wf;
extern int len0;
extern int* address;
extern int* temp_array;
extern int Max_Vtx;
extern long long Wbest;
extern int* TTbest;
extern int len_best;
extern tms start, finish;
extern double real_solve1;
extern double real_solve2;
extern int* adjaclen;
extern int** adjacMatrix;

int expand_tabu(int SelN);
int expand_scc(int SelN);
int expand_cc(int SelN);
int expand_tabu_with_scc(int SelN);
extern int expand_edge(int SelN);

extern int TABUL;
extern int STABUL;
extern int DTABUL;
extern int Mumi_Weigt();
int backtract_tabu (int SelN);
int backtract_dynamic_tabu (int SelN);
int backtract_scc (int SelN);
int backtract_cc (int SelN);
int backtract_scc_tabu (int SelN);
extern int backtract_tabu_co_d (int SelN);
extern int backtract_edge (int SelN);

extern int backtract_tabu_tabul2 (int SelN);
extern int backtract_tabu_tabul3 (int SelN);
extern int backtract_tabu_tabul3_C1 (int SelN);
extern int backtract_tabu_tabul3_randomC1 (int SelN);

extern double CO;
extern double CO_D;
extern double CO_S;
extern int randomInt (int n);
extern int edge_is (int m, int n);
extern int** Edge;
int plateau_tabu_constant(int SelN);
int plateau_dynamic_tabu(int SelN);
int plateau_tabu(int SelN);
int plateau_hscc(int SelN);
int plateau_scc(int SelN);
int plateau_cc(int SelN);
int plateau_scc_tabu(int SelN);
int plateau_hscc_ctabu(int SelN);
int plateau_scc_ctabu(int SelN);
int plateau_scc_with_tabu(int SelN);
extern int plateau_edge(int SelN);
extern int plateau_tabu_new(int SelN);
extern int plateau_tabu_co(int SelN);

extern int plateau_tabu_tabul2(int SelN);
extern int plateau_tabu_tabul2_C1(int SelN);
extern int plateau_tabu_tabul2_randomC1(int SelN);

extern int** edge_forbidden;

bool is_forbidden_tabu(int node);
bool is_forbidden_cc(int node);
bool is_forbidden_cc_tabu(int node);
bool is_forbidden_cc_or_tabu(int node);
extern bool is_forbidden_edge(int node);

extern void clearGamma();
extern int selectC0();
extern long long Waim;
extern long long lbest;
extern void free_memory();
extern int WselectC0();
extern int WselectC1();
extern int WselectC1_crafted_weight();
extern double time_limit;

extern int selectC0_random();
extern int selectC0_greedy();

extern int Mumi_Weigt_random();
extern int Mumi_Weigt_greedy();
extern int Mumi_Weigt_probability();
extern int Mumi_Weigt_greedy_crafted_weight();
extern int Mumi_Weigt_probability_crafted_weight();
long long tabu_constant(int Max_Iter);
//long long tabu_dynamic(int Max_Iter);
long long tabu_consecutive(int Max_Iter);
//long long tabu_dynamic_and_constant(int Max_Iter);

extern int get_degree(int node);
extern int get_weight(int node);

extern long long Max_num;

void update_hash_expand_constant(int m);
void update_hash_expand_dynamic(int m);
void update_hash_backtract_constant(int m1);
void update_hash_backtract_dynamic(int m1);
void update_hash_plateau_constant(int m, int m1);
void update_hash_plateau_dynamic(int m, int m1);

extern void verify();
extern void print_solution(); 
extern double paws_sp;
void update_weight_SWT_scheme();
void update_weight_PAWS_scheme();
extern double sum_crafted_weight;
extern double averaged_crafted_weight;
extern double swt_threshold;
#endif
