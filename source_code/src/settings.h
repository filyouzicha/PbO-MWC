#ifndef _SETTINGS_H
#define _SETTINGS_H

#include "basis.h"

extern int (*breaking_tie_ptr) (int [], int, int []);
extern bool (*is_forbidden_ptr) (int node);
extern void (*select_var_from_C1_ptr) (int& l1, int& l2, long long& w1, long long& w2);
extern void (*select_var_from_C1_crafted_weight_ptr) (int& l1, double& w1);
extern int (*expand_ptr) (int SelN);
extern int (*backtract_ptr) (int SelN); 
extern int (*plateau_ptr) (int SelN);
extern long long (*tabu_ptr) (int Max_Iter);
extern void (*update_weight_ptr) ();
extern int (*selectC0_ptr) ();
extern int (*Mumi_Weigt_ptr) ();
extern int (*Mumi_Weigt_crafted_ptr) ();

//extern void (*update_hash_expand_ptr) (int m);
//extern void (*update_hash_backtract_ptr) (int m1);
//extern void (*update_hash_plateau_ptr) (int m, int m1);
//extern void (*new_hash_ptr) (int num);
//extern void (*delete_hash_ptr) ();
//extern void (*reset_hash_ptr) ();
//extern void (new_hash_constant) (int num);
//extern void (delete_hash_constant) ();
//extern void (reset_hash_constant) ();
//extern void (new_hash_dynamic) (int num);
//extern void (delete_hash_dynamic) ();
//extern void (reset_hash_dynamic) ();

extern int (*select_by_greedy_ptr) (int node);

void default_algorithm_settings();

extern long long Waim;
extern int TABUL;
extern int BMS;
extern int len_improve;
extern double drop_random_prob;
extern double random_walk_prob;
extern double update_weight_prob;
extern int perform_random_walk;
extern int perform_update_weight;
extern int perform_perturb_step;
extern int perform_restart_step;
extern bool restart_dynamic;
extern double swt_p;
extern double swt_q;
extern double div_prob;
extern double p_prop;
extern double inc_v;
extern double pert_prob;
extern double res_prob;


#endif
