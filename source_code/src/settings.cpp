#include "settings.h"
#include "heuristic.h"

int (*breaking_tie_ptr) (int [], int, int []);

bool (*is_forbidden_ptr) (int node);

void (*select_var_from_C1_ptr) (int& l1, int& l2, long long& w1, long long& w2);
void (*select_var_from_C1_crafted_weight_ptr) (int& l1, double& w1);

int (*expand_ptr) (int SelN);
int (*backtract_ptr) (int SelN); 
int (*plateau_ptr) (int SelN); 
void (*update_weight_ptr) ();
long long (*tabu_ptr) (int Max_Iter);

int (*selectC0_ptr) ();

int (*Mumi_Weigt_ptr) ();
int (*Mumi_Weigt_crafted_ptr) ();

//void (*update_hash_expand_ptr) (int m);
//void (*update_hash_backtract_ptr) (int m1);
//void (*update_hash_plateau_ptr) (int m, int m1);

//void (*new_hash_ptr) (int num);
//void (*delete_hash_ptr) ();
//void (*reset_hash_ptr) ();
int (*select_by_greedy_ptr) (int node);
void algorithm_settings_structured_instances() {
    label_reduction = false;
    breaking_tie_ptr = breaking_tie_by_timestamp;
    is_forbidden_ptr = is_forbidden_cc;
    select_var_from_C1_ptr = select_var_from_C1_BMS;
    select_var_from_C1_crafted_weight_ptr = select_var_from_C1_crafted_weight_BMS;
    expand_ptr = expand_scc;
    backtract_ptr = backtract_scc;
    plateau_ptr = plateau_scc;
    tabu_ptr = tabu_constant;
    update_weight_ptr = update_weight_SWT_scheme;
    restart_dynamic = false;
    selectC0_ptr = selectC0_greedy;
    Mumi_Weigt_ptr = Mumi_Weigt_greedy;
    Mumi_Weigt_crafted_ptr = Mumi_Weigt_greedy_crafted_weight;
    TABUL = 7;
    BMS = 100;
    len_improve = 4000;
    drop_random_prob = 20;
    CO_D = 1.0;
    CO_S = 1.0;
    //update_hash_expand_ptr = update_hash_expand_constant;
    //update_hash_backtract_ptr = update_hash_backtract_constant;
    //update_hash_plateau_ptr = update_hash_plateau_constant;
    //new_hash_ptr = new_hash_constant; 
    //delete_hash_ptr = delete_hash_constant; 
    //reset_hash_ptr = reset_hash_constant;
    select_by_greedy_ptr = get_weight;
    
    perform_random_walk = 0;
    random_walk_prob = 0.01;
    
    perform_update_weight = 0;
    update_weight_prob = 0.01;
    paws_sp = 0.2;
    swt_p = 0.1;
    swt_q = 0;
    div_prob = 0.05;
    p_prop = 0.1;
    inc_v = 0.5;
    
    perform_perturb_step = 0;
    pert_prob = 0.001;
    
    perform_restart_step = 0;
    res_prob = 0.001;
}

void default_algorithm_settings() //Default settings
{
	algorithm_settings_structured_instances();
}
