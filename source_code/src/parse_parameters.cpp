#include "parse_parameters.h"
#include "settings.h"
#include "heuristic.h"

bool parse_parameters(int argc, char **argv) {
    int i = 0;
    int temp_para = 0;

    default_algorithm_settings();

    for (i=1; i<argc; i++) {
        if (0 == strcmp(argv[i], "-seed")) {
            i++;
            if (i >= argc) 
                return false;
            sscanf(argv[i], "%d", &seed);
        }
        else if (0 == strcmp(argv[i], "-inst")) {
            i++;
            if (i >= argc) 
                return false;
            File_Name = argv[i];
        }
        else if (0 == strcmp(argv[i], "-cutoff")) {
            i++;
            if (i >= argc) 
                return false;
            sscanf(argv[i], "%lf", &time_limit);
        }
        else if (0 == strcmp(argv[i], "-reduct")) {
            i++;
            if (i >= argc) 
                return false;
            sscanf(argv[i], "%d", &temp_para);
            if (0 == temp_para) {
                label_reduction = false;
            }
            else {
                label_reduction = true;
            }
        }
        else if (0 == strcmp(argv[i], "-bt")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                breaking_tie_ptr = breaking_tie_by_random;
            }
            else if (0 == strcmp(argv[i], "1")) {
                breaking_tie_ptr = breaking_tie_by_timestamp;
            }
            else if (0 == strcmp(argv[i], "2")) {
                breaking_tie_ptr = breaking_tie_by_degree;
            }
        }
        else if (0 == strcmp(argv[i], "-tabu")) {
            i++;
            if (i >= argc)
                return false;
            if ( 0 == strcmp(argv[i], "0")) {//scc
                is_forbidden_ptr = is_forbidden_cc;
                expand_ptr = expand_scc;
                backtract_ptr = backtract_scc;
                plateau_ptr = plateau_scc;
            }
            else if (0 == strcmp(argv[i], "1")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu;
                backtract_ptr = backtract_tabu;
                plateau_ptr = plateau_tabu;
            }
            else if (0 == strcmp(argv[i], "2")) { //tabu+scc
                is_forbidden_ptr = is_forbidden_cc_tabu;
                expand_ptr = expand_scc; //expand_scc_tabu is equal to expand_scc
                backtract_ptr = backtract_scc_tabu;
                plateau_ptr = plateau_scc_tabu;
            }
            else if (0 == strcmp(argv[i], "3")) {// cc
                is_forbidden_ptr = is_forbidden_cc;
                expand_ptr = expand_cc;
                backtract_ptr = backtract_cc;
                plateau_ptr = plateau_cc;
            }
            else if (0 == strcmp(argv[i], "4")) { // constant tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu;
                backtract_ptr = backtract_tabu;
                plateau_ptr = plateau_tabu_constant;
            }
            else if (0 == strcmp(argv[i], "5")) { // half scc
                is_forbidden_ptr = is_forbidden_cc;
                expand_ptr = expand_scc;
                backtract_ptr = backtract_scc;
                plateau_ptr = plateau_hscc;
            }
            else if (0 == strcmp(argv[i], "6")) { //constant tabu + half scc
                is_forbidden_ptr = is_forbidden_cc_tabu;
                expand_ptr = expand_scc;  //expand_hscc_ctabu is equal to expand_scc
                backtract_ptr = backtract_scc_tabu;
                plateau_ptr = plateau_hscc_ctabu;
            }
            else if (0 == strcmp(argv[i], "7")) { //constant tabu + scc
                is_forbidden_ptr = is_forbidden_cc_tabu;
                expand_ptr = expand_scc;  //expand_scc_ctabu is equal to expand_scc
                backtract_ptr = backtract_scc_tabu;
                plateau_ptr = plateau_scc_ctabu;
            }
            else if (0 == strcmp(argv[i], "8")) { //tabu or scc
                is_forbidden_ptr = is_forbidden_cc_or_tabu;
                expand_ptr = expand_scc; //expand_scc_tabu is equal to expand_scc
                backtract_ptr = backtract_scc_tabu;
                plateau_ptr = plateau_scc_tabu;
            }
            else if (0 == strcmp(argv[i], "9")) { //update scc using tabu
                is_forbidden_ptr = is_forbidden_cc;
                expand_ptr = expand_scc; //expand_scc_tabu is equal to expand_scc
                backtract_ptr = backtract_scc_tabu;
                plateau_ptr = plateau_scc_with_tabu;
            }
            else if (0 == strcmp(argv[i], "10")) { //update tabu using scc
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc; //expand_scc_tabu is equal to expand_scc
                backtract_ptr = backtract_tabu;
                plateau_ptr = plateau_tabu;
            }
            else if (0 == strcmp(argv[i], "11")) { //update tabu using scc
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc; //expand_scc_tabu is equal to expand_scc
                backtract_ptr = backtract_dynamic_tabu;
                plateau_ptr = plateau_dynamic_tabu;
            }
            else if (0 == strcmp(argv[i], "12")) { //update tabu using scc
                is_forbidden_ptr = is_forbidden_edge;
                expand_ptr = expand_edge; //expand_edge 
                backtract_ptr = backtract_edge;
                plateau_ptr = plateau_edge;
            }
            else if (0 == strcmp(argv[i], "13")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                //expand_ptr = expand_tabu;
                backtract_ptr = backtract_tabu;
                plateau_ptr = plateau_tabu_new;
            }
            else if (0 == strcmp(argv[i], "15")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu;
                plateau_ptr = plateau_tabu_co;
            }
            else if (0 == strcmp(argv[i], "16")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu_co_d;
                plateau_ptr = plateau_tabu_co;
            }
            else if (0 == strcmp(argv[i], "17")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu_tabul3;
                plateau_ptr = plateau_tabu_tabul2;
            }
            else if (0 == strcmp(argv[i], "18")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu_tabul3_randomC1;
                plateau_ptr = plateau_tabu_tabul2;
            }
            else if (0 == strcmp(argv[i], "19")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu_tabul3;
                plateau_ptr = plateau_tabu_tabul2_randomC1;
            }
            else if (0 == strcmp(argv[i], "20")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                backtract_ptr = backtract_tabu_tabul3;
                plateau_ptr = plateau_tabu_tabul2_C1;
            }
            else if (0 == strcmp(argv[i], "21")) {//tabu
                is_forbidden_ptr = is_forbidden_tabu;
                expand_ptr = expand_tabu_with_scc;
                //backtract_ptr = backtract_tabu_tabul2;
                backtract_ptr = backtract_tabu;
                //plateau_ptr = plateau_tabu_tabul2_randomC1;
                plateau_ptr = plateau_tabu;
            }
        }
        else if (0 == strcmp(argv[i], "-bms")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                select_var_from_C1_ptr = select_var_from_C1_traverse;         
                select_var_from_C1_crafted_weight_ptr = select_var_from_C1_crafted_weight_traverse;         
            }
            else if (0 == strcmp(argv[i], "1")) {
                select_var_from_C1_ptr = select_var_from_C1_BMS;
                select_var_from_C1_crafted_weight_ptr = select_var_from_C1_crafted_weight_BMS;
            }
        }
        else if (0 == strcmp(argv[i], "-res")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                tabu_ptr = tabu_constant;    
            }
            /*else if (0 == strcmp(argv[i], "1")) {
                tabu_ptr = tabu_dynamic;
                restart_dynamic = true;
                update_hash_expand_ptr = update_hash_expand_dynamic;
                update_hash_backtract_ptr = update_hash_backtract_dynamic;
                update_hash_plateau_ptr = update_hash_plateau_dynamic;
                new_hash_ptr = new_hash_dynamic;
                delete_hash_ptr = delete_hash_dynamic;
                reset_hash_ptr = reset_hash_dynamic;
            }*/
            else if (0 == strcmp(argv[i], "2")) {
                tabu_ptr = tabu_consecutive;
            }
            /*else if (0 == strcmp(argv[i], "3")) {
                tabu_ptr = tabu_dynamic_and_constant;
                restart_dynamic = true;
                update_hash_expand_ptr = update_hash_expand_dynamic;
                update_hash_backtract_ptr = update_hash_backtract_dynamic;
                update_hash_plateau_ptr = update_hash_plateau_dynamic;
                new_hash_ptr = new_hash_dynamic;
                delete_hash_ptr = delete_hash_dynamic;
                reset_hash_ptr = reset_hash_dynamic;
            }*/
        }
        else if (0 == strcmp(argv[i], "-update_weight_scheme")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                update_weight_ptr = update_weight_SWT_scheme;
            }
            else if (0 == strcmp(argv[i], "1")) {
                update_weight_ptr = update_weight_PAWS_scheme;
            }
        }
        else if (0 == strcmp(argv[i], "-cons")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                selectC0_ptr = selectC0_random;
            }
            else if (0 == strcmp(argv[i], "1")) {
                selectC0_ptr = selectC0_greedy;
                select_by_greedy_ptr = get_weight;
            }
            else if (0 == strcmp(argv[i], "2")) {
                selectC0_ptr = selectC0_greedy;
                select_by_greedy_ptr = get_degree;
            }
        }
        else if (0 == strcmp(argv[i], "-drop")) {
            i++;
            if (i >= argc)
                return false;
            if (0 == strcmp(argv[i], "0")) {
                Mumi_Weigt_ptr = Mumi_Weigt_greedy;
                Mumi_Weigt_crafted_ptr = Mumi_Weigt_greedy_crafted_weight;
            }
            else if (0 == strcmp(argv[i], "1")) {
                Mumi_Weigt_ptr = Mumi_Weigt_probability;
                Mumi_Weigt_crafted_ptr = Mumi_Weigt_probability_crafted_weight;
            }
            else if (0 == strcmp(argv[i], "2")) {
                Mumi_Weigt_ptr = Mumi_Weigt_random;
                Mumi_Weigt_crafted_ptr = Mumi_Weigt_random;
            }
        }
        else if (0 == strcmp(argv[i], "-Waim")) {
            i++;
            if (i >= argc)
                return false;
            Waim = atoll(argv[i]);
            //cout << "Waim: " << Waim << endl;
        }
        else if (0 == strcmp(argv[i], "-tabul")) {
            i++;
            if (i >= argc)
                return false;
            TABUL = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-co_d")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &CO_D);
        }
        else if (0 == strcmp(argv[i], "-co")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &CO);
        }
        else if (0 == strcmp(argv[i], "-co_s")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &CO_S);
        }
        else if (0 == strcmp(argv[i], "-bn")) {
            i++;
            if (i >= argc)
                return false;
            BMS = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-len")) {
            i++;
            if (i >= argc)
                return false;
            len_improve = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-lenconsecutive")) {
            i++;
            if (i >= argc)
                return false;
            len_improve = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-tabul2")) {
            i++;
            if (i >= argc)
                return false;
            STABUL = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-tabul3")) {
            i++;
            if (i >= argc)
                return false;
            DTABUL = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-rd_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &drop_random_prob);
        }
        else if (0 == strcmp(argv[i], "-rw_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &random_walk_prob);
        }
        else if (0 == strcmp(argv[i], "-uw_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &update_weight_prob);
        }
        else if (0 == strcmp(argv[i], "-p_rw")) {
            i++;
            if (i >= argc)
                return false;
            perform_random_walk = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-p_uw")) {
            i++;
            if (i >= argc)
                return false;
            perform_update_weight = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-p_ps")) {
            i++;
            if (i >= argc)
                return false;
            perform_perturb_step = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-p_rs")) {
            i++;
            if (i >= argc)
                return false;
            perform_restart_step = atoi(argv[i]);
        }
        else if (0 == strcmp(argv[i], "-paws_sp")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &paws_sp);
        }
        else if (0 == strcmp(argv[i], "-swt_p")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &swt_p);
        }
        else if (0 == strcmp(argv[i], "-swt_q")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &swt_q);
        }
        else if (0 == strcmp(argv[i], "-div_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &div_prob);
        }
        else if (0 == strcmp(argv[i], "-p_prop")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &p_prop);
        }
        else if (0 == strcmp(argv[i], "-inc_v")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &inc_v);
        }
        else if (0 == strcmp(argv[i], "-pert_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &pert_prob);
        }
        else if (0 == strcmp(argv[i], "-res_prob")) {
            i++;
            if (i >= argc)
                return false;
            sscanf(argv[i], "%lf", &res_prob);
        }
    }
    return true;
}
