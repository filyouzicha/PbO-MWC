#include "heuristic.h"
#include "settings.h"

int perform_random_walk = 1;
double random_walk_prob = 0.01;
double random_add_prob = 33;
double random_swap_prob = 67;
double random_drop_prob = 100;
double update_weight_prob = 0.01;
double paws_sp = 0;

int perform_update_weight = 1;
double inc_v = 0.5;
double swt_p = 0.1;
double swt_q = 0;
double averaged_crafted_weight = 0;
double sum_crafted_weight = 0;
double swt_threshold = 200;
double div_prob = 0.001;
double p_prop = 0.1;

double CO = 1.09;

int perform_perturb_step = 1;
double pert_prob = 0.0001;

int perform_restart_step = 1;
double res_prob = 0.0001;

double CO_D = 1.0;
double CO_S = 1.0;

void print_solution() {
    for (int v=0; v<Max_Vtx; v++) {
        if (TTbest[v])
            cout << v+1 << ",";
    }
    cout << endl;
}

int breaking_tie_by_random(int Index[], int length, int Node[]) {
    return Index[rand()%length];
}

int breaking_tie_by_timestamp(int Index[], int length, int Node[]) {
    int k = 0;
    k = Index[0];
    long long oldest_time = time_stamp[Node[k]];
    int index, node;
    long long time;
    for (int j=1; j<length; j++) {
        index = Index[j];
        node = Node[index];
        time = time_stamp[node];
        if (time < oldest_time) {
            oldest_time = time;
            k = index;
        }
    }
    return k;
}

int breaking_tie_by_degree(int Index[], int length, int Node[]) {
    int k = 0;
    k = Index[0];
    int max_degree = neighbor_len[Node[k]];
    int index, node;
    int degree;
    for (int j=1; j<length; j++) {
        index = Index[j];
        node = Node[index];
        degree = neighbor_len[node];
        if (degree > max_degree) {
            max_degree = degree;
            k = index;
        }
    }
    return k;
}
void smooth_weight_paws(long long w) {
    for (int i=0; i<Max_Vtx; i++) {
        if (We_penalty[i] >= w) {
            We_penalty[i] -= w;
        }
    }
    return;
}

//outside current clique
void increase_weight(long long w) {
    for (int i=0; i<Max_Vtx; i++) {
        if (0 == vectex[i]) {
            We_penalty[i] += w;
        }
    }
    return;
}
void increase_penalty_weight(long long w) {
    int i=0;
    for (i=0; i<len; i++) {
        We_penalty[cruset[i]] += w;
        //cout << We_penalty[cruset[i]] << endl;
    }
    return;
}
void update_weight_PAWS_scheme() {
    //double prob = rng.nextClosed();
    int prob = rand()%MY_RAND_MAX_INT;
    if (prob < paws_sp*100000) { // default: 0.8
        smooth_weight_paws(1); 
    }
    else {
        increase_penalty_weight(1);
    }
    return;
}

// update scheme: SWT
void smooth_weight_swt() {
    for (int i=0; i<Max_Vtx; i++) {
        sum_crafted_weight -= We_penalty[i];
        We_penalty[i] = swt_p*We_penalty[i] + swt_q*averaged_crafted_weight; // default: swt_p=0.1, swt_q=0
        sum_crafted_weight += We_penalty[i];
    }
    averaged_crafted_weight = sum_crafted_weight/Max_Vtx;
    //cout << averaged_crafted_weight << endl;
    return;
}
void increase_penalty_weight_swt(double w) {
    int i=0;
    for (i=0; i<len; i++) {
        We_penalty[cruset[i]] += w;
        sum_crafted_weight += w;
    }
    averaged_crafted_weight = sum_crafted_weight/Max_Vtx;
    //cout << averaged_crafted_weight << endl;
    return;
}
void update_weight_SWT_scheme() {
    if (averaged_crafted_weight <= swt_threshold) {
        increase_penalty_weight_swt(inc_v);
    }
    else {
        smooth_weight_swt(); 
    }
}

void update_hash_expand_constant(int m) {
    return;
}
void update_hash_backtract_constant(int m1) {
    return;
}
void update_hash_plateau_constant(int m, int m1) {
    return;
}

void update_hash_expand_dynamic(int m) {
    ptr_to_hashed_clique->update_hash_wrt_add(m);
    return;
}
void update_hash_backtract_dynamic(int m1) {
    ptr_to_hashed_clique->update_hash_wrt_drop(m1);
    return;
}
void update_hash_plateau_dynamic(int m, int m1) {
    ptr_to_hashed_clique->update_hash_wrt_add(m);
    ptr_to_hashed_clique->update_hash_wrt_drop(m1);
    return;
}

bool is_forbidden_tabu(int node) {
    return (tabuin[node] > Iter ? true:false);
}

bool is_forbidden_edge(int node) {
    int forbidden_num = 0;
    for (int i=0; i<len; i++) {
        if(edge_forbidden[node][cruset[i]] > 0) {
            forbidden_num ++;
        }
    }
    //if (double(double(forbidden_num*1.1) - double(len)) > 0.0 || tabuin[node] > Iter)
    if (double(double(forbidden_num*CO) - len) > 0.0)
        return true;
    return false;
}

// return false: not forbidden; true: forbidden
bool is_forbidden_cc(int node) {
    return (0 == conf_change[node] ? true:false);
}

bool is_forbidden_cc_tabu(int node) {
    return ((0 == conf_change[node] && tabuin[node] > Iter) ? true:false);
}

bool is_forbidden_cc_or_tabu(int node) {
    return ((0 == conf_change[node] || tabuin[node] > Iter) ? true:false);
}
// select the highest weight: We - We_penalty
void select_var_from_C1_crafted_weight_BMS(int& l1, double& w1) {
    int kkk=BMS, count=0, i=0, m=0, n=0;
    double wmn=0;
    if (len1 <= BMS) {
        select_var_from_C1_crafted_weight_traverse(l1, w1);
        return;
    }
    for (count=0; count<kkk; count++) {
        i = rand()%len1;
        m = C1[i];
        if (is_forbidden_ptr(m)) {
            continue;
        }
        n = BC[m];
        //wmn = We_penalty[n] - We_penalty[m];
        //wmn = We_penalty[n];
        wmn = (We[m]-p_prop*We_penalty[m]) - (We[n]-p_prop*We_penalty[n]);
        //wmn = (We[m]-We_penalty[m]) - (We[n]-We_penalty[n]);
        //wmn = (-We_penalty[m]) - (-We_penalty[n]);
        if (wmn > w1) {
            l1 = 0;
            w1 = wmn;
            FC1[l1++] = i;
        }
        else if (wmn == w1) {
            FC1[l1++] = i;
        }
    }
    return;
}

void select_var_from_C1_crafted_weight_traverse(int& l1, double& w1) {
    int m=0, n=0,i=0;
    double wmn=0;
    w1 = -9223372036854775807;
    for (i=0; i<len1; i++) {
        m = C1[i];
        if (!is_forbidden_ptr(m)) {
            n = BC[m];
            //wmn = We_penalty[n] - We_penalty[m];
            //wmn = We_penalty[n];
            //if (We_penalty[m] || We_penalty[n]) {
            //    cout << "error" << endl;
            //}
            wmn = (We[m]-p_prop*We_penalty[m]) - (We[n]-p_prop*We_penalty[n]);
            //wmn = (We[m]-We_penalty[m]) - (We[n]-We_penalty[n]);
            //wmn = (-We_penalty[m]) - (-We_penalty[n]);
            if (wmn > w1) {
                l1 = 0;
                w1 = wmn;
                FC1[l1++] = i;
            }
            else if (wmn == w1) {
                FC1[l1++] = i;
            }
        }
    }
    return;
}

void select_var_from_C1_BMS (int& l1, int& l2, long long& w1, long long& w2) {
    int kkk = BMS, i = 0, count = 0, m = 0, n = 0;
    long long wmn = 0;

    if (len1 <= BMS) {
        select_var_from_C1_traverse(l1, l2, w1, w2);
        return;
    }
    for (count=0; count<kkk; count++) {
        i = rand()%len1;
        m = C1[i]; // add
        n = BC[m]; // delete
        wmn = We[m] - We[n];
        if (!is_forbidden_ptr(m)) { // find the nodes that lead highest weight-increase for the current clique.
            if (wmn > w1) {
                l1 = 0;
                w1 = wmn;
                FC1[l1++] = i;
            }
            else if (wmn == w1) {
                FC1[l1++] = i;
            }
        }
        else {
            if (wmn > w2) {
                l2 = 0;
                w2 = wmn;
                TC1[l2++] = i;
            }
            else if (wmn == w2) {
                TC1[l2++] = i;
            }
        }
    }
    return;
}

void select_var_from_C1_traverse(int& l1, int& l2, long long& w1, long long& w2) {
    int m = 0, n = 0, i = 0;
    long long wmn = 0;
    for (i=0; i<len1; i++) {
        m = C1[i];
        n = BC[m];
        wmn = We[m] - We[n];
        if (!is_forbidden_ptr(m)) { // find the nodes that lead highest weight-increase for the current clique.
            if (wmn > w1) {
                l1 = 0;
                w1 = wmn;
                FC1[l1++] = i;
            }
            else if (wmn == w1) {
                FC1[l1++] = i;
            }
        }
        else {
            if (wmn > w2) {
                l2 = 0;
                w2 = wmn;
                TC1[l2++] = i;
            }
            else if (wmn == w2) {
                TC1[l2++] = i;
            }
        }
    }
    return;
}

int expand_tabu_with_scc(int SelN) {
    int i, j, k, k1, l, am, m, n, n1, node2;

    m = C0[ SelN ];
    cruset[ len++ ] = m;
    vectex[ m ] = 1;
    Wf = Wf + We[ m ];

    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;
    //time_stamp[m] = Iter;
    time_stamp[m] = Num_Iter;

    for(i=0; i<neighbor_len[m]; i++){
        node2 = neighbor[m][i];
        tabuin[node2] = Iter;
    }
    // modify by cy
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( funch[ n ] == 1 )
        {
            //delete vertex that is not adjacent to m
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            //add vertex to C1 that is only not adjacent to m, adjacent to other vertices in cruset
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;
        }
        else if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    } 
    //restart
    //update_hash_expand_ptr(m);

    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}
int expand_tabu(int SelN) {
    int i, j, k, k1, l, am, m, n, n1;

    m = C0[ SelN ];
    cruset[ len++ ] = m;
    vectex[ m ] = 1;
    Wf = Wf + We[ m ];

    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;
    //time_stamp[m] = Iter;
    time_stamp[m] = Num_Iter;

    // modify by cy
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( funch[ n ] == 1 )
        {
            //delete vertex that is not adjacent to m
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            //add vertex to C1 that is only not adjacent to m, adjacent to other vertices in cruset
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;
        }
        else if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    } 
    //restart
    //update_hash_expand_ptr(m);

    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}

int expand_edge(int SelN) {
    int i, j, k, k1, l, am, m, n, n1;

    m = C0[ SelN ]; // the node is m
    cruset[ len++ ] = m; // add m into the current set
    vectex[ m ] = 1; // set the flag?
    Wf = Wf + We[ m ]; // Wf is the weight of the current clique, i,e, weight found. update it.

    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;

    int node2;	
  
    time_stamp[m] = Num_Iter;
    for(i=0;i<neighbor_len[m];i++) {
        node2 = neighbor[m][i];
        edge_forbidden[m][node2] = edge_forbidden[node2][m] = 0;
        //conf_change[node2] = 1;
    }

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++; // WYY: funch[n] traces the number of nodes that are in the current clique

        if( funch[ n ] == 1 )
        {   // WYY: remove n from C0
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            // put it into C1
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;

            BC[ n ] = m; // WYY: BC[n] = m denotes that n is Being Connected by m.
        }
        else if( funch[ n ] == 2 )
        {
            // remove n it from C1
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    }

    //restart
    
    //update_hash_expand_ptr(m);

    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}
int expand_scc(int SelN) {
    int i, j, k, k1, l, am, m, n, n1;

    m = C0[ SelN ]; // the node is m
    cruset[ len++ ] = m; // add m into the current set
    vectex[ m ] = 1; // set the flag?
    Wf = Wf + We[ m ]; // Wf is the weight of the current clique, i,e, weight found. update it.

    /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
    /*if(DEBUG) {
      printf("\nin expand");
      printf("\nadd node %d", m);
      }*/
    //neighbor_add(m);

    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;

    int node2;	
    //for(i=0;i<Max_Vtx;i++)
    //    temp_array[i]=0;

    //conf_change[m] = 0; // keep this node from being removed immediately
    time_stamp[m] = Num_Iter;
    for(i=0;i<neighbor_len[m];i++){
        node2 = neighbor[m][i];
        //modify by cy
        conf_change[node2] = 1;
    }

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++; // WYY: funch[n] traces the number of nodes that are in the current clique

        if( funch[ n ] == 1 )
        {   // WYY: remove n from C0
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            // put it into C1
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;

            BC[ n ] = m; // WYY: BC[n] = m denotes that n is Being Connected by m.
        }
        else if( funch[ n ] == 2 )
        {
            // remove n it from C1
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    }

    //restart
    
    //update_hash_expand_ptr(m);

    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}

int expand_cc(int SelN) {
    int i, j, k, k1, l, am, m, n, n1;

    m = C0[ SelN ]; // the node is m
    cruset[ len++ ] = m; // add m into the current set
    vectex[ m ] = 1; // set the flag?
    Wf = Wf + We[ m ]; // Wf is the weight of the current clique, i,e, weight found. update it.

    len0--;
    n1 = C0[ len0 ];
    k1 = address[ m ];
    C0[ k1 ] = n1;
    address[ n1 ] = k1;

    int node2;	
    time_stamp[m] = Num_Iter;
    for(i=0;i<neighbor_len[m];i++){
        node2 = neighbor[m][i];
        //modify by cy
        conf_change[node2] = 1;
    }

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++; // WYY: funch[n] traces the number of nodes that are in the current clique

        if( funch[ n ] == 1 )
        {   // WYY: remove n from C0
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            // put it into C1
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;

            BC[ n ] = m; // WYY: BC[n] = m denotes that n is Being Connected by m.
        }
        else if( funch[ n ] == 2 )
        {
            // remove n it from C1
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    }
    //restart
    //update_hash_expand_ptr(m);

    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}

int backtract_dynamic_tabu(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + ceil(CO_D*(len0 + 2));//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_edge(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    len--;
    cruset[ SelN ] = cruset[ len ];
    //conf_change[m1] = 0; 
    int node2=0;
    for (i=0; i<len; i++) {
        node2 = cruset[i];
        edge_forbidden[node2][m1] = edge_forbidden[m1][node2] = 1;
    }
    edge_forbidden[m1][m1] = 0;
   
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_tabu_tabul3_C1(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    tabuin[ m1 ] = Iter + DTABUL + len1;
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_tabu_tabul3_randomC1(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    tabuin[ m1 ] = Iter + DTABUL + randomInt(len1+2);
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_tabu_tabul2(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    tabuin[ m1 ] = Iter + STABUL;// + randomInt(len0+2);
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_tabu_tabul3(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    tabuin[ m1 ] = Iter + DTABUL;// + randomInt(len0+2);
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}
int backtract_tabu_co_d(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    //tabuin[ m1 ] = Iter + TABUL;
    tabuin[ m1 ] = Iter + TABUL + CO_D*(len1+2);
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}

int backtract_tabu(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;
    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL;
    //tabuin[ m1 ] = Iter + TABUL + randomInt(len0 + 2);//(neighbor_len[m1]-len + 2);
    len--;
    cruset[ SelN ] = cruset[ len ];
    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;
    time_stamp[m1] = Num_Iter; //modify by cy

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}

int backtract_scc(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;

    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;

    len--;
    cruset[ SelN ] = cruset[ len ];

    /* WYY: functions of neighborhood updating */
    //neighbor_drop(m1);
    conf_change[m1] = 0;
    time_stamp[m1] = Num_Iter; //modify by cy

    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}

int backtract_cc(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1, node2;
    if( SelN == -1 )
        return -1;

    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;

    len--;
    cruset[ SelN ] = cruset[ len ];

    /* WYY: functions of neighborhood updating */
    conf_change[m1] = 0;
    time_stamp[m1] = Num_Iter; //modify by cy

    /*for(i=0;i<neighbor_len[m];i++){
      node2 = neighbor[m][i];
    //modify by cy
    conf_change[node2] = 1;
    }*/

    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
}
int backtract_scc_tabu(int SelN) {
    int i, j, k, l, m, m1, n, k1, n1;
    if( SelN == -1 )
        return -1;

    m1 = cruset[ SelN ];
    Wf = Wf - We[ m1 ];
    vectex[ m1 ] = 0;

    len--;
    cruset[ SelN ] = cruset[ len ];

    conf_change[m1] = 0;
    time_stamp[m1] = Num_Iter; //modify by cy
    tabuin[ m1 ] = Iter + TABUL;

    C0[ len0 ] = m1;
    address[ m1 ] = len0;
    len0++;

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_backtract_ptr(m1);
    return 0;
}

int plateau_tabu_constant(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    //cout << "len1 = " << len1 << endl;
    tabuin[ m1 ] = Iter + TABUL; + randomInt( len1+2 );
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_dynamic_tabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + ceil(CO_S*( len1+2 ));
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}

int plateau_tabu_tabul2_randomC1(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];
    tabuin[ m1 ] = Iter + STABUL + randomInt(len1+2);
    
    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        
    }
    return 1;   
}
int plateau_tabu_tabul2(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    tabuin[ m1 ] = Iter + STABUL ;
    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        
    }
    return 1;   
}
int plateau_tabu_tabul2_C1(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    tabuin[ m1 ] = Iter + STABUL + len1;
    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        
    }
    return 1;   
}
int plateau_tabu_co(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL + CO*(len1+2);
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_tabu_new(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + STABUL;
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_tabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if(0 == Edge[m1][m])
            //if(1 == edge_is(m1, m))
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m];i++){
      n=neighbor[m][i];
      temp_array[n]=1;
      }

      for(auto i:remaining_vertex) {
      if(i==m)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m ]; i++ ) {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL + randomInt( len1+2 );
    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;
    time_stamp[m] = Num_Iter;   //add
    time_stamp[m1] = Num_Iter;  //delete

    /*memset(temp_array, 0, sizeof(int)*Max_Vtx);
      for(i=0;i<neighbor_len[m1];i++){
      n=neighbor[m1][i];
      temp_array[n]=1;
      }
      for(auto i:remaining_vertex) {
      if(i==m1)continue;
      if(temp_array[i]==1)continue;
      n=i;
      */
    for( i = 0; i < adjaclen[ m1 ]; i++ ) {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}

int plateau_hscc(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        //if( edge_is(m1, m)== 1 )
        //    break;
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];
    int node2 = 0;
    for(i=0;i<neighbor_len[m];i++){
        node2 = neighbor[m][i];
        conf_change[node2] = 1;
    }
    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete
    
    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_edge(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        //if( edge_is(m1, m)== 1 )
        //    break;
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    int node2=0;
    for (i=0; i<len; i++) {
        node2 = cruset[i];
        edge_forbidden[node2][m1] = edge_forbidden[m1][node2] = 1;
    }
    edge_forbidden[m1][m1] = 0;
    
    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;
    
    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;


    /*for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
      for(i=0;i<neighbor_len[m];i++){

      n=neighbor[m][i];
      temp_array[n]=1;
      }

    //for(i=0;i<Max_Vtx;i++)
    for(auto i:remaining_vertex)
    {
    if(i==m)continue;
    if(temp_array[i]==1)continue;
    n=i;
    */   
    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    //tabuin[ m1 ] = Iter + TABUL + randomInt( len1+2 );
    /* WYY: neighborhood updating */
    /*if(DEBUG) {
      printf("\nin plateau: remove node %d", m1);
      dump_neighborhood();
      }*/
    //neighbor_drop(m1);
    //conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete
    /*for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
      for(i=0;i<neighbor_len[m1];i++){

      n=neighbor[m1][i];
      temp_array[n]=1;
      }

    //for(i=0;i<Max_Vtx;i++)
    for(auto i:remaining_vertex)
    {
    if(i==m1)continue;
    if(temp_array[i]==1)continue;
    n=i;
    */ 
    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}

int plateau_scc(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        //if( edge_is(m1, m)== 1 )
        //    break;
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    /* WYY: set the nodes in cruset as neighbors of m, and vice versa */
    /*if(DEBUG) {
      printf("\nin plateau: add node %d", m);
      dump_cur_clique();
      }*/
    //neighbor_add(m); // Attention: here, we don't change conf_change values for m's neighbors

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;


    /*for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
      for(i=0;i<neighbor_len[m];i++){

      n=neighbor[m][i];
      temp_array[n]=1;
      }

    //for(i=0;i<Max_Vtx;i++)
    for(auto i:remaining_vertex)
    {
    if(i==m)continue;
    if(temp_array[i]==1)continue;
    n=i;
    */   
    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    /* WYY: neighborhood updating */
    /*if(DEBUG) {
      printf("\nin plateau: remove node %d", m1);
      dump_neighborhood();
      }*/
    //neighbor_drop(m1);
    conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete
    /*for(i=0;i<Max_Vtx;i++)
      temp_array[i]=0;
      for(i=0;i<neighbor_len[m1];i++){

      n=neighbor[m1][i];
      temp_array[n]=1;
      }

    //for(i=0;i<Max_Vtx;i++)
    for(auto i:remaining_vertex)
    {
    if(i==m1)continue;
    if(temp_array[i]==1)continue;
    n=i;
    */ 
    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}

int plateau_cc(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti, node2;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    //conf_change[m] = 0;
    for(i=0;i<neighbor_len[m];i++){
        node2 = neighbor[m][i];
        //modify by cy
        conf_change[node2] = 1;
    }
    /*for(i=0;i<neighbor_len[m1];i++){
      node2 = neighbor[m1][i];
    //modify by cy
    conf_change[node2] = 1;
    }*/

    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }
    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}
int plateau_scc_ctabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];
    
    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL;// + randomInt( len1+2 );

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_hscc_ctabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];
    int node2 = 0;
    for(i=0;i<neighbor_len[m];i++){
        node2 = neighbor[m][i];
        conf_change[node2] = 1;
    }
    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL;// + randomInt( len1+2 );

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int plateau_scc_tabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, n, n1, mm1, ti;

    m = C1[ SelN  ];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    for(ti = 0; ti < len; ti++)
    {
        m1 = cruset[ ti ];
        if( 0 == Edge[m1][m] )
            break;
    }

    Wf = Wf + We[ m ] - We[ m1 ];

    //the expand process, put m into the current independent set
    vectex[ m ] = 1;
    cruset[ len++ ] = m;

    //delete m from C1
    k1 = address[ m ];
    len1--;
    n1 = C1[ len1 ];
    C1[ k1 ] = n1;
    address[ n1 ] = k1;

    for( i = 0; i < adjaclen[ m ]; i++ )
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            //cout << "tt k1 = " << k1 << "len0 = " << len0 << "n = " << n << "m = " << m << " m1 = " << m1 << endl;
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }        
    } 

    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL + randomInt( len1+2 );

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
        /*for( i = 0; i < Max_Vtx; i++ )
          {
          Tbest[ i ] = vectex[ i ];
          }*/
    }
    return 1;   
}
int forbidden_tabu(int node2) {
    return is_forbidden_tabu(node2)?conf_change[node2]:1;
}

int plateau_scc_with_tabu(int SelN) {
    int i, j, k, k1, l, m0, m, m1, m2, n, n1, mm1, ti, node2;
    m = C1[SelN];
    // WYY: swap(m1, m), where m is to be added, and m1 to be removed. m and m1 have no edge.
    m1 = BC[m];
    for(ti=0; ti<len; ti++) {
        m2 = cruset[ti];
        if(m2 != m1) {
            continue;
        }
        else {
            break;
        }
    }
    Wf = Wf + We[m] - We[m1];

    //the expand process, put m into the current independent set
    vectex[m] = 1;
    cruset[len++] = m;

    //delete m from C1
    k1 = address[m];
    len1--;
    n1 = C1[len1];
    C1[k1] = n1;
    address[n1] = k1;

    for(i = 0; i < adjaclen[ m ]; i++)
    {
        n = adjacMatrix[ m ][ i ];
        funch[ n ]++;
        if( (funch[ n ] == 1) && ( vectex[ n ] == 0 ) )
        {
            k1 = address[ n ];
            len0--;
            n1 = C0[ len0 ];
            C0[ k1 ] = n1;
            address[ n1 ] = k1;

            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
            BC[ n ] = m;

            //getchar();
        }
        else if( funch[ n ] == 2 )
        {
            len1--;
            n1 = C1[ len1 ];
            k1 = address[ n ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;
        }
    } 
    for(i=0;i<neighbor_len[m];i++) {
        node2 = neighbor[m][i];
        /*if (is_forbidden_tabu(node2)) {
            continue;
        }
        else {
            conf_change[node2] = 1;
        }*/

    conf_change[node2]=is_forbidden_tabu(node2)?conf_change[node2]:1;
        //conf_change[node2] = 1;
    }
     
    //the backtrack process, delete m1 from the current independent set
    vectex[ m1 ] = 0;
    tabuin[ m1 ] = Iter + TABUL + randomInt( len1+2 );

    len--;
    cruset[ ti ] = cruset[ len ];
    C1[ len1 ] = m1;
    address[ m1 ] = len1;
    len1++;

    conf_change[m1] = 0;
    
    time_stamp[m] = Num_Iter;  //add
    time_stamp[m1] = Num_Iter; //delete

    for( i = 0; i < adjaclen[ m1 ]; i++ )
    {
        n = adjacMatrix[ m1 ][ i ];
        funch[ n ]--;
        if( (funch[ n ] == 0) && (vectex[ n ] == 0) )
        {
            k1 = address[ n ];           
            len1--;
            n1 = C1[ len1 ];
            C1[ k1 ] = n1;
            address[ n1 ] = k1;

            C0[ len0 ] = n;
            address[ n ] = len0;
            len0++;
        }
        else if( funch[ n ] == 1 )
        {
            C1[ len1 ] = n;
            address[ n ] = len1;
            len1++;
        }
    }

    //restart
    //update_hash_plateau_ptr(m, m1);
    if( Wf > Wbest )
    {
        times(&finish);
        real_solve2 = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        real_solve2 = round(real_solve2 * 100)/100.0; 
        Wbest = Wf;
        len_best = len;
    }
    return 1;   
}
// WYY: C0 is the set of nodes that can be added? and C1 is the set nodes that can be swapped with?
/*int selectC0_random( )
  {
  int i, j, k, l, m;
  l = 0;
  if( len0 > 30 )
  {
  k = randomInt( len0 );
  return k;
  }
  for( i = 0; i < len0; i++ )
  {
  k = C0[ i ];
  if( tabuin[ k ] <= Iter )
  TC1[ l++ ] = i;
  }
  if( l == 0 )
  return -1;
  else
  {
  k = randomInt( l );
  k = TC1[ k ];
  return k;
  }
  }
  */
int selectC0_random() 
{
    if (0 != len0) {
        return randomInt(len0); 
    }
    else 
        return -1;
}
int get_degree(int node) {
    return neighbor_len[node]; 
}
int get_weight(int node) {
    return We[node];
}
int selectC0_greedy()
{
    int i, l = 0, w1 = 0, k;
    if (0 == len0)
        return -1;
    for( i = 0; i < len0; i++ )
    {
        k = C0[ i ];
        if (select_by_greedy_ptr(k) > w1) {
            l = 0;
            w1 = select_by_greedy_ptr(k);
            FC1[l++] = i;
        }
        else if (select_by_greedy_ptr(k) == w1) {
            FC1[l++] = i;
        }
    }
    if (l > 1) {
        k = randomInt(l);
        k = FC1[k];
        return k;
    }
    else {
        return FC1[0];
    }
}
//inline void clique_construction(int &l)
inline void clique_construction()
{
    int am = 0;
    while( 1 )
    {
        if (0 == len) {
            am = selectC0_random();
        }
        else {
            am = selectC0_ptr();
        }
        if( am != -1 )
        {
            Num_Iter++;
            Iter++;
            expand_ptr( am );
        }
        else 
            break;
    }

    times(&finish);
    double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
    finish_time = round(finish_time * 100)/100.0;
    if(Wbest > lbest) {
        lbest = Wbest;
        real_solve1 = real_solve2;
        printf("c %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
        for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
        //printf("o %d\n", Max_num-lbest);
    }
    if(finish_time > time_limit){
#ifdef NDEBUG
        verify();
        print_solution();
#endif
        free_memory();
        exit(0);
    }
}

long long tabu_constant(int Max_Iter) {
    int i, j, k, l, bestlen = 0, am, am1, ti, m1, m, n; 
    double ww, ww1, ww2;
    int add_index = 0, swap_index = 0, drop_index = 0;
    Iter = 0;
    clearGamma(); 
    clique_construction(); //construction
    double prob = 0, prob1 = 0, prob2 = 0;
    int p_operation = 0;
    int total_operation_len = 0;
    //while (Iter < Max_Iter) {
    while (1) {
        if (Wbest > lbest) {
            times(&finish);
            double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
            finish_time = round(finish_time * 100)/100.0;
            lbest = Wbest;
            real_solve1 = real_solve2;
            printf("c %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
            for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
            if(finish_time>time_limit){
                printf("o %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
                verify();
                print_solution();
#endif
                free_memory();
                exit(0);
            }
        }
        else {
            if (0 == Iter % 1000) {
                times(&finish);
                double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
                finish_time = round(finish_time * 100)/100.0;
                if(finish_time>time_limit){
                    printf("o %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
                    verify();
                    print_solution();
#endif
                    free_memory();
                    exit(0);
                }
            }
        }
        // if perform random walk
        if (perform_random_walk) {
            //prob = rand()%MY_RAND_MAX_INT;
            prob = rng.nextClosed();
            //cout << "prob: " << prob << endl;
            //if (prob < random_walk_prob*100000) { // random walk default:0.01 
            if (prob < random_walk_prob) { // random walk default:0.01 
                //cout << "prob: " << prob << " random walk prob: " <<random_walk_prob << endl; 
                p_operation = rand()%100;
                if (p_operation < random_add_prob && 0 != len0) { // 33
                    Num_Iter++;
                    Iter++;
                    add_index = rand()%len0;
                    l = expand_ptr(add_index);
                    if(Wbest == Waim)
                        return Wbest;
                }
                else if (p_operation < random_swap_prob && 0 != len1) { // 34
                    Num_Iter++;
                    swap_index = rand()%len1;
                    m = C1[swap_index]; // add
                    n = BC[m]; // drop
                    if ((vectex[n] != 1) || (0 != Edge[m][n])) {
                        for(j=0; j<len; j++) {
                            k = cruset[j];
                            if(0 == Edge[m][k])
                                break;
                        }
                        BC[m] = k;
                    }
                    Iter++;
                    l = plateau_ptr(swap_index);
                    if(Wbest == Waim)
                        return Wbest;
                }
                else if (0 != len) { // 33
                    Num_Iter++;
                    Iter++;
                    drop_index = rand()%len;
                    l = backtract_ptr(drop_index);  // drop
                }
                else {
                    return Wbest;
                }
                continue;
            } // end perform random walk
        } // end if perform random walk
        
        // intensification 
        am = WselectC0();
        am1 = WselectC1();
        if ((am != -1) && (am1 != -1)) {  // score > 0
            ww = We[ C0[ am ] ]; // add_score
            ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ]; // swap_score
            if (ww > ww1) { // add
                Num_Iter++;
                Iter++;
                l = expand_ptr(am);   
                if (Wbest == Waim) {
                    return Wbest;
                }
            }
            else { // swap
                Num_Iter++;
                Iter++;
                l = plateau_ptr(am1); 
                if (Wbest == Waim) {
                    return Wbest; 
                }
            }
            continue;
        }
        else if ((am != -1) && (am1 == -1)) { // score > 0
            Num_Iter++;
            Iter++;
            l = expand_ptr(am);  // add
            if (Wbest == Waim) {
                return Wbest;
            }
            continue;
        }
        else if ((am == -1) && (am1 != -1)) {
            ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
            if (ww1 > 0) {
                Num_Iter++;
                Iter++;
                l = plateau_ptr(am1);  // swap
                if (Wbest == Waim) {
                    return Wbest;
                }
                continue;
            }
        }
        // update weight
        if (perform_update_weight) {
            update_weight_ptr();
        }
        // diversification
        //prob = rng.nextClosed();
        prob = rand()%MY_RAND_MAX_INT;
        prob1 = rand()%MY_RAND_MAX_INT;
        //update by cy 2019-07-17 
        //prob2 = rand()%MY_RAND_MAX_INT;
        prob2 = rng.nextClosed();
        if (perform_update_weight && prob < div_prob*100000) { // default: 0.05
            am1 = WselectC1_crafted_weight();
            if (am1 != -1) {
                ti = Mumi_Weigt_greedy_crafted_weight();
                if (-1 == ti)
                    return Wbest;
                m1 = cruset[ti];
                ww1 = (We[C1[am1]]-p_prop*We_penalty[C1[am1]]) - (We[BC[C1[am1]]]-p_prop*We_penalty[BC[C1[am1]]]);
                
                ww2 = -(We[m1]-p_prop*We_penalty[m1]);
                
                if (ww1 > ww2) {
                    Num_Iter++;
                    Iter++;
                    l = plateau_ptr(am1);
                }
                else {
                    Num_Iter++;
                    Iter++;
                    k = backtract_ptr(ti);
                    if (k == -1) {
                        return Wbest;
                    }
                }
            }
        }
        else if (perform_perturb_step && prob1 < pert_prob*100000) {
            int rand_add_v = 0, in_clique = 1, add_index = 0, sign_add = 0;
            while (in_clique) {
                rand_add_v = rand()%Max_Vtx;
                in_clique = vectex[rand_add_v];
            }
            int len_drop=0;
            for (int i=0; i<len; i++) {
                if (0 == Edge[rand_add_v][cruset[i]]) { // cruset[i] is not adjacent to rand_add_v 
                    FC1[len_drop++] = cruset[i];
                }
            }
            if (len_drop > 0) {
                for (int i=0; i<len_drop; i++) {
                    for (int j=0; j<len; j++) {
                        if (FC1[i] == cruset[j]) {
                            k = backtract_ptr(j);
                            break;
                        }
                    }
                }
            }
            for (int i=0; i<len0; i++) {
                if (C0[i] == rand_add_v) {
                    add_index = i;
                    sign_add = 1;
                }
            }
            if (sign_add) {
                expand_ptr(add_index);
            }
            else {
                cout << "perturb error!" << endl;
            }
        }
        //else if (perform_restart_step && prob2 < res_prob*100000) {
        else if (perform_restart_step && prob2 < res_prob) {
            //cout << "restart" << endl;
            return Wbest;
        }
        else {
            if (am1 != -1) {
                ti = Mumi_Weigt_greedy();
                m1 = cruset[ ti ];
                ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
                ww2 = - We[ m1 ];
                if( ww1 > ww2 ) {
                    //cout << "sssswap" << endl;
                    Num_Iter++;
                    Iter++;
                    l = plateau_ptr( am1 );  // swap
                }
                else {
                    //cout << "drop" << endl;
                    Num_Iter++;
                    Iter++;
                    k = backtract_ptr(ti);  // drop
                    if( k == -1 )
                        return Wbest;
                }
            }
            else {
                ti = Mumi_Weigt_ptr();
                Num_Iter++;
                Iter++;
                k = backtract_ptr(ti);  // drop
                if( k == -1 )
                    return Wbest;
            }
        }
    } // end while
    return Wbest;
}

inline void drop_clear()
{
    while(len != 0) backtract_ptr(0);
}
inline int rand_index_in_clique()
{
    if(len) return rand() % len;
    return -1;
}
/*
   long long tabu_dynamic(int Max_Iter) {
   int i, j, k, l, bestlen = 0, am, am1, ti, m1, m, n; 
   long long ww, ww1, ww2;
   Iter = 0;
   clearGamma();
//restart
#ifdef clique_hash_mode
last_step_improved = 1;
#endif
int add_index = 0;
    int p_search = 0;
    int p_operation = 0;
    int total_operation_len = 0;
    while(1)
    {
        if(len == 0)
        {
            clique_construction();
            //restart
#ifdef clique_hash_mode
            last_step_improved = 1;
#endif
        }

        if(Wbest > lbest) {
            times(&finish);
            double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
            finish_time = round(finish_time * 100)/100.0;
            lbest = Wbest;
            real_solve1 = real_solve2;
            printf("c %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
            for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
            if(finish_time>time_limit){
#ifdef NDEBUG
                verify();
                print_solution();
#endif
                free_memory();
                exit(0);
            } 
        }
        else if(0 == Iter%1000){
            times(&finish);
            double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
            finish_time = round(finish_time * 100)/100.0;
            if(finish_time>time_limit){
#ifdef NDEBUG
                verify();
                print_solution();
#endif
                free_memory();
                exit(0);
            }
        }
        // random walk
        if (perform_random_walk) {
            p_search = rand()%MY_RAND_MAX_INT;
            if (p_search < random_walk_prob*100000) { // random walk 
                p_operation = rand()%100;
                if (p_operation < random_add_prob && 0 != len0) {
                    Num_Iter++;
                    add_index = rand()%len0;
                    l = expand_ptr(add_index);
                    Iter++;
                    if(Wbest == Waim)
                        return Wbest;
                }
                else if (p_operation < random_swap_prob && 0 != len1) {
                    Num_Iter++;
                    int swap_index = rand()%len1;
                    m = C1[swap_index];
                    n = BC[m];
                    if((vectex[n] != 1) || (0 != Edge[m][n])) {
                        for(j=0; j<len; j++) {
                            k = cruset[j];
                            if(0 == Edge[m][k])
                                break;
                        }
                        BC[m] = k;
                    }
                    l = plateau_ptr(swap_index);
                    Iter++;
                    if(Wbest == Waim)
                        return Wbest;
                }
                else if (0 != len) {
                    Num_Iter++;
                    int drop_index = rand()%len;
                    l = backtract_ptr(drop_index);  // drop
                    Iter++;
                }
                else {
                    return Wbest;
                }
                continue;
            }
        }
        // intensification 
        am = WselectC0();
        if(len == 1)//when |C|=1, forbid swapping, base-1
        {
            am1 = -1;
        }
        else
        {
            am1 = WselectC1();
        }

#ifdef first_improved_revisit_restart_mode
        int is_first_improved_step = 0;
#endif

        ww = We[ C0[ am ] ];
        ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
        //should test whether restart is needed

        if(am != -1)
        {
#ifdef first_improved_revisit_restart_mode
            if(!last_step_improved)
            { 
                is_first_improved_step = 1;
            }
#endif
            if( am1 != -1 )
            {
                if( ww > ww1 )
                {
                    Num_Iter++;
                    expand_ptr( am );  
                    Iter++;
                }
                else
                {
                    Num_Iter++;
                    plateau_ptr( am1 );
                    Iter++;
                }
            }
            else if( am1 == -1 )//add-set is not empty; swap-set is empty
            {
                Num_Iter++;
                expand_ptr( am );
                Iter++;
            }
            //restart
#ifdef clique_hash_mode
            last_step_improved = 1;
#endif
        }
        else
        {
            //restart
#ifdef clique_hash_mode
            if(am1 == -1 || ww1 <= 0)
            {
                last_step_improved = 0;
            }
            else
            {
#ifdef first_improved_revisit_restart_mode
                if(!last_step_improved)
                {
                    is_first_improved_step = 1;
                }
#endif
                last_step_improved = 1;
            }	
#endif
            ti = Mumi_Weigt_greedy();					 
            if( am1 != -1 )//add-set is empty; swap-set is not empty
            {
                m1 = cruset[ ti ];
                ww2 = - We[ m1 ];
                if( ww1 > ww2 )
                {
                    Num_Iter++;
                    plateau_ptr( am1 );
                    Iter++;
                }
                else
                {
                    Num_Iter++;
                    backtract_ptr(ti);
                    Iter++;
                }
            }
            else//add-set is empty; swap-set is empty
            {
                Num_Iter++;
                backtract_ptr(Mumi_Weigt_ptr());
                Iter++;
            }
        }

#ifdef first_improved_revisit_restart_mode
        if(is_first_improved_step)
        {
            ptr_to_hashed_clique->mark_hash_entry();
            if(ptr_to_hashed_clique->curr_hash_entry_marked_too_frequently())
            {
                return Wbest;
            }
        }
#endif
    }

    return Wbest;
}
*/
long long tabu_consecutive(int Max_Iter) {
    int i, j, k, l, bestlen = 0, am, am1, ti, m1; 
    long long ww, ww1, ww2;
    Iter = 0;
    clearGamma(); 
    clique_construction();

    while( Iter < Max_Iter )
    {
        am = WselectC0();
        am1 = WselectC1();
        if( (am != -1) && (am1 != -1) )
        {
            ww = We[ C0[ am ] ];
            ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];

            if( ww > ww1 )
            {
                Num_Iter++;
                l = expand_ptr( am );  // add 
                Iter++;
                if( Wbest == Waim )
                    return Wbest;
            }
            else
            {
                Num_Iter++;
                l = plateau_ptr( am1 ); // swap
                Iter++;
                if( Wbest == Waim )
                    return Wbest; 
            }
        }
        else if( (am != -1) && (am1 == -1) )
        {
            //modify by cy
            Num_Iter++;
            l = expand_ptr( am );  // add
            Iter++;
            if( Wbest == Waim )
                return Wbest;
        }
        else if( (am == -1) && (am1 != -1) )
        {
            ti = Mumi_Weigt_greedy();
            m1 = cruset[ ti ];
            ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
            ww2 = - We[ m1 ];
            if( ww1 > ww2 )
            {
                Num_Iter++;
                l = plateau_ptr( am1 );  // swap
                Iter++;
                if( Wbest == Waim )
                    return Wbest; 
            }
            else
            {
                Num_Iter++;
                k = backtract_ptr(ti);  // drop
                Iter++;
                if( k == -1 )
                    return Wbest;
            }
        }
        else if( (am == -1) && (am1 == -1) )
        {
            ti = Mumi_Weigt_ptr();
            Num_Iter++;
            k = backtract_ptr(ti);  // drop
            Iter++;
            if( k == -1 )
                return Wbest;
        }
        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        finish_time = round(finish_time * 100)/100.0;
        if(Wbest > lbest){
            lbest = Wbest;
            real_solve1 = real_solve2;
            printf("c %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
            for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
            Iter = 0; // consecutive
            //printf("o %d\n", Max_num-lbest);
        }
        if(finish_time>time_limit){
            printf("o %.2f %lld %lld\n", real_solve1,lbest,Num_Iter);
            //printf("o %d\n", Max_num-lbest);
#ifdef NDEBUG
            verify();
            print_solution();
#endif
            free_memory();
            exit(0);
        } 
    }
    return Wbest;
}
/*
long long tabu_dynamic_and_constant(int Max_Iter) {
    int i, j, k, l, bestlen = 0, am, am1, ti, m1;
    long long ww, ww1, ww2;
    Iter = 0;
    clearGamma();
    //restart
#ifdef clique_hash_mode
    //if (restart_dynamic) {
    last_step_improved = 1;
    //}
#endif

    while(1)
    {
        if (Iter >= Max_Iter) {
            return Wbest;
        }
        if(len == 0)
        {
            clique_construction();
            //restart
#ifdef clique_hash_mode
            last_step_improved = 1;
#endif
        }
        am = WselectC0();

        if(len == 1)//when |C|=1, forbid swapping, base-1
        {
            am1 = -1;
        }
        else
        {
            am1 = WselectC1();
        }

#ifdef first_improved_revisit_restart_mode
        int is_first_improved_step = 0;
#endif

        ww = We[ C0[ am ] ];
        ww1 = We[ C1[ am1 ] ] - We[ BC[ C1[ am1 ] ] ];
        //should test whether restart is needed

        if(am != -1)
        {
#ifdef first_improved_revisit_restart_mode
            if(!last_step_improved)
            { 
                is_first_improved_step = 1;
            }
#endif
            if( am1 != -1 )
            {
                if( ww > ww1 )
                {
                    Num_Iter++;
                    expand_ptr( am );  
                    Iter++;
                }
                else
                {
                    Num_Iter++;
                    plateau_ptr( am1 );
                    Iter++;
                }
            }
            else if( am1 == -1 )//add-set is not empty; swap-set is empty
            {
                Num_Iter++;
                expand_ptr( am );
                Iter++;
            }
            //restart
#ifdef clique_hash_mode
            last_step_improved = 1;
#endif
        }
        else
        {
            //restart
#ifdef clique_hash_mode
            if(am1 == -1 || ww1 <= 0)
            {
                last_step_improved = 0;
            }
            else
            {
#ifdef first_improved_revisit_restart_mode
                if(!last_step_improved)
                {
                    is_first_improved_step = 1;
                }
#endif
                last_step_improved = 1;
            }	
#endif
            ti = Mumi_Weigt_greedy();					 
            if( am1 != -1 )//add-set is empty; swap-set is not empty
            {
                m1 = cruset[ ti ];
                ww2 = - We[ m1 ];
                if( ww1 > ww2 )
                {
                    Num_Iter++;
                    plateau_ptr( am1 );
                    Iter++;
                }
                else
                {
                    Num_Iter++;
                    backtract_ptr(ti);		               
                    Iter++;
                }
            }
            else//add-set is empty; swap-set is empty
            {
                Num_Iter++;
                backtract_ptr(Mumi_Weigt_ptr());
                Iter++;
            }
        }

#ifdef first_improved_revisit_restart_mode
        if(is_first_improved_step)
        {
            ptr_to_hashed_clique->mark_hash_entry();
            if(ptr_to_hashed_clique->curr_hash_entry_marked_too_frequently())
            {
                return Wbest;
                //drop_clear();
                //Iter++;
                //continue;
            }
        }
#endif
        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        finish_time = round(finish_time * 100)/100.0;
        if(Wbest>lbest){
            lbest=Wbest;
            real_solve1=real_solve2;
            printf("c	%.2f	%lld  %lld\n", real_solve1,lbest,Num_Iter);
            //printf("o %d\n", Max_num-lbest);
        }
        if(finish_time>time_limit){
            free_memory();
            exit(0);
        } 
    }

    return Wbest;
}
*/
