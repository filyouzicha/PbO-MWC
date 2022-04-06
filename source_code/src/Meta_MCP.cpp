// =====================================================================================
//
//       Filename:  .cpp
//
//    Description:  This is a solver for weighted maximum clique problem based on SCC and BMS
//
//        Version:  1.0
//        Created:  
//       Revision:  none
//       Compiler:  g++
//
//         Author:  

//         For example：
// =====================================================================================



#include "basis.h"
#include "parse_parameters.h"
#include "settings.h"
#include "cliqueHash.h"
#include "heuristic.h"
#define	MAXV	15000000
#define MAXE	80000000

int BMS=100;
//long long Max_num = 3179769846193;
long long Max_num = 3000000000000;
tms start, finish;
double time_limit;
double real;
long long lbest;
long long Num_Iter = 0;
struct Edge1{
    int v1;
    int v2;
};
int** Edge;
int** edge_forbidden;
//Edge1 edge[MAXE];  
Edge1* edge = NULL;
//int v_degree_tmp[MAXV];
int* v_degree_tmp = NULL;

//int adjaclen[MAXV];
int* adjaclen = NULL;
int** adjacMatrix = NULL;

double real_solve1=-1;
double real_solve2=-1;
vector< vector<int> > neighbor;
int* neighbor_len = NULL;
int* conf_change = NULL;

long long* time_stamp = NULL;

int* temp_array = NULL;

char * File_Name;

int* vectex = NULL;
int* funch = NULL;
int* address = NULL;
long long* tabuin = NULL;

int Max_Vtx = 0,  Max_Iter; 
//int f;
//int fbest;

int* cruset = NULL;
int len;
int tm1;
int tm2;
int* C0 = NULL; // the candidates node for ADD operation?
int* C1 = NULL; // the Candidate nodes for SWAP operation?
long long* We = NULL; // Weights of nodes
double* We_penalty = NULL; // Weights of nodes
int* BC = NULL;

int len0; // the length of C0
int len1; // the length of C1

int* TC1 = NULL; // Temporal candidate nodes?

long long Iter; // the number of iterations taken
int TABUL = 7;
int STABUL = 7;
int DTABUL = 7;
long long Wf;
long long Wbest;
int* FC1 = NULL;
int* Tbest = NULL;
int* TTbest = NULL;
long long Waim;

int len_best = 0;
int len_W;
int Iteration[ 100 ];
double time_used[ 100 ];
int len_used[ 100 ];
int W_used[ 100 ];
char outfilename[30];
int len_improve = 4000;
int len_time;
int Wmode;
double drop_random_prob = 0.2;
//int TABUL0 = 5;

//add by cy
long long* vertex_neighbor_weight = NULL;
bool is_reduction = true;
Remaining_vertex remaining_vertex;
bool label_reduction = false;
//#ifdef GRAPH_REDUCTION 
int size_threshold = 100000;
int last_reduction_best = 0;

//add by cy restart rrwl
CliqueHash *ptr_to_hashed_clique;
int last_step_improved = 1;
bool restart_dynamic = false;
//#endif

/************************************************************************/
/*   WYY: Variables for Configuration Checking                    */
/************************************************************************/
/* neighbor[i][j] = n means that i and n are connceted, and n is the jth 
 *  neighbor of n.
 */

/* time_stamp[m]=j means the conf_change value of node m recently 
 * changes at jth iteration.
 */

//struct  timeval start;
//struct  timeval end;

void dump_neighborhood();
#define DEBUG 0
///////////////////////////
int edge_is(int m, int n)
{
    int i;
    int index;
    for(i=0;i<neighbor_len[m];i++){
        index=neighbor[m][i];
        if(index==n)return 0;
    }
    return 1;
}

/*void new_hash_constant (int num) {
    return;
}
void delete_hash_constant () {
    return;
}
void reset_hash_constant () {
    return;
}
void new_hash_dynamic (int num) {
    ptr_to_hashed_clique = new CliqueHash(num);
    return;
}
void delete_hash_dynamic () {
    delete ptr_to_hashed_clique;
    return;
}
void reset_hash_dynamic () {
    ptr_to_hashed_clique->reset_hash_entry();
    return;
}
*/
void malloc_memory(int countv, int counte) {
    int max_vtx = countv + 2;
    int max_edge = counte + 2;
    /*v_degree_tmp = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == v_degree_tmp) {
        cout << "malloc v_degree_tmp failed." << endl;
        exit(0);
    }*/
    adjaclen = (int *)calloc(max_vtx, sizeof(int));
    if (NULL == adjaclen) {
        cout << "malloc adjaclen failed." << endl;
        exit(0);
    }
    neighbor_len = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == neighbor_len) {
        cout << "malloc neighbor_len failed." << endl;
        exit(0);
    }
    conf_change = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == conf_change) {
        cout << "malloc conf_change failed." << endl;
        exit(0);
    }
    time_stamp = (long long *)malloc(sizeof(long long) * max_vtx);
    if (NULL == time_stamp) {
        cout << "malloc time_stamp failed." << endl;
        exit(0);
    }
    temp_array = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == temp_array) {
        cout << "malloc temp_array failed." << endl;
        exit(0);
    }
    vectex = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == vectex) {
        cout << "malloc vectex failed." << endl;
        exit(0);
    }
    funch = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == funch) {
        cout << "malloc funch failed." << endl;
        exit(0);
    }
    address = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == address) {
        cout << "malloc address failed." << endl;
        exit(0);
    }
    tabuin = (long long *)malloc(sizeof(long long) * max_vtx);
    if (NULL == tabuin) {
        cout << "malloc tabuin failed." << endl;
        exit(0);
    }
    cruset = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == cruset) {
        cout << "malloc cruset failed." << endl;
        exit(0);
    }
    C0 = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == C0) {
        cout << "malloc C0 failed." << endl;
        exit(0);
    }
    C1 = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == C1) {
        cout << "malloc C1 failed." << endl;
        exit(0);
    }
    We = (long long *)malloc(sizeof(long long) * max_vtx);
    if (NULL == We) {
        cout << "malloc We failed." << endl;
        exit(0);
    }
    We_penalty = (double *)calloc(max_vtx, sizeof(double));
    if (NULL == We_penalty) {
        cout << "malloc We_penalty failed." << endl;
        exit(0);
    }
    BC = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == BC) {
        cout << "malloc BC failed." << endl;
        exit(0);
    }
    TC1 = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == TC1) {
        cout << "malloc TC1 failed." << endl;
        exit(0);
    }
    FC1 = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == FC1) {
        cout << "malloc FC1 failed." << endl;
        exit(0);
    }
    Tbest = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == Tbest) {
        cout << "malloc Tbest failed." << endl;
        exit(0);
    }
    TTbest = (int *)malloc(sizeof(int) * max_vtx);
    if (NULL == TTbest) {
        cout << "malloc TTbest failed." << endl;
        exit(0);
    }
    vertex_neighbor_weight = (long long *)malloc(sizeof(long long) * max_vtx);
    if (NULL == vertex_neighbor_weight) {
        cout << "malloc vertex_neighbor_weight failed." << endl;
        exit(0);
    }
    edge = (Edge1 *)malloc(sizeof(Edge1) * max_edge);
    if (NULL == edge) {
        cout << "malloc edge failed." << endl;
        exit(0);
    }
    Edge = (int**) malloc(sizeof(int*) * max_vtx);
    if (NULL == Edge) {
        cout << "malloc Edge failed." << endl;
        exit(0);
    }
    edge_forbidden = (int**) malloc(sizeof(int*) * max_vtx);
    if (NULL == edge_forbidden) {
        cout << "malloc edge_forbidden failed." << endl;
        exit(0);
    }
    adjacMatrix = (int**) malloc(sizeof(int*) * max_vtx);
    if (NULL == adjacMatrix) {
        cout << "malloc adjacMatrix failed." << endl;
        exit(0);
    }
    //restart
    //new_hash_ptr(max_vtx);
}

void free_memory() {
    //free(v_degree_tmp);
    free(adjaclen);
    free(neighbor_len);
    free(conf_change);
    free(time_stamp);
    free(temp_array);
    free(vectex);
    free(funch);
    free(address);
    free(tabuin);
    free(cruset);
    free(C0);
    free(C1);
    free(We);
    free(We_penalty);
    free(BC);
    free(TC1);
    free(FC1);
    free(Tbest);
    free(TTbest);
    free(edge);
    for (int i=0; i<Max_Vtx; i++) {
        free(edge_forbidden[i]);
    }
    free(edge_forbidden);
    //restart
    //delete_hash_ptr();
}

// section 0, initiaze
void Initializing()
{
    ifstream FIC;
    FILE *fp;
    FIC.open(File_Name);
    if ( FIC.fail() )
    {
        cout << "### Erreur open, File_Name " << File_Name << endl;
        getchar();
        exit(0);
    }
    char StrReading[100];
    //Max_Vclique=300;
    // FIC >> StrReading;
    if ( FIC.eof() )
    {
        cout << "### Error open, File_Name " << File_Name << endl;
        exit(0);
    }
    int nb_vtx=0, nb_edg=-1, max_edg=0;
    int x1;
    long long x2;

    string line;
    istringstream is;
    string p, tmp;

    do {
        getline(FIC, line);
        is.clear();
        is.str(line);
        is >> p >> tmp >> Max_Vtx >> nb_edg;
    } while (p != "p");
    cout << Max_Vtx << " " << nb_edg << endl;
    malloc_memory(Max_Vtx, nb_edg);
    cout << "------- end malloc memory" << endl;
    //FIC >> tt >> Max_Vtx >> nb_edg;
    //neighbor=(int **)malloc(Max_Vtx*sizeof(int*));//neighbor set
    neighbor.resize(Max_Vtx);//neighbor set
    /*working_vertex.resize(Max_Vtx);
    next_working_vertex.resize(Max_Vtx);
    is_pending.resize(Max_Vtx);*/
    remaining_vertex.init(Max_Vtx);
    for( int x = 0 ; x < Max_Vtx ; x++ ) 
    {
        conf_change[x] = 1; // initialize
        time_stamp[x] = 0;
        //adjaclen[x]=0;
        neighbor_len[x]=0;
        vectex[x]  = x;
        address[x] = x;
        Edge[x] = (int*)calloc(Max_Vtx, sizeof(int));
        edge_forbidden[x] = (int*)calloc(Max_Vtx, sizeof(int));
        Edge[x][x] = 1;
        adjacMatrix[x] = (int*)calloc(Max_Vtx, sizeof(int));
#ifdef FORMAT_CLQ
        We[x] = (x+1) % Wmode + 1;
        swt_threshold += We[x];
        //We[x] = 1;
#endif
    }
    int ii;
/*#ifdef FORMAT_CLQ
    for (int x=0; x<Max_Vtx; x++) {
        //We[x] = (x+1) % Wmode + 1;
        We[x] = 1;
    }
#endif
*/
    ii = 0;
#ifdef FORMAT_WCLQ
    FIC >> tmp >> x1 >> x2;
    while ("n" == tmp) {
        x1--;
        We[x1] = x2;
        swt_threshold += x2;
        FIC >> tmp >> x1 >> x2;
    }
    x1--;
    x2--;
    max_edg++;
    neighbor_len[x1]++;
    neighbor_len[x2]++;
    vertex_neighbor_weight[x1] += We[x2];
    vertex_neighbor_weight[x2] += We[x1];
    edge[ii].v1 = x2;
    edge[ii].v2 = x1;
    ii++;
#endif
    /*for( int x = 0; x < Max_Vtx; x++ )
      {
      FIC >>tmp>> x1 >> x2;
      x1--;
      We[ x1 ] =x2;       
      }
      */
    swt_threshold = (swt_threshold/Max_Vtx);
    cout << swt_threshold << endl;
    while(FIC >> tmp >> x1 >> x2)
    {
        x1--; x2--;
        if ( x1<0 || x2<0 || x1>=Max_Vtx || x2 >=Max_Vtx )
        {
            cout << "### Error of node : x1="
                << x1 << ", x2=" << x2 << endl;
            exit(0);
        }
        max_edg++;
        neighbor_len[x1]++;
        neighbor_len[x2]++;
        vertex_neighbor_weight[x1] += We[x2];
        vertex_neighbor_weight[x2] += We[x1];
        edge[ii].v1 = x2;
        edge[ii].v2 = x1;   
        ii++;
    }
    int v;
    /*for (v=0; v<Max_Vtx; v++) {
        //adjaclen[v]=Max_Vtx-1-neighbor_len[v];
        Edge[v] = (int*)calloc(Max_Vtx, sizeof(int));
        Edge[v][v] = 1;
        adjacMatrix[v] = (int*)calloc(Max_Vtx, sizeof(int));
    }
    */

    //for (v=0; v<Max_Vtx; v++) neighbor[v]=(int *)malloc( neighbor_len[v]*sizeof(int));

    //for(v=0; v<Max_Vtx; v++) v_degree_tmp[v]=0; 
    int e,v1,v2;  
    for (e=0; e<nb_edg; e++)
    {
        v1=edge[e].v1;
        v2=edge[e].v2;
        neighbor[v1].push_back(v2);
        neighbor[v2].push_back(v1);
        //neighbor[v1][v_degree_tmp[v1]] = v2;
        //neighbor[v2][v_degree_tmp[v2]] = v1;
        //v_degree_tmp[v1]++;
        //v_degree_tmp[v2]++;
        Edge[v1][v2] = Edge[v2][v1] = 1;
    }
    int i;


    if ( 0 && max_edg != nb_edg )
    {
        cout << "### Error de lecture du graphe, nbre aretes : annonce="
            << nb_edg << ", lu=" << max_edg  << endl;
        exit(0);
    }


    for( int x = 0; x < Max_Vtx; x++ )
    {
        // We[ x ] = (x+1)%Wmode + 1;
        BC[ x ] = 0;
        //We[ x ] = 1;
        //We[ x ] = ( rand() % 10 ) + 1;
        for (int y = 0; y < Max_Vtx; y ++) {
            if (0 == Edge[x][y]) {
                adjacMatrix[x][adjaclen[x]] = y;
                adjaclen[x] ++;
            }
        }
    }
    //modify by cy
    //reset_hash_ptr();
    FIC.close();
}



// WYY
void dump_conf_change() {
    printf("\nconf_change:\n");
    for(int i = 0; i < Max_Vtx; i++) {
        printf("%4d(%d) ", i, conf_change[i]);
    }
    printf("\n");
}

// WYY
void dump_neighborhood() {
    printf("Neighborhood:\n");
    for(int i = 0; i < Max_Vtx; i++) {
        printf(": ");
        for(int j = 0; j < neighbor_len[i]; j++)
            printf("%d ", neighbor[i][j]);
        printf("\n");	
    }
    return;
}

// WYY
void dump_cur_clique() {
    return;
    int n;
    printf("\ncurrent clique includes %d nodes:", len);
    for(int i = 0; i < len; i++) {
        n = cruset[i];
        printf("%d(%d) ", n, vectex[n]);
    }
    printf("\n");
}

int randomInt( int n )
{
    return rand() % n;
}

void clearGamma()
{
    int i, j, k, l;
    memset( vectex, 0, tm1 );
    memset( funch, 0, tm1 );
    memset( tabuin, 0, Max_Vtx*sizeof(long long) );
    memset( C1, 0, tm1 );
    memset( BC, 0, tm1 );
    int index = 0;
    
    for( i = 0; i < Max_Vtx; i++ )
    {
        C0[ i ] = i;
        address[ i ] = i;
        // wyy: clear configuration Information for the next restart
        conf_change[i] = 1;
        //time_stamp[i] = 0; //modify by cy according to NuMVC
    }
    len0 = Max_Vtx;
    len1 = 0;
    len = 0;
    Wf = 0;
    Wbest = 0;
    //restart
#ifdef clique_hash_mode
    //ptr_to_hashed_clique->reset_hash_entry();
    //reset_hash_ptr();
#endif
}

// WYY: Select from C0, a node, which is not in tabu and has the max weight, 
// or satisfy the aspiration rule though in tabu.
int WselectC0() {
    int i, j, k, l1, l2, m;
    l1 = 0;
    l2 = 0;
    long long w1 = 0, w2 = 0;

    for (i=0; i<len0; i++) {
        k = C0[i];
        // WYY:store nodes that are not in tabu list and with the maximum weight, in FC1
        if (!is_forbidden_ptr(k)) {
            if (We[ k ] > w1) {
                l1 = 0;
                w1 = We[k];
                FC1[l1++] = i;
            }
            else if (We[ k ] == w1) {
                FC1[l1++] = i;
            }
        }
        else {  // WYY: stores nodes that are being in tabu but with the maximum weight, in TC1
            if (We[ k ] > w2) {
                l2 = 0;
                w2 = We[k];
                TC1[l2++] = i;
            }
            else if (We[ k ] == w2) {
                TC1[l2++] = i;
            }
        }
    }

    // WYY: to check first if the aspiration rule is applicable.
    // If not, select a nodes which have the highest weithgts; break ties randomly.
    if( (l2 > 0) && ((w2+Wf)>Wbest) && (w2>w1)) { 
        // WYY: Select the node with the oldest age
        // modify by cy
        //if (TABUL > 3)
        //    TABUL--;
        return breaking_tie_ptr(TC1, l2, C0);
        //cout << "add: yes in aspiration w2+Wf = " << w2+Wf << endl;
    }  
    else if( l1 > 0 ) {
        // WYY
        // modify by cy
        return breaking_tie_ptr(FC1, l1, C0);
    }
    else {
        return -1;
    }
}

// return the index of the node selected from C1, select the highest one
int WselectC1_crafted_weight() {
    if (0 == len1) {
        return -1;
    }
    int l1 = 0;
    double w1 = -9223372036854775807;
    select_var_from_C1_crafted_weight_ptr(l1, w1);
    if (l1 > 0) {
        return breaking_tie_ptr(FC1, l1, C1);
    }
    else {
        return -1;
    }
}

int WselectC1() {
    if( 0 == len1) {
        return -1;
    }
    int i, j, k, l, l1, l2, m, n;
    long long wmn = 0, w1, w2;
    l1 = 0;
    l2 = 0;
    w1 = -9223372036854775807;
    w2 = -9223372036854775807;
    for (i=0; i<len1; i++) {
        m = C1[i];
        n = BC[m];
        if ((1 != vectex[n]) || (0 != Edge[m][n])) {
            for (j=0; j<len; j++) {
                k = cruset[j];
                if (0 == Edge[m][k])
                    break;
            }
            BC[ m ] = k;
        }
    }
   
    select_var_from_C1_ptr(l1, l2, w1, w2);
    if ( (l2 > 0) && (w2 > w1) && ((w2+Wf) > Wbest) ) {
        //cout << "swap: yes in aspiration w2+Wf = " << w2+Wf << endl;
        return breaking_tie_ptr(TC1, l2, C1);
    }  
    else if (l1 > 0) {
        return breaking_tie_ptr(FC1, l1, C1);
    }
    else {
        return -1;
    }
}

// WYY: find nodes with minimum weight in the current clique to remove
int Mumi_Weigt_probability() {
    double p = rng.nextClosed();
    if (p < drop_random_prob) {
        //cout << "random" << endl;
        return Mumi_Weigt_random();
    }
    else {
        //cout << "greedy" << endl;
        return Mumi_Weigt_greedy();
    }
}

int Mumi_Weigt_random(){
    if (len <= 0)
        return -1;
    else 
        return randomInt(len);
}
int Mumi_Weigt_greedy() {
    int i, j, k, l1, m;
    long long w1 = 9223372036854775807;
    l1 = 0;
    // WYY: find in cruset the nodes with lowest weights, chose one of it randomly
    for( i = 0; i < len; i++ ) {
        k = cruset[ i ];
        if( We[ k ] < w1 ) {
            l1 = 0;
            w1 = We[ k ];
            FC1[ l1++ ] = i;
        }
        else if ( We[ k ] == w1 ) {
            FC1[ l1++ ] = i;
        }
    }

    if( l1 == 0 ) {
        return -1;
    }
    return FC1[randomInt(l1)];
}
int Mumi_Weigt_probability_crafted_weight() {
    double p = rng.nextClosed();
    if (p < drop_random_prob) {
        //cout << "random" << endl;
        return Mumi_Weigt_random();
    }
    else {
        //cout << "greedy" << endl;
        return Mumi_Weigt_greedy_crafted_weight();
    }
}

int Mumi_Weigt_greedy_crafted_weight() {
    int i, j, k, l1, m;
    double w1 = 9223372036854775807, ww2 = 0;
    l1 = 0;
    // WYY: find in cruset the nodes with lowest weights, chose one of it randomly
    for (i = 0; i < len; i++) {
        k = cruset[i];
        //ww2 = -We_penalty[k];
        ww2 = We[k]-p_prop*We_penalty[k];
        //ww2 = We[k]-We_penalty[k];
        //ww2 = 0-We_penalty[k];
        if (ww2 < w1) {
            l1 = 0;
            w1 = ww2;
            FC1[l1++] = i;
        }
        else if (ww2 == w1) {
            FC1[l1++] = i;
        }
    }

    if( l1 == 0 ) {
        return -1;
    }
    return FC1[randomInt(l1)];
}
void verify()
{
    int i, j, k1, k2, l, m;
    for( i = 0; i < Max_Vtx; i++ )
    {
        if( TTbest[ i ] == 1 )
        {
            for( j = i+1; j < Max_Vtx; j++ )
                if( (TTbest[ j ] == 1) && ( 0 == Edge[i][j]) )
                    cout << "hello there is something wrong" << endl;
        }
    }
    cout << "verified" << endl;
}

// WYY: Validate that the cruset is indeed a clique
void validate() {
    int i, j, k1, k2, l, m;
    for( i = 0; i < Max_Vtx; i++ )
    {
        if( vectex[ i ] == 1 )
        {
            for( j = i + 1; j < Max_Vtx; j++ )
                //if( (vectex[ j ] == 1) && ( edge_is(i, j)== 1 ) ) {
                if( (vectex[ j ] == 1) && (0 == Edge[i][j]) ) {
                    cout << "hello there is something wrong" << endl;
                    getchar();
                }
        }
    }
}

void Output()
{
    int i , j, k, l, sum; 
    FILE *fp ;
    int len = strlen(File_Name);

    strcpy(outfilename,File_Name) ;
    outfilename[len]='.';
    outfilename[len+1]='o';
    outfilename[len+2]='u';
    outfilename[len+3]='t';
    outfilename[len+4]='\0';

    fp = fopen(outfilename, "a+"); 
    for( i = 0; i < 1; i++ )
    {
        fprintf(fp, "sum = %6d, iter = %6d, len = %5d,  time = %8lf \n", 
                W_used[ i ], Iteration[ i ], len_used[ i ],  time_used[ i ] ); 
    }

    fclose(fp); // WYY
    return;

    fprintf(fp, "\n\n the total information: \n");
    int wavg, iteravg, lenbb, success;
    wavg = iteravg = lenbb = success = 0;
    int best_v = 0;
    double timeavg = 0; 
    for( i = 0; i < 100; i++ )
        if( W_used[ i ] > best_v )
        {
            best_v = W_used[ i ];  
            lenbb = len_used[ i ];
        }

    int count = 0;
    fprintf(fp, "\n The best weight value for the maximum weighted problem is %6d \n", best_v);
    for( i = 0; i < 100; i++ )
    {
        wavg = wavg + W_used[ i ];
    }  
    double twavg = (double (wavg)) / 100 ; 
    for( i = 0; i < 100; i++ )
        if( W_used[ i ] == best_v )
        {
            count++;
            iteravg = iteravg + Iteration[ i ];
            timeavg = timeavg + time_used[ i ];
        }

    iteravg =  int ( (double (iteravg)) / count );
    timeavg = timeavg / (count*1000);
    fprintf(fp, "avg_sum = %10lf, succes = %6d, len = %5d, avg_iter = %6d,  time = %8lf \n", 
            twavg, count, lenbb,  iteravg, timeavg );
    fclose(fp) ;
    return ;
}

long long Max_Tabu()
{
    int i, j, k, m;
    long long l = 0;
    lbest = 0;
    int lenbest = 0;
    times(&start);
    while(1)
    {
        l = tabu_ptr(len_improve);
        /*if (is_reduction) {
          graph_reduction_iterative(l);
          is_reduction = false;
          }*/
        if( l > lbest )
        {
            real_solve1=real_solve2;
            lbest = l; 
            len_W = len_best;      
            printf("c %.2f %lld %lld\n", real_solve1, lbest, Num_Iter);
#ifdef NDEBUG
            for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
            //printf("o %d\n", Max_num-lbest);
        }

        if( lbest >= Waim )
            return lbest;

        times(&finish);
        double finish_time = double(finish.tms_utime - start.tms_utime + finish.tms_stime - start.tms_stime)/sysconf(_SC_CLK_TCK);
        finish_time = round(finish_time * 100)/100.0;
        if(finish_time>time_limit) {
            printf("o	%.2f	%lld  %lld\n", real_solve1,lbest, Num_Iter);
           //printf("o %d\n", Max_num-lbest);
#ifdef NDEBUG
            verify();
            print_solution();
#endif
            free_memory();
            exit(0);
        }
//#ifdef GRAPH_REDUCTION
        /*if (label_reduction  && last_reduction_best < lbest) {
            if (remaining_vertex.size() > size_threshold) {
                graph_reduction(lbest);
            }
            else {
                graph_reduction_iterative(lbest);
            }
            last_reduction_best = lbest;
            if (remaining_vertex.empty()) {
                cout << "find optimal solution" << endl;
                exit(0);
            }
        }*/
//#endif
    }
    return lbest;
}

int seed;
int main(int argc, char **argv)
{
    Waim=9223372036854775807;
    //File_Name = argv[1];
    Wmode=200;
    //time_limit=atof(argv[2]);
    //seed = atoi(argv[3]);
    if (argc<7) {
        cout << "./Meta -inst <inst> -seed <seed> -cutoff <cutoff_time>" << endl;
        exit(0);
    }

    parse_parameters(argc, argv);
    printf("c	%s\n",argv[1]);
    srand(seed);
    rng.seed(seed);
    cout << "------- start Init" << endl;
    Initializing();
    cout << "------- end Init" << endl;
    tm1 = Max_Vtx*sizeof( int );
    int i;
    long long l = 0;
    l = Max_Tabu();
    
    if( l > lbest )
    {
        real_solve1=real_solve2;
        lbest = l; 
        len_W = len_best;      
        printf("c	%.2f	%lld  %lld\n", real_solve1,lbest,Num_Iter);
#ifdef NDEBUG
            for(int i = 0; i < Max_Vtx; i++) TTbest[i] = vectex[i];
#endif
        //printf("o %d\n", Max_num-lbest);
    }
#ifdef NDEBUG
            verify();
            print_solution();
#endif
    free_memory();
    cout << "finished" << endl;
    return 0;
}
