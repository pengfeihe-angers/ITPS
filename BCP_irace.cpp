//============================================================================
// Name        : BCP_PR.cpp
// Author      : Jintong
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>
#include <time.h>
#include <ctime>
#include <vector>
#include <cmath>
#include <math.h>
#include <complex>
#include <assert.h>
#include <iomanip>
#include <sys/times.h>
#include <unistd.h>

using namespace std;

#define min(a,b) ((a)<(b)?(a):(b))

// Variables time
struct tms glo_start;
struct tms glo_end;
struct tms midt;
double clockTicksPerSecond;
double startTimeSeconds;
double cpuTimeTheBest;



/*data structure of a solution*/
struct Structure_Solution{
	int *permutation;
	int *permutationNew;
	int *weightcount;
	double cbmp;
	long double EF2;
	long double EF3;
};
Structure_Solution *CurrentS, *BestS, *Thebest,*LastS;

/*structure of element bucket*/
struct list_bucket{
	int node[2];
	int flag_empty;
	int flag_tail;
	list_bucket *next;
	int count;
};
list_bucket *glo_list_wc;
list_bucket *element_edge;

/*input parameter*/
char *name_fiche;
int iteration;
int flag_method;
int left_r;
double alreadybest;
char *type_method;
//int seed;
char *rep;

double deg_rot;
int ialpha;
int ibeta;
int ip;
int igama;
int iran;
double dalpha;

double dp;
double dgama;
double dselect=0.0;
long double *para_f2;
long double *para_f3;


//double *time_local;
//double *time_total;
//char *name_final_result;
//char *name_all_results;
char name_all_results[256];
char name_final_result[256];

/*definition of variable*/
int *seq_node;
int **linked_list;
int **start_end;
int **edge;
int *degree_node;
//int *affected;
int *aff_plus;
int *aff_moins;
//int *BCN;
//int *T;				// Used as an auxiliary vector
//int *O;				// Indicates the neighborhood exploration order
//int **B;			// Best moves found in the neighborhood
//int *n;				// Neighbors of each node of the graph
//int *aux;           // Used to update the weightCount vector efficiently
int *exchange;		// Used to exchange the label in a permutation specially in the Initialisation
int *assit;
int **simple_edge;
int **Bestmove;

int num_read=0;
int limite=0;
int max_vertex1,max_vertex2,num_edge;
int max_degree=0;

int L;
int L_0;
int L_max;
int T;
double P_0;
int **TL;
int num_TL;
int taille_TL;
int gama;
int fp00;

int *CV;

/***************************************/
int **tabu_move;
int list_aj[15]={1,2,1,4,1,2,1,8,1,2,1,4,1,2,1};

/*definition of function*/
void read_fiche(void);
int *get_vector(int size);
int **get_matrix(int num_row,int num_col);
void generate_linked_list(int node1, int node2);
void setdatastructure(void);
void Hillclimbing(void);
void InitialSolution(void);
void GeneratePermutation(void);
void get_fitness(Structure_Solution *s);
void CopySolution(Structure_Solution *s1, Structure_Solution *s2);
int choosebestmove(Structure_Solution *s,int f_p,int ite);
void creatrandompermutation(int *s);
int op_swap_recalculate_local(int node1, int node2,Structure_Solution *C, int *newwc);
void clearmatrix(int *m,int r);
void makemove(int node1, int node2, Structure_Solution *C);
void freedatastructure(void);
void generate_element_edge(int cd, int node1, int node2);
void update_solution(int node1, int node2,Structure_Solution *C, int *newwc);
void export_bucket_list(void);
void op_leftrotation(Structure_Solution *C,int ite);
void freebucket(void);
int check_TL(int v, int u,int ite, int f_p);
void add_TL(int v, int u,int ite, int f_p);
timespec diff(timespec start, timespec end);
int judge_wc(int *newwc, int *otherwc);
void copy_wc(int *fromwc, int *towc);
double calcul_EF3(int CB_max, int *wc);
int judge_one(int *newwc, int *otherwc);
int judge_or(int *newwc, int *otherwc);
void get_close(Structure_Solution *C,int f_p,int ite);
int calcul_tenure(int ite);
int calcul_cyc(int u, int v);

/****************************************************/
/*new_input_to_tuning*/
char *name_instance;		//-i
int seed;					//-seed
int 	part_ls;				//-pls
double	part_th;				//-pth
double	pro_best;			//-prb
double	pro_nf1;			//-prc
double	taux_nf2;			//-tnf
int dep_lr;					//-dep

/****************************/
char resultsFile[100];
char benchmark[100];
/*Program*/




int parameters(int count, char *arguments[])  {
	char *temp, filename[80] = "no";
	char *nf=filename;
//	char *wewant=filename;
//	char *token;
	strcpy(resultsFile, filename);
	strcpy(benchmark, filename);
	strcpy(nf, filename);
	pro_best=-1.0;
	pro_nf1=-1.0;
	taux_nf2=-1.0;
	dep_lr=-1;
	part_ls=-1;
	part_th=-1;

	while (count != 1) {
		temp = arguments[count - 2];
		if (strcmp(temp,"-i") == 0) strcpy(benchmark, arguments[count - 1]);
		else
			if (strcmp(temp,"--seed") == 0) seed= atof(arguments[count - 1]);
		else
			if (strcmp(temp,"-rep") == 0) rep=arguments[count - 1];
		else
			if (strcmp(temp,"-alb") == 0) alreadybest= atoi(arguments[count - 1]);
		else
			if (strcmp(temp,"-pls") == 0) part_ls = atoi(arguments[count - 1]);
		else
			if (strcmp(temp,"-pth") == 0) part_th = atof(arguments[count - 1]);
		else
			if (strcmp(temp,"-prb") == 0) pro_best = atof(arguments[count - 1]);
		else
			if (strcmp(temp,"-prc") == 0) pro_nf1 = atof(arguments[count - 1]);
		else
			if (strcmp(temp,"-tnf") == 0) taux_nf2 = atof(arguments[count - 1]);
		else
			if (strcmp(temp,"-dep") == 0) dep_lr = atoi(arguments[count - 1]);

		else {  // unknow parameter
			return 0;
		}
		count = count - 2;
	}
	if (strcmp(benchmark, "no") == 0 || pro_best == -1 || pro_nf1 == -1 || taux_nf2 == -1 || dep_lr == -1 || part_ls == -1 || part_th == -1) {
		printf("enter error\n");
		exit(-1);
	}
	strcpy(nf, benchmark);

//	token=strtok(nf,"/");
//	while(token!=NULL){
//		wewant=token;
//		token=strtok(NULL,"/");
//	}
//	if(strcmp(wewant,"path100.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"cycle650.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"cycle1000.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"T-dwt__592.mtx.rnd") == 0) alreadybest=7;
//	else if(strcmp(wewant,"cycle200.rnd") == 0) alreadybest=1;
//
//	else if(strcmp(wewant,"caterpillar29.rnd") == 0) alreadybest=24;
//	else if(strcmp(wewant,"hypercube11.rnd") == 0) alreadybest=526;
//	else if(strcmp(wewant,"cycle300.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"path825.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"path200.rnd") == 0) alreadybest=1;
//
//	else if(strcmp(wewant,"X-can__715.mtx.rnd") == 0) alreadybest=52;
//	else if(strcmp(wewant,"mesh2D20x50.rnd") == 0) alreadybest=20;
//	else if(strcmp(wewant,"Q-494_bus.mtx.rnd") == 0) alreadybest=5;
//	else if(strcmp(wewant,"W-685_bus.mtx.rnd") == 0) alreadybest=6;
//	else if(strcmp(wewant,"U-662_bus.mtx.rnd") == 0) alreadybest=5;
//
//	else if(strcmp(wewant,"path650.rnd") == 0) alreadybest=1;
//	else if(strcmp(wewant,"mesh2D5x25.rnd") == 0) alreadybest=5;
//	else if(strcmp(wewant,"tree2x9.rnd") == 0) alreadybest=57;
//	else if(strcmp(wewant,"P-can__445.mtx.rnd") == 0) alreadybest=6;
//	else if(strcmp(wewant,"mesh3D12x12x12.rnd") == 0) alreadybest=114;
//	else alreadybest=1;

//	cout<<alreadybest<<endl;
//	cout<<benchmark<<endl;
//
//	exit(0);
	return 1;
}


void read_fiche(){

	ifstream FIC;
		FIC.open(benchmark);
//	     FIC.open(name_fiche);

	     if ( FIC.fail() ){
	    	 cout << "### No way,check your fiche " << benchmark << endl;
	    	 exit(-1);

	     }
	     char line1[100];

	     FIC.getline(line1,100,'\n'); /*ignore the first line*/
	     FIC>>max_vertex1>>max_vertex2>>num_edge; /*get the number of vertex and edge*/
	     if ( FIC.eof() ){
	    	 cout << "### Your fiche is empty " << benchmark << endl;
	    	 exit(-1);
	     }
//	     cout<<max_vertex1<<max_vertex2<<num_edge<<endl;
	     seq_node=(int*)get_vector(max_vertex1);
	     linked_list=(int**)get_matrix(2*num_edge,2);


	     int i;
	     for (i=0;i<max_vertex1;i++) seq_node[i]=-1;

	     for (i=0;i<2*num_edge;i++) {
	    	 linked_list[i][0]=-1;
	    	 linked_list[i][1]=-1;
	     }

	     int x1,x2,t;
	     int cnt = 0;
	     while(!FIC.eof() && cnt < num_edge){
	    	 FIC>>x1>>x2;
	    	 /*wonder if the node is out of range*/
	    	 if ((x2<1)||(x1<1)||(x1>max_vertex1)||(x2>max_vertex1)) {
	    		 cout<<"the number of node is out of range"<<endl;
	    		 exit(-1);
	    	 }
	    	 /*ensure the x1<x2*/
	    	if(x1!=x2){
	    	 if (x1>x2){
	    		 t=x1;
	    		 x1=x2;
	    		 x2=t;
	    	 }
	    	 	 generate_linked_list(x1,x2);
	    	 	 generate_linked_list(x2,x1);
	    	}
	    	cnt++;
	     }
	     assert(cnt == num_edge);

	     if(num_read/2!=num_edge){
	    	 cout<<"the number of edge is not correct"<<endl;
	    	 exit(-1);
	     }
	     FIC.close();

}
int *get_vector(int size){
	int *Poi;
	Poi= (int*)malloc(size*sizeof(int));
	if (Poi==NULL){
		cout<<"Memory error in get_vector"<<endl;
		exit(-1);
	}
	return Poi;
}

int **get_matrix(int num_row,int num_col){
	int **Poi,i;
	Poi=(int**)malloc(num_row*sizeof(int*));
	if (!Poi){
			cout<<"Memory error in get_matrix"<<endl;
			exit(-1);
		}
	for (i=0;i<num_row;i++){
		Poi[i]=(int*)malloc(num_col*sizeof(int));
		if (!Poi[i]){
					cout<<"Memory error in get_matrix"<<endl;
					exit(-1);
		}
	}
	return Poi;
}

void generate_linked_list(int node1, int node2){
	int j=0,k=0,temp=0;
	int *p_now=NULL;
	j=node1;
	k=node2;
	temp=num_read;
	if (seq_node[j-1]==-1){
		seq_node[j-1]=num_read;
		linked_list[num_read][0]=k-1;
		linked_list[num_read][1]=-1;
		num_read++;
	}
	else{
		p_now=&seq_node[j-1];
		while (*p_now!=-1){
			if (linked_list[*p_now][0]>(k-1)){
				linked_list[num_read][0]=k-1;
				linked_list[num_read][1]=*p_now;
				*p_now=temp;
				num_read++;
				break;
				}
			else{
				p_now=&linked_list[*p_now][1];
			}
		}

		if (*p_now==-1){
				linked_list[num_read][0]=k-1;
				linked_list[num_read][1]=-1;
				*p_now=temp;
				num_read++;
				}
	}

}

void setdatastructure(){
	int line_read=0,i=0,j=0,temp=0;
	int p_now=0;
	limite=max_vertex1/2;

	exchange=(int*)get_vector(max_vertex1);
	assit=(int*)get_vector(max_vertex1);
	Bestmove=(int**)get_matrix(max_vertex1,2);
	start_end=(int**)get_matrix(max_vertex1,2);
	edge=(int**)get_matrix(2*num_edge,2);
	simple_edge=(int**)get_matrix(num_edge,2);
	degree_node=(int*)get_vector(max_vertex1);
	TL=(int**)get_matrix(10*max_vertex1,3);
	tabu_move=(int**)get_matrix(max_vertex1,max_vertex1);
	CV=(int*)get_vector(max_vertex1);

	for(i=0;i<max_vertex1;i++) CV[i]=0;
	for(i=0;i<max_vertex1;i++){
		for(j=0;j<max_vertex1;j++) tabu_move[i][j]=0;
	}

	for (i=0;i<=max_vertex1-1;i++){
		Bestmove[i][0]=0;
		Bestmove[i][1]=1;
	}

	para_f2=(long double*)malloc(max_vertex1*sizeof(long double));
	if (para_f2==NULL){
		cout<<"Memory error in get_vector"<<endl;
		exit(-1);
	}
	memset(para_f2,0,max_vertex1);
	for (i=0;i<=max_vertex1-1;i++){
		para_f2[i]=exp(-(i+1)*10.0/max_vertex1);
//		cout<<i<<" "<<para_f2[i]<<endl;
	}
//	exit(0);

	para_f3=(long double*)malloc((limite+1)*sizeof(long double));
	if (para_f3==NULL){
		cout<<"Memory error in get_vector"<<endl;
		exit(-1);
	}
	memset(para_f3,0,limite+1);
	for (i=1;i<=limite;i++){
		para_f3[i]=1.0/(max_vertex1);
	}
//	para_f3[1]=1.0/(max_vertex1*2);
//	for (i=2;i<=limite;i++){
//		para_f3[i]=para_f3[i-1]/2;
//		cout<<i<<" "<<para_f3[i]<<endl;
//	}

//	BCN=(int*)get_vector(max_vertex1);
//	memset(degree_node,0,max_vertex1);
	CurrentS=(Structure_Solution*)malloc(sizeof(Structure_Solution));
	CurrentS->permutation=(int*)get_vector(max_vertex1*sizeof(int));
	CurrentS->permutationNew=(int*)get_vector(max_vertex1*sizeof(int));
	CurrentS->weightcount=(int*)get_vector((limite+1)*sizeof(int));
	if(CurrentS==NULL){
		cout<<"Memory error in CurrentS"<<endl;
		exit(-1);
	}
	BestS=(Structure_Solution*)malloc(sizeof(Structure_Solution));
	BestS->permutation=(int*)get_vector(max_vertex1*sizeof(int));
	BestS->permutationNew=(int*)get_vector(max_vertex1*sizeof(int));
	BestS->weightcount=(int*)get_vector((limite+1)*sizeof(int));
	if(BestS==NULL){
		cout<<"Memory error in BestS"<<endl;
		exit(-1);
	}
	Thebest=(Structure_Solution*)malloc(sizeof(Structure_Solution));
	Thebest->permutation=(int*)get_vector(max_vertex1*sizeof(int));
	Thebest->permutationNew=(int*)get_vector(max_vertex1*sizeof(int));
	Thebest->weightcount=(int*)get_vector((limite+1)*sizeof(int));
	if(Thebest==NULL){
		cout<<"Memory error in Thebest"<<endl;
		exit(-1);
	}
	LastS=(Structure_Solution*)malloc(sizeof(Structure_Solution));
	LastS->permutation=(int*)get_vector(max_vertex1*sizeof(int));
	LastS->permutationNew=(int*)get_vector(max_vertex1*sizeof(int));
	LastS->weightcount=(int*)get_vector((limite+1)*sizeof(int));
	if(LastS==NULL){
		cout<<"Memory error in LastS"<<endl;
		exit(-1);
	}
	glo_list_wc=(list_bucket*)malloc((limite+1)*sizeof(list_bucket));
	if(glo_list_wc==NULL){
		cout<<"Memory error in glo_list_wc"<<endl;
		exit(-1);
	}
	for(i=0;i<limite+1;i++){
		glo_list_wc[i].flag_empty=1;
		glo_list_wc[i].flag_tail=1;
		glo_list_wc[i].count=0;
//		glo_list_wc[i].next=(list_bucket*)malloc(sizeof(list_bucket));
	}
	for (i=0;i<max_vertex1;i++){
		degree_node[i]=0;
	}

	for(i=0;i<max_vertex1;i++){
		p_now=seq_node[i];
		if (p_now!=-1){
		edge[line_read][0]=i;
		edge[line_read][1]=linked_list[p_now][0];
		start_end[i][0]=line_read;
		line_read++;
		p_now=linked_list[p_now][1];
		while(p_now!=-1){
			edge[line_read][0]=i;
			edge[line_read][1]=linked_list[p_now][0];
			line_read++;
			p_now=linked_list[p_now][1];
		}
		if (p_now==-1){
			start_end[i][1]=line_read-1;
		}
		}
		else {
			cout<<"there is a independent vertex"<<endl;
			exit(-1);
		}
	}
	j=0;
	for (i=0;i<2*num_edge;i++){
		if(edge[i][0]<edge[i][1]){
			simple_edge[j][0]=edge[i][0];
			simple_edge[j][1]=edge[i][1];
			degree_node[edge[i][0]]++;
			degree_node[edge[i][1]]++;
//			cout<<simple_edge[j][0]<<" "<<simple_edge[j][1]<<endl;
			j++;
		}
	}

	for (i=0;i<max_vertex1;i++){
		temp=start_end[i][1]-start_end[i][0]+1;
		if (temp>max_degree){
			max_degree=temp;
		}
	}

	aff_plus=(int*)get_vector(2*max_degree);
	aff_moins=(int*)get_vector(2*max_degree);

	for (i=0;i<num_edge;i++){
		free(linked_list[i]);
		linked_list[i]=NULL;
	}
	free(seq_node);seq_node=NULL;

	for(i=0;i<2*num_edge;i++){
				free(linked_list[i]);
				linked_list[i]=NULL;
			}

	free(linked_list);linked_list=NULL;
}

void Hillclimbing(){
	int flag_move __attribute__((unused));
	double flip __attribute__ ((unused)), endTimeSeconds = 0.0;
	int flag_bruit=0;
	int flag_TL=0;
	int min_TL=0;
	int part;
	int part_TL;
	long long ite=0;
	int i;
	int f_p;
	int f_lr;
//	double cb_last;
	int fp05;
//	long double little=1e-9;
//	long double time_local=0.0;
	long double time_total=0.0;
	long long f_stop=0;
//	timespec time_start;
//	timespec time_end;
	int judge=0;
	int judge_last=0;
//	double method_random;
	int count_left=0;
	double fp05ub=0.0;
//	double fp05lb=0.0;
//	double xlub=0.0;
//	double xllb=0.0;
//	int flag_fp=-1;
//	int ite_bruit;
//	int flag_close=0;
//	clock_t time_start,time_end;			/*time*/

//	ite_bruit=0.1*max_vertex1;

	flag_move=0;
//	min_TL=3*max_vertex1;
	min_TL=part_th*max_vertex1;				/*part of threshold*/
//	min_TL=20;

//	part_TL=max_vertex1*part_ls;			/*part of local search*/
	part_TL=part_ls;
//	part_TL=1/dgama;
//	part_TL=max_vertex1*0.25;
//	cout<<part_TL<<endl;
	part=part_TL+min_TL;


	taille_TL=7*max_vertex1;
	f_lr=1;

    iran=max_vertex1*dgama;
    num_TL=0;
ç

//  if (rot_p==0) rot_p=1;
    if (iran==0) iran=1;

	InitialSolution();
	for (i=0;i<max_vertex1;i++){
		exchange[i]=i;
	}

//	cout<<CurrentS->cbmp<<endl;
/*export the result to file*/
	ofstream caout(name_final_result,ios::out|ios::app);
	if (caout.is_open()){
		caout<<ite<<" ";
    	caout<<CurrentS->cbmp<<" ";
    	caout<<"0"<<endl;
    	caout.close();
   	 }
	/***************************************************/
//	    xllb=0.3;
//	    xlub=0.3;
//	    fp00lb=(1/3-xllb*(CurrentS->cbmp-alreadybest)/(limite-alreadybest))*CurrentS->cbmp;
//	    fp00ub=(0.5-xlub*(CurrentS->cbmp-alreadybest)/(limite-alreadybest))*CurrentS->cbmp;
//	    fp00=fp00ub;
	 fp05ub=(limite-alreadybest)*0.25;
//	 fp05lb=(limite-alreadybest)*0.1;
	 fp05=fp05ub;
	 fp00=1;
	/***************************************************/
	CopySolution(CurrentS,Thebest);
	CopySolution(CurrentS, LastS);
//	cb_last=CurrentS->cbmp+1;
//	fp00=cb_last*0.1;
//	fp00=50;
	f_p=fp00;


	times(&glo_start);
	startTimeSeconds = glo_start.tms_utime/clockTicksPerSecond;

//	while (ite<=iteration && Thebest->cbmp>alreadybest){
	while (Thebest->cbmp>alreadybest){
//		method_random=rand()/(RAND_MAX+1.0);

//		if (method_random>=0.7) flag_method=1;
//		else flag_method=3;
//		if(Thebest->cbmp<=105)
//		{
//			fp00=5;
//			fp05=5;
//		}
//		clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&time_start); /*Start time*/

		ite++;
		f_stop++;
		dselect=rand()/(RAND_MAX+1.0);

/*Point to check*/
		flag_move=choosebestmove(CurrentS,f_p,ite);

		/***********************************************************/
		if (flag_method==1){
			judge=judge_wc(CurrentS->weightcount,Thebest->weightcount);
			judge_last=judge_wc(CurrentS->weightcount,LastS->weightcount);
		}
		else if(flag_method==3){
			judge=judge_or(CurrentS->weightcount,Thebest->weightcount);
			judge_last=judge_or(CurrentS->weightcount,LastS->weightcount);
		}
		else {
			judge=judge_one(CurrentS->weightcount,Thebest->weightcount);
			judge_last=judge_one(CurrentS->weightcount,LastS->weightcount);
		}

		/*store the record*/
		if(CurrentS->cbmp<Thebest->cbmp){
			count_left=0;
			CopySolution(CurrentS, Thebest);
			times(&glo_end);
			endTimeSeconds = glo_end.tms_utime/clockTicksPerSecond;
			time_total= endTimeSeconds - startTimeSeconds;
			/********************************************************************/
			num_TL=0;
//			fp05ub=ceil((Thebest->cbmp-alreadybest)*0.5);
//			fp05lb=ceil((Thebest->cbmp-alreadybest)*0.3);
//			if (fp05ub<=3) fp05ub=3;

			/********************************************************************/
//			cout<<endl;
//			cout<<Thebest->cbmp<<endl;
			cout<<"iteration="<<ite<<" "<<"cbmp="<<Thebest->cbmp<<" time="<<time_total<<" f_p="<<f_p<<endl;
			ofstream caout(name_final_result,ios::out|ios::app);
			if (caout.is_open()){
				caout<<ite<<" ";
				caout<<Thebest->cbmp<<" ";
				caout<<time_total<<" ";
				caout<<" cbmp="<<Thebest->cbmp;
				caout<<endl;
				caout.close();
			}
		}
/*hold on the state of decent*/
		if(judge<0){
//		if(CurrentS->cbmp<Thebest->cbmp){
			CopySolution(CurrentS, Thebest);
			f_p=fp00;
//			num_TL=0;
			flag_bruit=0;
			flag_TL=0;
			f_stop=0;
		}
/*state of LR*/
		else if(flag_bruit>part){
			flag_bruit=0;
			flag_TL=0;
			f_p=fp00;
			f_lr=0;
			num_TL=0;
			if(count_left!=0 && count_left%dep_lr==0)
			{
				for (i=0;i<1;i++) op_leftrotation(CurrentS,ite);
			}
			count_left++;
		}
/*state of TL*/
		else if (((f_lr==1 && judge>=0)||(f_lr==0 && judge_last>=0))){

			flag_bruit++;
			flag_TL++;
			f_lr=1;
			if (flag_TL==part_TL){
				fp05=ceil(Thebest->cbmp*(1/(0.00891104*Thebest->cbmp+0.52663736)+0.16331589));
				f_p=fp05;
			}
		}

/***********************************************************/


		CopySolution(CurrentS, LastS);

		times(&midt);
		endTimeSeconds = midt.tms_utime/clockTicksPerSecond;
		cpuTimeTheBest = endTimeSeconds - startTimeSeconds;
		time_total=cpuTimeTheBest;

		if (time_total>600) {
//			cout<<"stop by the limite of time 600"<<endl;
			break;
		}
		if (f_stop>=5e7){
//			cout<<"time="<<time_total<<" ite="<<ite<<endl;
//			cout<<"stop by the non-improve 5e7"<<endl;
			break;
		}
	}
//	for(i=0;i<max_vertex1;i++) cout<<Thebest->permutation[i]<<" ";
	get_fitness(Thebest);
	cout<<Thebest->cbmp<<endl;
//	printf("%f\n",Thebest->cbmp);


	// ofstream caout(name_final_result,ios::out|ios::app);
	// if (caout.is_open()){
	// 	caout<<Thebest->cbmp;
	// 	caout.close();
	// }

	cout<<"program finished"<<endl;
}
void InitialSolution(){
	int i;
	GeneratePermutation();
	get_fitness(CurrentS);
	for (i=0;i<max_vertex1;i++){
	CurrentS->permutationNew[CurrentS->permutation[i]]=i;
	}
}

void GeneratePermutation(){
	int i,temp=0;
//	cout<<endl;
	for (i=0;i<max_vertex1;i++){
		CurrentS->permutation[i]=i;
//		cout<<CurrentS->permutation[i]<<" ";
	}
//	cout<<endl;

	for (i=0;i<max_vertex1;i++){
		exchange[i]=rand()%max_vertex1;
	}

	for (i=0;i<max_vertex1;i++){
		temp=CurrentS->permutation[i];
		CurrentS->permutation[i]=CurrentS->permutation[exchange[i]];
		CurrentS->permutation[exchange[i]]=temp;
	}
//	for (i=0;i<max_vertex1;i++){
//		cout<<CurrentS->permutation[i]<<" ";
//	}
//	exit(0);

}
void get_fitness(Structure_Solution *s){
	int absoluteDifference,cyclicDifference;
	int i;
	int limite=max_vertex1/2;
	int CB_max;
	s->cbmp=0.0;
	s->EF2=0.0;
	s->EF3=0.0;
	memset(s->weightcount,0,(limite+1)*sizeof(int));
	for (i=0;i<num_edge;i++){
		absoluteDifference=abs(s->permutation[simple_edge[i][0]]-s->permutation[simple_edge[i][1]]);
		cyclicDifference=min(absoluteDifference, max_vertex1-absoluteDifference);
		s->weightcount[cyclicDifference]++;

		generate_element_edge(cyclicDifference,simple_edge[i][0],simple_edge[i][1]);
		glo_list_wc[cyclicDifference].count++;
		if (cyclicDifference>s->cbmp)
			s->cbmp=cyclicDifference;
	}
	/*In fact EF2 is the hashlist*/
	for (i=0;i<=max_vertex1-1;i++){
		s->EF2=s->EF2+i*para_f2[s->permutation[i]];
	}
	/*hash list fini*/
	/*EF3*/
	CB_max=s->cbmp;
//	if (CB_max>=2){
//		s->EF3=CB_max+1.0/((s->weightcount[CB_max-1]+1)*max_vertex1);
//	}
//	else {
//		cout<<"achive the best CB=1"<<endl;
//		exit(0);
//	}
	s->EF3=calcul_EF3(CB_max,s->weightcount);

}

void CopySolution(Structure_Solution *s1, Structure_Solution *s2){
	memcpy(s2->permutation, s1->permutation, max_vertex1*sizeof(int));
	memcpy(s2->permutationNew, s1->permutationNew, max_vertex1*sizeof(int));
	memcpy(s2->weightcount, s1->weightcount, (limite+1)*sizeof(int));
	s2->cbmp = s1->cbmp;
	s2->EF2=s1->EF2;
	s2->EF3=s1->EF3;

}

int choosebestmove(Structure_Solution *s,int f_p,int ite){
	int i=0,j=0,k=0,l=0,h=0;
	int m=0;
//	int n=0;
	int o=0;
	int flag=0;
	int u,v;
	int newCB=0;
	int bestCB;
	int numBest=0;
	int nummove;
	int taille_bucket;
	list_bucket *temp_read;
	int local_vector_read[2*num_edge];
//	int local_vector_read2[2*num_edge];
	int taille_node;
	int flag_arret=0;
//	int flag_continue=0;
	int BC_min;
	int BC_var=0;
	int veri_count=0;
	int ne_nodes[2*num_edge];
	int taille_nei;
	int flag_s;
	int flag_e;
	int dif=0,dif_max=0;
	int node_insert;
	int label_insert;
	int lu,ru,midu;
//	int nodec;
	int set_s[max_vertex1];
	int taille_s;
	int num_incr  __attribute__((unused));
//	int min_incr=10000;
	int temp;
	int flag_tabu;
	static int sel_NF=1;
//	static int count=0;
	static int c_p=0;
//	long double sw_EF2=0.0;
	long double EF3=0.0;
	long double EF3best=0.0;
	long double EF3_wc=0.0;
	int max_CB=0;
	int newwc[limite+1];
//	int bestwc[limite+1];
	int judge_best=0;
	int oldwc[limite+1];
	int judge=0;
	int oldwc_b[limite+1];
	int judge_b=0;
//	int oldwc_w[limite+1];
//	int judge_w=0;
//	int temp_v;
//	int temp_u;
//	int Best_CB;
	int tenure;
	int dis_u;



//	int threshold;
	double p_chooseb;
	int count_neigh=0;
//	int ran_per[max_vertex1];

//	for (i=0;i<max_vertex1;i++) ran_per[i]=i;
	//bestCB=s->cbmp+f_p;
	bestCB=s->cbmp+f_p;
	memset(aff_plus,0,2*max_degree);
	memset(aff_moins,0,2*max_degree);
	taille_bucket=(int)s->cbmp;
	BC_min=taille_bucket*dalpha;
//	Best_CB=Thebest->cbmp;

	max_CB=s->cbmp;
	i=0;
	p_chooseb=exp(-(c_p/50));
	if (p_chooseb<=0.4) p_chooseb=0.4;
	p_chooseb=pro_best;
	EF3best=s->EF3+f_p;
	EF3_wc=s->EF3+f_p;
//	/*bestwc*/
//	for (o=0;o<=limite;o++){
//		bestwc[o]=Thebest->weightcount[o];
//	}

	/*initialise the oldwc: the threshold of NF*/
//	if(max_CB<=limite-f_p-1) {
//		for (o=0;o<=limite;o++){
//			oldwc[o]=glo_list_wc[o].count;
//		}
//		oldwc[max_CB+f_p]=1+glo_list_wc[max_CB+f_p].count;
//	}
//	else {
//		for (o=0;o<=limite;o++){
//			oldwc[o]=glo_list_wc[o].count;
//		}
//	}

	for(i=0;i<max_vertex1;i++) CV[i]=0;

	/*another initialisation*/
	if(max_CB<=limite-f_p-1) {
		for (o=0;o<=limite;o++){
			oldwc[o]=glo_list_wc[o].count;
		}
		oldwc[max_CB+f_p]=max_vertex1+glo_list_wc[max_CB+f_p].count;
	}
	else {
		for (o=0;o<=limite;o++){
			oldwc[o]=glo_list_wc[o].count;
		}
	}
//	/*oldwc2*/
//	threshold=Thebest->cbmp;
//	if(threshold<=limite-f_p-1) {
//
//		for (o=0;o<=limite;o++){
//			oldwc[o]=Thebest->weightcount[i];
//		}
//		oldwc[threshold+f_p]=1+Thebest->weightcount[threshold+f_p];
//	}
//	else {
//		for (o=0;o<=limite;o++){
//			oldwc[o]=Thebest->weightcount[i];
//		}
//	}
	i=0;
	/*browse the neighbor of near-critical nodes*/
	for (BC_var=BC_min;BC_var<=taille_bucket;BC_var++){
		temp_read=&glo_list_wc[BC_var];
		veri_count=temp_read->count+veri_count;
		if (temp_read->count!=0) {
			do{
				temp_read=temp_read->next;
				local_vector_read[i]=temp_read->node[0];

				CV[local_vector_read[i]]=1;					//new adding

				i++;
				local_vector_read[i]=temp_read->node[1];

				CV[local_vector_read[i]]=1;					//new adding
				i++;
				if (i>2*veri_count){
					cout<<"out of range"<<endl;
					exit(-1);
				}

			}while(temp_read->flag_tail==0);
		}

	}
	taille_node=i;

/*browse the neighborhood of the critical vertex*/
	for (i=0;i<taille_node;i++){
//		cout<<"get the Ne1"<<endl;

		if(CV[local_vector_read[i]]==0) continue;
	/*---------avoid the double 3--------*/
//		flag_arret=0;
//		for(j=0;j<i;j++){
//			if (local_vector_read[i]==local_vector_read[j]){
//				flag_arret=1;
//				break;
//			}
//		}
//		if (flag_arret==1) continue;
	/***************************************/
		u = local_vector_read[i];
		CV[u]=0;
	/*Neiborhood 1*/
		if (sel_NF==1){
	/*Create sequence of label of all the neighborhood verties of u*/
			taille_nei=0;
			flag_s=start_end[u][0];
			flag_e=start_end[u][1];
			for (j=flag_s;j<=flag_e;j++){
				node_insert=edge[j][1];
				label_insert=s->permutation[node_insert];
				if (taille_nei==0){
					ne_nodes[taille_nei]=label_insert;
				}
				else {
					for (k=0;k<taille_nei;k++){
						if (k==taille_nei-1){
							if (label_insert>ne_nodes[k]){
								ne_nodes[k+1]=label_insert;
								break;
							}
							else{
								ne_nodes[k+1]=ne_nodes[k];
								ne_nodes[k]=label_insert;
								break;
							}
						}
						else{
							if (label_insert>ne_nodes[k]){
								continue;
							}
							else{
								for (l=taille_nei-1;l>=k;l--){
									ne_nodes[l+1]=ne_nodes[l];
								}
									ne_nodes[k]=label_insert;
									break;
							}
						}
					}
				}
				taille_nei++;
			}
		ne_nodes[taille_nei]=ne_nodes[0]+max_vertex1;
		taille_nei++;
	/*look for the lu and ru, then calcul the midu and then looke for the vertex closer to midu*/
				lu=0;
				ru=0;
				midu=0;
//				nodec=0;
				dif=0;
				dif_max=0;
				for (j=0;j<taille_nei-1;j++){
					dif=ne_nodes[j+1]-ne_nodes[j];
					if (dif>dif_max){
						lu=ne_nodes[j+1];
						ru=ne_nodes[j];
						dif_max=dif;
					}
				}
				midu=(lu+ru+max_vertex1)/2;
				midu=midu%max_vertex1;
				taille_s=0;

				dis_u=calcul_cyc(s->permutation[u], midu);
				for(j=1;j<dis_u;j++){
					if((midu+j+max_vertex1)%max_vertex1!=s->permutation[u]) {
						set_s[taille_s]=s->permutationNew[(midu+j+max_vertex1)%max_vertex1];
						taille_s++;
					}
					if((midu-j+max_vertex1)%max_vertex1!=s->permutation[u]) {
						set_s[taille_s]=s->permutationNew[(midu-j+max_vertex1)%max_vertex1];
						taille_s++;
					}
				}

//				for (j=0;j<taille_nei-1;j++){
//					if(abs(midu-ne_nodes[j])<abs(midu-s->permutation[u])){
//						set_s[taille_s]=s->permutationNew[ne_nodes[j]];
//						taille_s++;
//					}
//				}
					/*--------------browse the voisinage-------------------*/
				for (j=0;j<taille_s;j++){
					v=set_s[j];
//					if (tabu_move[v][u]>=ite  && f_p==fp00) continue; /*tabu move*/
					if (tabu_move[v][u]>=ite) continue;
					if (tabu_move[u][v]>=ite) continue;
					num_incr=op_swap_recalculate_local(v,u,s,newwc);
					count_neigh++;

					for (h=limite;newwc[h]==0;h--);
					newCB=h;
					/*EF3*/
					EF3=calcul_EF3(newCB,newwc);
//					if (newCB>=2){
//						EF3=newCB+1.0/(newwc[newCB-1]+1)*para_f3[newCB-1];
//					}
//					else EF3=0;
					/*Fonction EF3*/
					if (flag_method==0){
						flag_tabu=0;
						if (flag_tabu==0){
							judge=EF3-EF3_wc;
							if(judge<0){
								if (numBest==0){
									EF3best=EF3;
									Bestmove[numBest][0]=v;
									Bestmove[numBest][1]=u;
									numBest++;
									flag=1;
								}
								else {						/*Bestmove is not empty*/
									judge_b=EF3-EF3best; /*judge if better than best neighbor*/
									if (judge_b<=0){					/*yes*/
										if(numBest<max_vertex1-1){
											EF3best=EF3;
											//EF3_wc=EF3;
											Bestmove[numBest][0]=Bestmove[0][0];
											Bestmove[numBest][1]=Bestmove[0][1];
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
											numBest++;
										}
										else{
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
										}

									}
								}
							}
						}
					}
					/*Fonction 1*/
					if (flag_method==1){
						flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
						if (flag_tabu==0){
							judge=judge_wc(newwc,oldwc);
							if(judge<=0){
								if (numBest==0){			/*Bestmove is empty*/
									copy_wc(newwc,oldwc_b);
									Bestmove[numBest][0]=v;
									Bestmove[numBest][1]=u;
									numBest++;
									flag=1;
								}
								else {						/*Bestmove is not empty*/
									judge_b=judge_wc(newwc,oldwc_b); /*judge if better than best neighbor*/
									if (judge_b<=0){					/*yes*/
										if(numBest<max_vertex1-1){
											copy_wc(newwc,oldwc_b);
											//copy_wc(newwc,oldwc); /*try*/
											Bestmove[numBest][0]=Bestmove[0][0];
											Bestmove[numBest][1]=Bestmove[0][1];
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
											numBest++;
										}
										else{
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
										}

									}
									else {
										if(numBest<max_vertex1-1){
											Bestmove[numBest][0]=v;
											Bestmove[numBest][1]=u;
											numBest++;
										}
									}
								}


							}
						}
					}
					/*function 1 finis*/
					/*function 2 start*/
					if(flag_method==2){
						flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
						if (flag_tabu==0){
							judge=judge_one(newwc,oldwc);
							if(judge<=0){
								if (numBest==0){			/*Bestmove is empty*/
									copy_wc(newwc,oldwc_b);
									Bestmove[numBest][0]=v;
									Bestmove[numBest][1]=u;
									numBest++;
									flag=1;
								}
								else{
									judge_b=judge_one(newwc,oldwc_b); /*judge if better than best neighbor*/
									if (judge_b<=0){					/*yes*/
										if(numBest<max_vertex1-1){
											copy_wc(newwc,oldwc_b);
											//copy_wc(newwc,oldwc); /*try*/
											Bestmove[numBest][0]=Bestmove[0][0];
											Bestmove[numBest][1]=Bestmove[0][1];
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
											numBest++;
										}
										else{
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
										}

									}
									else {
										if(numBest<max_vertex1-1){
											Bestmove[numBest][0]=v;
											Bestmove[numBest][1]=u;
											numBest++;
										}
									}
								}
							}
						}
					}
					/*function 2 finis*/
					/*function 3 strat*/
					if(flag_method==3){
						flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
						if (flag_tabu==0){
							//judge=judge_or(newwc,oldwc);
							judge=h-bestCB;
							if(judge<=0){
							//if(1){
								if (numBest==0){			/*Bestmove is empty*/
									copy_wc(newwc,oldwc_b);
									bestCB=h;
									Bestmove[numBest][0]=v;
									Bestmove[numBest][1]=u;
									numBest++;
									flag=1;
								}
								else{
									//judge_b=judge_or(newwc,oldwc_b); /*judge if better than best neighbor*/
									judge_b=h-bestCB;
									if (judge_b<=0){					/*yes*/
										if(numBest<max_vertex1-1){
											copy_wc(newwc,oldwc_b);
											//copy_wc(newwc,oldwc); /*try*/
											bestCB=h;
											Bestmove[numBest][0]=Bestmove[0][0];
											Bestmove[numBest][1]=Bestmove[0][1];
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
											numBest++;
										}
										else{
											Bestmove[0][0]=v;
											Bestmove[0][1]=u;
										}

									}
//									else {
//										if(numBest<max_vertex1-1){
//											Bestmove[numBest][0]=v;
//											Bestmove[numBest][1]=u;
//											numBest++;
//										}
//									}
								}
							}
						}
					}
					/*function 3 finis*/
				}
		}
		/*Neiborhood 2*/
		else if(sel_NF==2){
//			cout<<"get the Ne2"<<endl;
			/*Generate random set for exchanging*/
			for (m=0;m<max_vertex1;m++){
				assit[m]=rand()%max_vertex1;
			}
			for (m=0;m<max_vertex1;m++){
				temp=exchange[m];
				exchange[m]=exchange[assit[m]];
				exchange[assit[m]]=temp;
			}
			/*iran is the maximam number of vertex*/
			for(j=0;j<iran;j++){
				v=exchange[j];
//				if (tabu_move[v][u]>=ite && f_p==fp00) continue;  /*tabu move*/
				if (tabu_move[v][u]>=ite) continue;
				if (tabu_move[u][v]>=ite) continue;
				num_incr=0;
				num_incr=op_swap_recalculate_local(v,u,s,newwc);
				count_neigh++;
				for (h=limite;newwc[h]==0;h--);
				newCB=h;
				/*EF3*/
				EF3=calcul_EF3(newCB,newwc);
//				if (newCB>=2){
//					EF3=newCB+1.0/(newwc[newCB-1]+1)*para_f3[newCB-1];
//				}
//				else EF3=0;
				/*Fonction EF3*/
				if (flag_method==0){
					flag_tabu=0;
					if (flag_tabu==0){
						judge=EF3-EF3_wc;
						if(judge<0){
							if (numBest==0){
								EF3best=EF3;
								Bestmove[numBest][0]=v;
								Bestmove[numBest][1]=u;
								numBest++;
								flag=1;
							}
							else {						/*Bestmove is not empty*/
								judge_b=EF3-EF3best; /*judge if better than best neighbor*/
								if (judge_b<=0){					/*yes*/
									if(numBest<max_vertex1-1){
										EF3best=EF3;
										//EF3_wc=EF3;
										Bestmove[numBest][0]=Bestmove[0][0];
										Bestmove[numBest][1]=Bestmove[0][1];
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
										numBest++;
									}
									else{
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
									}

								}
							}
						}
					}
				}
				/*Fonction 1*/
				if (flag_method==1){
					flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
					if (flag_tabu==0){
						judge=judge_wc(newwc,oldwc);
						if(judge<=0){
							if (numBest==0){			/*Bestmove is empty*/
								copy_wc(newwc,oldwc_b);
								Bestmove[numBest][0]=v;
								Bestmove[numBest][1]=u;
								numBest++;
								flag=1;
							}
							else {						/*Bestmove is not empty*/
								judge_b=judge_wc(newwc,oldwc_b); /*judge if better than best neighbor*/
								if (judge_b<=0){					/*yes*/
									if(numBest<max_vertex1-1){
										copy_wc(newwc,oldwc_b);
										//copy_wc(newwc,oldwc); /*try*/
										bestCB=h;
										Bestmove[numBest][0]=Bestmove[0][0];
										Bestmove[numBest][1]=Bestmove[0][1];
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
										numBest++;
									}
									else{
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
									}

								}
								else {
									if(numBest<max_vertex1-1){
										Bestmove[numBest][0]=v;
										Bestmove[numBest][1]=u;
										numBest++;
									}
								}
							}


						}
					}
				}
				/*function 1 finis*/
				/*function 2 start*/
				if(flag_method==2){
					flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
					if (flag_tabu==0){
						judge=judge_one(newwc,oldwc);
						if(judge<=0){
							if (numBest==0){			/*Bestmove is empty*/
								copy_wc(newwc,oldwc_b);
								Bestmove[numBest][0]=v;
								Bestmove[numBest][1]=u;
								numBest++;
								flag=1;
							}
							else{
								judge_b=judge_one(newwc,oldwc_b); /*judge if better than best neighbor*/
								if (judge_b<=0){					/*yes*/
									if(numBest<max_vertex1-1){
										copy_wc(newwc,oldwc_b);
										//copy_wc(newwc,oldwc); /*try*/
										Bestmove[numBest][0]=Bestmove[0][0];
										Bestmove[numBest][1]=Bestmove[0][1];
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
										numBest++;
									}
									else{
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
									}

								}
								else {
									if(numBest<max_vertex1-1){
										Bestmove[numBest][0]=v;
										Bestmove[numBest][1]=u;
										numBest++;
									}
								}
							}
						}
					}
				}
				/*function 2 finis*/
				/*function 3 strat*/
				if(flag_method==3){
					flag_tabu=0; /*check if in the state of tabu and in th TL if in the state of tabu*/
					if (flag_tabu==0){
						//judge=judge_or(newwc,oldwc);
						judge=h-bestCB;
						//judge_b=h-bestCB;
						if(judge<=0){
						//if(1){
							if (numBest==0){			/*Bestmove is empty*/
								copy_wc(newwc,oldwc_b);
								Bestmove[numBest][0]=v;
								Bestmove[numBest][1]=u;
								numBest++;
								flag=1;
							}
							else{
								//judge_b=judge_or(newwc,oldwc_b); /*judge if better than best neighbor*/
								judge_b=h-bestCB;
								if (judge_b<=0){					/*yes*/
									if(numBest<max_vertex1-1){
										copy_wc(newwc,oldwc_b);
										//copy_wc(newwc,oldwc); /*try*/
										bestCB=h;
										Bestmove[numBest][0]=Bestmove[0][0];
										Bestmove[numBest][1]=Bestmove[0][1];
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
										numBest++;
									}
									else{
										Bestmove[0][0]=v;
										Bestmove[0][1]=u;
									}

								}
//								else {
//									if(numBest<max_vertex1-1){
//										Bestmove[numBest][0]=v;
//										Bestmove[numBest][1]=u;
//										numBest++;
//									}
//								}
							}
						}
					}
				}
				/*function 3 finis*/

			}
		}
	}

			/*NF2 finis*/
			/*NF3 start*/
			/*NF3 finis*/


//	cout<<"finish the opswap"<<endl;

	/*make move*/

	if (flag){
		// if(f_p==fp00) nummove=0;
		// else{
			if(p_chooseb<=(rand()/(RAND_MAX+1.0))) nummove=0;
			else nummove=rand()%numBest;
		// }

		tenure=calcul_tenure(ite);
//		if (f_p==fp00){
		if(1){
//			tabu_move[Bestmove[nummove][0]][Bestmove[nummove][1]]=ite+gama;
//			tabu_move[Bestmove[nummove][1]][Bestmove[nummove][0]]=ite+gama;
			tabu_move[Bestmove[nummove][0]][Bestmove[nummove][1]]=ite+tenure;
			tabu_move[Bestmove[nummove][1]][Bestmove[nummove][0]]=ite+tenure;
		}

		makemove(Bestmove[nummove][0], Bestmove[nummove][1], s);
	}
	if (flag==0) nummove=0;

	/*judge_best*/
if(flag_method==1){
	judge_best=judge_wc(s->weightcount,Thebest->weightcount);
}
else if(flag_method==3){
	judge_best=judge_or(s->weightcount,Thebest->weightcount);
}

else{
	judge_best=judge_one(s->weightcount,Thebest->weightcount);
}
 	/*sortir all the results*/
// 		ofstream allout(name_all_results,ios::out|ios::app);
// 		if (allout.is_open()){
// 			allout<<ite<<'\t'<<"NB="<<numBest<<'\t'<<"c_n="<<count_neigh<<'\t'<<"f_p="<<f_p<<'\t'<<"f="<<flag<<"  NF="<<sel_NF;
// 			if (judge_best<0) allout<<" "<<"judge=1";
// 			else allout<<" "<<"judge=0";
// 			allout<<"  "<<count;
// 			allout<<"  taille_node="<<taille_node;
// //			for (o=limite;oldwc[o]==0 && o>=0;o--);
// 			allout<<" cbmp="<<s->cbmp<<" CB=";
// 			max_CB=s->cbmp;
// 			for (m=max_CB;m>=max_CB-10;m--) {
// 				if(m>=0){
// 				allout<<s->weightcount[m]<<" ";
// 				}
// 				else break;
// 			}
// 			allout<<setprecision(12)<<"EF2="<<s->EF2;
// 			allout<<" nummove="<<nummove<<" move="<<Bestmove[nummove][0]<<" "<<Bestmove[nummove][1]<<" "<<num_TL;
//
// 			allout<<endl;
// 	    	allout.close();
// 	   	 }
 	/***************************************/
	/*select the next iteration NF in mode LS*/
		if(f_p==fp00){
			if(numBest==0) sel_NF=2;
			else sel_NF=1;
//			if (judge_best<0){
//				sel_NF=1;
//				count=0;
//			}
//			else if(count==100 || flag==0 || numBest<=1){
//				sel_NF=2;
//			}
//			else if (sel_NF==2){
//				sel_NF=1;
//				count=0;
//			}
//			else count++;
		}
	/*select the next iteration NF in mode BLS*/
		else{
			if (dselect>=dp) sel_NF=1;
			else  sel_NF=2;
	//		else sel_NF=3;
			if (flag==0) sel_NF=2;
		}



	/*gerer c_p*/
		if (judge_best<0) c_p=0;
		else {
			c_p++;
			// tabu_move[Bestmove[nummove][0]][Bestmove[nummove][1]]=ite+tenure;
			// tabu_move[Bestmove[nummove][1]][Bestmove[nummove][0]]=ite+tenure;
		}

	return flag;
}

void creatrandompermutation(int *s){
	int *temp;
	int i,exc;
	temp=(int*)get_vector(max_vertex1);
	for (i=0;i<max_vertex1;i++){
		exchange[i]=i;
	}
	for (i=0;i<max_vertex1;i++){
		temp[i]=rand()%max_vertex1;
	}
	for (i=0;i<max_vertex1;i++){
		exc=exchange[i];
		exchange[i]=exchange[temp[i]];
		exchange[temp[i]]=exc;
	}

}

int op_swap_recalculate_local(int node1,int node2,Structure_Solution *C,int *newwc){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i,j=0,k=0;
//	int w;
//	int m;
	int list_start,list_end;
	int abd,cycd;
//	long double sw_EF2;

//	sw_EF2=C->EF2;
	clearmatrix(aff_plus,2*max_degree);
	clearmatrix(aff_moins,2*max_degree);

	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];

	for (i=0;i<limite+1;i++){
		newwc[i]=C->weightcount[i];
	}


	list_start=start_end[node1][0];
	list_end=start_end[node1][1];
	for (i=list_start;i<list_end+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){
			/*get the new cycd*/
			abd=abs(lab_oldnode1-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]--;
			aff_moins[j]=cycd; j++;


			abd=abs(lab_newnode1-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]++;
			aff_plus[k]=cycd; k++;
		}
	}

	list_start=start_end[node2][0];
	list_end=start_end[node2][1];
	for (i=list_start;i<list_end+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){

			/*calcul the cycd*/
			abd=abs(lab_oldnode2-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]--;
			aff_moins[j]=cycd; j++;

			abd=abs(lab_newnode2-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]++;
			aff_plus[k]=cycd; k++;
		}
	}
	return 0;
}

void clearmatrix(int *m, int r){
	int i;
	for(i=0;i<r;i++){
		m[i]=0;
	}
}

void makemove(int node1, int node2, Structure_Solution *C){
	int temp_node;
	int i;
	int newwc[limite+1];

	update_solution(node1, node2, C, newwc);

	for (i=0;i<limite+1;i++){
		C->weightcount[i]=newwc[i];
	}
	for (i=limite;newwc[i]==0;i--);
	C->cbmp=i;
	C->EF3=calcul_EF3(i,newwc);
//	if (i>=2){
//		C->EF3=i+1.0/(newwc[i-1]+1)*para_f3[i-1];
//	}
//	else C->EF3=0;

	temp_node=C->permutation[node1];
	C->permutation[node1]=C->permutation[node2];
	C->permutation[node2]=temp_node;

	C->permutationNew[C->permutation[node1]]=node1;
	C->permutationNew[C->permutation[node2]]=node2;

//	cout<<"move finished"<<endl;
//	cout<<endl;
}

void freedatastructure(){
	int i;
	list_bucket *temp_read,*last_read;

	for(i=0;i<2*num_edge;i++){
		free(edge[i]);
		edge[i]=NULL;
	}
	for(i=0;i<num_edge;i++){
			free(simple_edge[i]);
			simple_edge[i]=NULL;
		}
	for(i=0;i<max_vertex1;i++){
				free(Bestmove[i]);
				Bestmove[i]=NULL;
			}
	for(i=0;i<max_vertex1;i++){
				free(start_end[i]);
				start_end[i]=NULL;
			}
	for(i=0;i<max_vertex1;i++){
				free(TL[i]);
				TL[i]=NULL;
			}
	for(i=0;i<max_vertex1;i++){
				free(tabu_move[i]);
				tabu_move[i]=NULL;
			}
	free(start_end); 					start_end=NULL;
	free(TL);							TL=NULL;
	free(tabu_move);					tabu_move=NULL;
	free(edge); 						edge=NULL;
	free(simple_edge); 					simple_edge=NULL;
	free(Bestmove); 					Bestmove=NULL;
	free(aff_plus); 					aff_plus=NULL;
	free(aff_moins); 					aff_moins=NULL;
	free(exchange); 					exchange=NULL;
	free(assit); 						assit=NULL;
	free(degree_node);					degree_node=NULL;
	free(para_f2);						para_f2=NULL;

	free(CurrentS->permutation); 		CurrentS->permutation=NULL;
	free(CurrentS->permutationNew); 	CurrentS->permutationNew=NULL;
	free(CurrentS->weightcount); 		CurrentS->weightcount=NULL;
	free(CurrentS);						CurrentS=NULL;

	free(BestS->permutation); 			BestS->permutation=NULL;
	free(BestS->permutationNew); 		BestS->permutationNew=NULL;
	free(BestS->weightcount); 			BestS->weightcount=NULL;
	free(BestS);						BestS=NULL;

	free(Thebest->permutation); 		Thebest->permutation=NULL;
	free(Thebest->permutationNew); 		Thebest->permutationNew=NULL;
	free(Thebest->weightcount); 		Thebest->weightcount=NULL;
	free(Thebest);						Thebest=NULL;

	free(LastS->permutation); 			LastS->permutation=NULL;
	free(LastS->permutationNew); 		LastS->permutationNew=NULL;
	free(LastS->weightcount); 			LastS->weightcount=NULL;
	free(LastS);						LastS=NULL;
	free(CV);							CV=NULL;

	for(i=0;i<limite+1;i++){
		temp_read=&glo_list_wc[i];
		if (temp_read->flag_tail==0){
			if(temp_read->flag_tail==0){
				last_read=temp_read;
				temp_read=temp_read->next;
			}
			while(temp_read->flag_tail==0){
					last_read=temp_read;
					temp_read=temp_read->next;
					free(last_read);
					last_read=NULL;;
				}
			free(temp_read);
			temp_read=NULL;
		}

	}
	free(glo_list_wc);
	glo_list_wc=NULL;
}

void generate_element_edge(int cd, int node1, int node2){

	element_edge=(list_bucket*)malloc(sizeof(list_bucket));
	if (element_edge==NULL){
		cout<<"error memory in element_edge"<<endl;
		exit(-1);
	}
//	element_edge->next=(list_bucket*)malloc(sizeof(list_bucket));

	element_edge->node[0]=node1;
	element_edge->node[1]=node2;
	element_edge->flag_empty=0;

	if (glo_list_wc[cd].flag_empty==1){
		glo_list_wc[cd].flag_empty=0;
		glo_list_wc[cd].flag_tail=0;
		glo_list_wc[cd].next=element_edge;
		element_edge->flag_tail=1;
	}
	else if(glo_list_wc[cd].flag_empty!=1){
		element_edge->flag_tail=0;
		element_edge->next=glo_list_wc[cd].next;
		glo_list_wc[cd].next=element_edge;
	}
}

void update_solution(int node1, int node2,Structure_Solution *C, int *newwc){
	int lab_oldnode1,lab_oldnode2;
	int lab_newnode1,lab_newnode2;
	int i=0,j=0,k=0;
	int list_start,list_end;
	int abd,cycd;
	int v1,v2,tempv;

	list_bucket *last_read;
	list_bucket *temp_read;
	list_bucket *element_move=NULL;
//	int local_flag_tail;
	int position=0;

//	cout<<"arrive at the update solution"<<endl;

	clearmatrix(aff_plus,2*max_degree);
	clearmatrix(aff_moins,2*max_degree);
	C->EF2=C->EF2-node1*para_f2[C->permutation[node1]]-node2*para_f2[C->permutation[node2]]+node1*para_f2[C->permutation[node2]]+node2*para_f2[C->permutation[node1]];

	lab_newnode2=lab_oldnode1=C->permutation[node1];
	lab_newnode1=lab_oldnode2=C->permutation[node2];

	for (i=0;i<limite+1;i++){
		newwc[i]=C->weightcount[i];
	}

	list_start=start_end[node1][0];
	list_end=start_end[node1][1];
	for (i=list_start;i<list_end+1;i++){
		if (lab_newnode1!=C->permutation[edge[i][1]]){

			v1=node1;
			v2=edge[i][1];
			if (v1>v2){				//ensure the v1<v2
				tempv=v2;
				v2=v1;
				v1=tempv;
			}
			/*calcul the old wc*/
			abd=abs(lab_oldnode1-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]--;
			aff_moins[j]=cycd; j++;

			/*delete the element_edge from the origin row*/
//			export_bucket_list();
			glo_list_wc[cycd].count--;
//			cout<<"cycd-= "<<cycd<<endl;
//			cout<<"v1= "<<v1<<" v2= "<<v2<<endl;
			temp_read=&glo_list_wc[cycd];
			position=0;

			do{
				last_read=temp_read;
				temp_read=temp_read->next;
				position++;
//				cout<<"position= "<<position<<endl;
				if(position>newwc[cycd]+1){
					cout<<"---------------------------------"<<endl;
					cout<<"position out of range1"<<endl;
					cout<<"---------------------------------"<<endl;
					exit(-1);
				}
				if (temp_read->node[0]==v1 && temp_read->node[1]==v2){
					element_move=temp_read;
					break;
				}

			}while(temp_read->flag_tail==0);

			if (position==1 && element_move->flag_tail==1){ 		//only one (single)
				glo_list_wc[cycd].flag_empty=1;
				glo_list_wc[cycd].flag_tail=1;
			}
			if (position!=1 && element_move->flag_tail==1){ 		//not first one but the last one
				last_read->flag_tail=1;
			}
			if (position!=1 && element_move->flag_tail!=1){ 		//not first one and not last one
				last_read->next=element_move->next;
			}
			if (position==1 && element_move->flag_tail!=1){ 		//first one and not the last one
				last_read->next=element_move->next;
			}
//			export_bucket_list();
			/*calcul the new wc*/
			abd=abs(lab_newnode1-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]++;
			aff_plus[k]=cycd; k++;

			glo_list_wc[cycd].count++;
//			cout<<"cycd+= "<<cycd<<endl;
			/*add the element_move to the new place*/

			if (glo_list_wc[cycd].flag_empty==1){
				glo_list_wc[cycd].flag_empty=0;
				glo_list_wc[cycd].flag_tail=0;
				glo_list_wc[cycd].next=element_move;
				element_move->flag_tail=1;
				}
			else if(glo_list_wc[cycd].flag_empty!=1){
				element_move->flag_tail=0;
				element_move->next=glo_list_wc[cycd].next;
				glo_list_wc[cycd].next=element_move;
				}
//			export_bucket_list();
		}
//		export_bucket_list();
	}
	list_start=start_end[node2][0];
	list_end=start_end[node2][1];
	for (i=list_start;i<list_end+1;i++){
		if (lab_newnode2!=C->permutation[edge[i][1]]){

			v1=node2;
			v2=edge[i][1];
			if (v1>v2){				//ensure the v1<v2
				tempv=v2;
				v2=v1;
				v1=tempv;
			}
			/*calcul the old wc*/
			abd=abs(lab_oldnode2-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]--;
			aff_moins[j]=cycd; j++;
			glo_list_wc[cycd].count--;

			position=0;
			temp_read=&glo_list_wc[cycd];

			do{
				last_read=temp_read;
				temp_read=temp_read->next;
				position++;
//				cout<<"position= "<<position<<endl;
				if(position>newwc[cycd]+1){
					cout<<"---------------------------------"<<endl;
					cout<<"position out of range2"<<endl;
					cout<<"---------------------------------"<<endl;
					exit(-1);
				}
				if (temp_read->node[0]==v1 && temp_read->node[1]==v2){
					element_move=temp_read;
					break;
				}

			}while(temp_read->flag_tail==0);

			if (position==1 && element_move->flag_tail==1){ 		//only one (single)
				glo_list_wc[cycd].flag_empty=1;
				glo_list_wc[cycd].flag_tail=1;
			}
			if (position!=1 && element_move->flag_tail==1){ 		//not first one but the last one
				last_read->flag_tail=1;
			}
			if (position!=1 && element_move->flag_tail!=1){ 		//not first one and not last one
				last_read->next=element_move->next;
			}
			if (position==1 && element_move->flag_tail!=1){ 		//first one and not the last one
				last_read->next=element_move->next;
			}
			/*calcul the new wc*/
			abd=abs(lab_newnode2-C->permutation[edge[i][1]]);
			cycd=min(abd, max_vertex1-abd); newwc[cycd]++;
			aff_plus[k]=cycd; k++;

			glo_list_wc[cycd].count++;
//			cout<<"cycd+= "<<cycd<<endl;
			/*add the element_move*/
			if (glo_list_wc[cycd].flag_empty==1){
				glo_list_wc[cycd].flag_empty=0;
				glo_list_wc[cycd].flag_tail=0;
				glo_list_wc[cycd].next=element_move;
				element_move->flag_tail=1;
				}
			else if(glo_list_wc[cycd].flag_empty!=1){
				element_move->flag_tail=0;
				element_move->next=glo_list_wc[cycd].next;
				glo_list_wc[cycd].next=element_move;
				}

//			export_bucket_list();
		}

	}
//	cout<<"finished the update"<<endl;

}

void export_bucket_list(){
	int i;
	list_bucket *temp_read;
	/*-----------------------------------------------------------*/
	cout<<"-------------bucket-------------"<<endl;
	for(i=limite;i>0;i--){
		cout<<"wc="<<i<<" "<<"#"<<glo_list_wc[i].count<<" ";
		temp_read=&glo_list_wc[i];
		if (temp_read->flag_empty==1){
			cout<<endl;
			continue;
		}
		do{
			temp_read=temp_read->next;
			cout<<temp_read->node[0]<<"/"<<temp_read->node[1]<<" ";

		}while(temp_read->flag_tail==0);
		cout<<endl;
	}
	cout<<"----------------------------------------------------"<<endl;
	cout<<endl;
	/*--------------------------------------------------------*/

}
void op_leftrotation(Structure_Solution *C,int ite){
	int temp_vertex;
	int times_shift;
	int i;
	int big,small;
	int temp_u;
	int temp_v;

	int taille_edge;
	int BC_min;
	int BC_var=0;
	int taille_bucket;
	list_bucket *temp_read;
	int veri_count=0;
	int local_vector_read[2*num_edge];
	int p_chosen;
	double psel;

	taille_bucket=(int)C->cbmp;

//	times_shift=rand()%taille_bucket+1;
	times_shift=taille_bucket;

	BC_min=taille_bucket*dalpha;
	i=0;
	/*browse the neighbor of near-critical nodes*/
	for (BC_var=BC_min;BC_var<=taille_bucket;BC_var++){
		temp_read=&glo_list_wc[BC_var];
		veri_count=temp_read->count+veri_count;
		if (temp_read->count!=0) {
			do{
						temp_read=temp_read->next;
						local_vector_read[i]=temp_read->node[0];
						i++;
						local_vector_read[i]=temp_read->node[1];
						i++;
						if (i>2*veri_count){
							cout<<"out of range in leftrotation"<<endl;
							exit(-1);
						}

				}while(temp_read->flag_tail==0);
		}

	}
	taille_edge=i/2;
	p_chosen=rand()%taille_edge;
	psel=rand()/(rand()/(RAND_MAX+1.0));
	p_chosen=p_chosen*2;
	big=local_vector_read[p_chosen];
	small=local_vector_read[p_chosen+1];
	if (C->permutation[big]<C->permutation[small]) {
		temp_vertex=small;
		small=big;
		small=temp_vertex;
	}
	/*****************************************/
	if(abs(C->permutation[big]-C->permutation[small])<max_vertex1/2){
		if(psel>0.5){
			temp_u=big;
			for(i=0;i<times_shift;i++){
				temp_v=C->permutationNew[C->permutation[temp_u]-1];
				makemove(temp_u,temp_v,C);
				// tabu_move[temp_u][temp_v]=ite+calcul_tenure(ite);
				// tabu_move[temp_v][temp_u]=ite+calcul_tenure(ite);
//				temp_u=temp_v;
			}
		}
		else {
			temp_u=small;
			for(i=0;i<times_shift;i++){
				temp_v=C->permutationNew[C->permutation[temp_u]+1];
				makemove(temp_u,temp_v,C);
				// tabu_move[temp_u][temp_v]=ite+calcul_tenure(ite);
				// tabu_move[temp_v][temp_u]=ite+calcul_tenure(ite);
//				temp_u=temp_v;
			}
		}
	}
	else{
		if(psel>0.5){
			temp_u=big;
			for (i=0;i<times_shift;i++){
				temp_v=C->permutationNew[(C->permutation[temp_u]+1+max_vertex1)%max_vertex1];
				makemove(temp_u,temp_v,C);
				// tabu_move[temp_u][temp_v]=ite+calcul_tenure(ite);
				// tabu_move[temp_v][temp_u]=ite+calcul_tenure(ite);				
//				temp_v=temp_u;
			}
		}
		else{
			temp_u=small;
			for (i=0;i<times_shift;i++){
				temp_v=C->permutationNew[(C->permutation[temp_u]-1+max_vertex1)%max_vertex1];
				makemove(temp_u,temp_v,C);
				// tabu_move[temp_u][temp_v]=ite+calcul_tenure(ite);
				// tabu_move[temp_v][temp_u]=ite+calcul_tenure(ite);				
//				temp_v=temp_u;
			}
		}
	}



}
void freebucket(){

}
int check_TL(int v, int u,int ite,int f_p){
	int i;
	int x1,x2;
	if (v<u) {
		x1=v;
		x2=u;
	}
	else {
		x1=u;
		x2=v;
	}
//	if (num_TL!=0){
	if (f_p==fp00 && num_TL!=0){
		for (i=0;i<num_TL;i++){
			if ((x1==TL[i][0]) && (x2==TL[i][1]) && (ite<=TL[i][2])) return 1;
//			else return 0;
		}
	}
	return 0;
}

void add_TL(int v, int u,int ite,int f_p){
	int x1,x2;
	if (v<u) {
		x1=v;
		x2=u;
	}
	else {
		x1=u;
		x2=v;
	}
//	if (num_TL<taille_TL){
	if (num_TL<taille_TL && f_p==fp00){
		TL[num_TL][0]=x1;
		TL[num_TL][1]=x2;
		TL[num_TL][2]=ite+gama;
		num_TL++;
	}
//	if (num_TL==taille_TL){
	if (num_TL==taille_TL && f_p==fp00){
		num_TL=0;
		TL[num_TL][0]=x1;
		TL[num_TL][1]=x2;
		TL[num_TL][2]=ite+gama;
	}

//	if (num_TL==taille_TL-1){
//		cout<<"num_TL not enough"<<endl;
//	}
/*	if (f_p!=fp00){
		ofstream allout(name_all_results,ios::out|ios::app);
				if (allout.is_open()){
					allout<<ite<<'\t'<<TL[num_TL-1][0];
					allout<<'\t'<<TL[num_TL-1][1];
					allout<<'\t'<<TL[num_TL-1][2];
					allout<<endl;
			    	allout.close();
			   	 }
	}*/

}

timespec diff(timespec start, timespec end){
	 timespec temp;
	 if ((end.tv_nsec-start.tv_nsec)<0) {
	 temp.tv_sec = end.tv_sec-start.tv_sec-1;
	 temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	 } else {
	 temp.tv_sec = end.tv_sec-start.tv_sec;
	 temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	 }
	 return temp;
}

int judge_wc(int *newwc, int *otherwc){
	int o,m;
	int deltawc[limite+1];
	for (o=0;o<=limite;o++){
		deltawc[o]=newwc[o]-otherwc[o];
	}
	for (m=limite;deltawc[m]==0 && m>=0;m--);
	return deltawc[m];
}
void copy_wc(int *fromwc, int *towc){
	int o;
	for (o=0;o<=limite;o++){
		towc[o]=fromwc[o];
	}
}
double calcul_EF3(int CB_max, int *wc){
	double EF3;
//	int bestcb;
//	bestcb=Thebest->cbmp;
	if (CB_max>=2){
		EF3=CB_max+1.0/((wc[CB_max-1]+1)*max_vertex1);
	}
	else {
		EF3=0;
	}
	return EF3;
}
int judge_one(int *newwc, int *otherwc){
	int o,m;
	int h;
	int key_c;
	int deltawc[limite+1];
	for (o=0;o<=limite;o++){
		deltawc[o]=newwc[o]-otherwc[o];
	}
	for (h=limite;otherwc[h]==0;h--);
	key_c=h;
	for (m=limite;deltawc[m]==0 && m>=key_c;m--);
	return deltawc[m];
}

int judge_or(int *newwc, int *otherwc){
	int o;
	int cb_new;
	int cb_old;

	for(o=limite;newwc[o]==0;o--);
	cb_new=o;

	for(o=limite;otherwc[o]==0;o--);
	cb_old=o;

	return (cb_new-cb_old);
}

int calcul_tenure(int ite){
	int aj;
	int tenure;
	int index_j;
	index_j=floor((ite%1500)/100);
	aj=list_aj[index_j];
	tenure=(aj-1)*100+rand()%100;
	return tenure;
}

void get_close(Structure_Solution *C,int f_p,int ite){
	int i;
	int taille_node;
	int BC_min;
	int BC_var=0;
	int taille_bucket;
	list_bucket *temp_read;
	int veri_count=0;
	int local_vector_read[2*num_edge];
	int u;
	int num;
	int list_nei[max_vertex1];
	int list_ex[max_vertex1];
	int u_deg;
	int deg05;
	int start;
	int end;
	int j,k,m;
	int label_u;
	int flag_tabu;

	taille_bucket=(int)C->cbmp;
	BC_min=taille_bucket*dalpha;
	i=0;
	/*browse the neighbor of near-critical nodes*/
	for (BC_var=BC_min;BC_var<=taille_bucket;BC_var++){
		temp_read=&glo_list_wc[BC_var];
		veri_count=temp_read->count+veri_count;
		if (temp_read->count!=0) {
			do{
						temp_read=temp_read->next;
						local_vector_read[i]=temp_read->node[0];
						i++;
						local_vector_read[i]=temp_read->node[1];
						i++;
						if (i>2*veri_count){
							cout<<"out of range in getclose"<<endl;
							exit(-1);
						}

				}while(temp_read->flag_tail==0);
		}

	}
	taille_node=i;
	num=rand()%taille_node;
	u=local_vector_read[num];
	label_u=C->permutation[u];

	start=start_end[u][0];
	end=start_end[u][1];
	u_deg=end-start+1;
	deg05=ceil(u_deg/2);

	j=0;
	k=0;
	for(i=start;i<=end;i++){
		list_nei[j]=edge[i][1];
		j++;
	}

	for(m=1;m<=deg05;m++){
		list_ex[k]=C->permutationNew[(label_u+m+max_vertex1)%max_vertex1];
		k++;
		list_ex[k]=C->permutationNew[(label_u-m+max_vertex1)%max_vertex1];
		k++;
	}
	for(i=0;i<u_deg;i++){
		flag_tabu=check_TL(list_nei[i],list_ex[i],ite,f_p);
		if(flag_tabu==0){
			makemove(list_nei[i],list_ex[i],C);
			add_TL(list_nei[i],list_ex[i],ite,f_p);
		}

	}

}

int calcul_cyc(int u, int v){
	int abd,cycd;
	abd=abs(u-v);
	cycd=min(abd, max_vertex1-abd);//C->weightcount[cycd]--;
	return cycd;
}

int main(int argc, char **argv){


	char final[50]="F";
//	char all[50]="A";
	int n_02;
	int n_09;
//	int n_05;
	int n_04;
	char tm[100]="2";

	type_method=tm;
	if(parameters(argc, argv)==0){
		cout<<benchmark<<" "<<seed<<" "<<part_ls<<" "<<part_th<<" "<<pro_best<<" "<<pro_nf1<<" "<<taux_nf2<<" "<<dep_lr<<endl;
		exit(-1);
	}
	cout<<benchmark<<" "<<seed<<" "<<part_ls<<" "<<part_th<<" "<<pro_best<<" "<<pro_nf1<<" "<<taux_nf2<<" "<<dep_lr<<endl;
//	else{
//		cout<<benchmark<<" "<<seed<<" "<<part_ls<<" "<<part_th<<" "<<pro_best<<" "<<pro_nf1<<" "<<taux_nf2<<" "<<dep_lr<<endl;
//		exit(-1);
//	}
//    if (argc == 12)
//    {
//          name_fiche = argv[1];
//          iteration = atoi(argv[2]);
//          flag_method=atoi(argv[3]);
//          type_method=argv[3];
//          seed=atoi(argv[4]);
//          rep=argv[5];
//          ialpha=atoi(argv[6]);
//          ibeta=atoi(argv[7]);
//          ip=atoi(argv[8]);
//          igama=atoi(argv[9]);
//          left_r=atoi(argv[10]);
//          alreadybest=atof(argv[11]);
//    }
//	if (argc == 5)
//     {
//           name_fiche = argv[1];
//           iteration = atoi(argv[2]);
//           name_final_result = argv[3];
//           name_all_results = argv[4];
//     }
//     else {
//           cout << endl << "Syntaxes :" << endl;
//           cout << "LOL, number of input is incorrect" << endl;
//           cout << "need only the name of map, iteration, name of final result and name of all results" << endl << endl;
//           exit(0);
//     }
     srand(seed);
//   srand((unsigned int)time(NULL));

//   alreadybest=1;
     flag_method=2;
//  deg_rot=left_r*0.01;
//   dalpha=ialpha*0.01;
     dalpha=1.0;

     dp=pro_nf1;
     dgama=taux_nf2;

     clockTicksPerSecond = (double)sysconf(_SC_CLK_TCK); // Get the clock ticks per second

//
//     cout << endl << "Name of fiche = " << benchmark << endl;
//     cout << "iteration = " << iteration << endl;
//     cout<<"---------------------------"<<endl;
     strcpy(name_final_result,benchmark);
     strcat(name_final_result,type_method);
     strcat(name_final_result,final);
     strcat(name_final_result,rep);

//     strcpy(name_all_results,name_fiche);
//     strcat(name_all_results,type_method);
//     strcat(name_all_results,all);
//     strcat(name_all_results,rep);

     read_fiche();

     n_02=0.2*max_vertex1;
     n_09=0.9*max_vertex1;
//     n_05=0.5*max_vertex1;
     n_04=0.4*max_vertex1;

     gama=n_09+rand()%n_02;
     L_0=0.05*max_vertex1;
     L_max=n_04+rand()%n_02;
     P_0=0.75;

     setdatastructure();
     Hillclimbing();
     freedatastructure();
     return 0;
}













