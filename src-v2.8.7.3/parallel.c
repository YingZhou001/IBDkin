#include "head.h"
#include "tools.h"
#include "print.h"
#include "parallel.h"

pthread_mutex_t mutex1;
pthread_mutex_t mutex2;



void read_buff(void)
{
    double a = my_wallclock();
    long int b = clock();
//fprintf(stderr, "\n\nthread-id%d, buff\n", (unsigned int)pthread_self());
    int i = 0;
    while(i < BUFF2 && (gzgets(ifp, strbuff2[i], BUFF-1) != Z_NULL)){
        i++;
    }
    buffi2 = i;
    Tgzread += clock() - b;
    t_gzread += my_wallclock() - a;
    return ;
}



void copy_buff(void)
{
    double a = my_wallclock();
    long int b = clock();

    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tAssign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, copy_buff_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n ERROR code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    for ( i = 0; i < Nthreads; i++ )pthread_join(thread_id[i], NULL);
    //fprintf(stderr, "done!\n");
    Tcpbuff += clock() - b;
    t_cpbuff += my_wallclock() - a;

    return ;
}

void *copy_buff_by_thread(void *args)
{
    int th_index = *((int *) args);
    int i;

    for(i = 0; i < buffi; i++){
        if(i % Nthreads != th_index)continue;
	strcpy(strbuff[i], strbuff2[i]);
    }
    return NULL;
}


void write_buff(void)
{
    double a = my_wallclock();
    long int b = clock();

    gzFile fp = gzopen(outfile, "a");
    gzbuffer(fp, BUFF3);
    int i;
    for(i = 0; i < buffi; i++){
	if(strbuff[i][0]!='\0')gzputs(fp, strbuff[i]);
    }
    gzclose(fp);
    Tgzputs += clock() - b;
    t_gzputs += my_wallclock() - a;
    return ;
}


void cal_pair_num(void)
{
    double a;
    long int b;
    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }

    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    rc = pthread_mutex_init(&mutex1, NULL);
    init_lock();
    if(rc){
	printf("\n ERROR code for pthread_mutex_init() is %d \n", rc);
	exit(1);
    }

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tAssign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, cal_pair_num_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n ERROR code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }

    pairNum = 0; segNum = 0;
    for ( i = 0; i < Nthreads; i++ ){
	pthread_join(thread_id[i], NULL);
	pairNum += pairNumByThread[i]; pairNumByThread[i] = 0;
	segNum += segNumByThread[i]; segNumByThread[i] = 0;
    }

    pthread_mutex_destroy(&mutex1);

    fprintf(stderr, "\t %ld IBD segments + %ld ID pairs are kept.\n", segNum, pairNum);
    //fprintf(stderr, "done!\n");
    if(checkTime == 1){
	t_cal_pair_num += my_wallclock() - a;
	Tcal_pair_num += clock() - b;
    }
    return ;
}

int cal_seg_num(hapIBD_t *head)
{
int segnum = 0;
while(head!=NULL){
segnum += 1;
head = head -> next;
}

//if(segnum>1)fprintf(stderr, "%d\n", segnum);
return segnum;
}

void *cal_pair_num_by_thread(void *args)
{
    int th_index = *((int *) args);
    Pair_t *cur;

    int i, j;

    for(i = 0; i < idNum; i++){

	pthread_mutex_lock(&mutex1);
	if(lock[i] == '1'){pthread_mutex_unlock(&mutex1);continue;}
	lock[i] = '1';
	pthread_mutex_unlock(&mutex1);

	pairNumByThread[th_index] += idhead[i]->num;
	if(idhead[i]->head != NULL){
	    for(j = 0; j < idhead[i]->num; j++){	
	    cur = (idhead[i]->head)[j];
	    segNumByThread[th_index] += cal_seg_num(cur->ibdhead);
	    }
	}
    }
    return NULL;
}





void shrink_id_pair(void)
{
    double a;
    long int b;
    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }

    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    rc = pthread_mutex_init(&mutex2, NULL);
    init_lock();
    if(rc){
	printf("\n ERROR code for pthread_mutex_init() is %d \n", rc);
	exit(1);
    }

//fprintf(stderr, "starting shrink_id_pair\n");

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tAssign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, shrink_id_pair_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n ERROR code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    for ( i = 0; i < Nthreads; i++ )pthread_join(thread_id[i], NULL);
//    fprintf(stderr, "done!\n");
    if(checkTime == 1){
	Tshrinke_id_pair += clock() - b;
	t_shrinke_id_pair += my_wallclock() - a;
	p_std_time(t_shrinke_id_pair);
	fprintf(stderr,"shrink_id_pair() Time = %s\n", timestr);

    }
    pthread_mutex_destroy(&mutex2);

    return ;
}

void *shrink_id_pair_by_thread(void *args)
{
    int th_index = *((int *) args);
    int *stack;
    int *value, *count;
    Pair_t *tmp;
    Pair_t **tmpp;
    int i, j, jj, k,num;

    for(i = 0; i < idNum; i++){

	pthread_mutex_lock(&mutex2);
	if(lock[i] == '1'){pthread_mutex_unlock(&mutex2);continue;}
	lock[i] = '1';
	pthread_mutex_unlock(&mutex2);


	stack = idhead[i]->stack;
	memSizeByThread[th_index] -= (stack[0]) * sizeof(int);
	num = idhead[i]->num;
	if(num == 0){
	    free(idhead[i]->stack);
	    idhead[i]->stack = NULL;
	    continue;
	}
	//stack = (int *)stack + 1;
	qsort(&(stack[1]), num, sizeof(int), int_cmp);


	value = (int *)calloc(num, sizeof(int));
	if(value == NULL){
	    my_error("fail to alloc mem for 'value'");
        }

	count = (int *)calloc(num, sizeof(int));
        if(count == NULL){
            my_error("fail to alloc mem for 'count'");
        }

	value[0] = stack[1];
	count[0] = 1;
	k = 0;
	for(j = 2; j < num+1; j++){
	    if(value[k] == stack[j])count[k] ++;
	    else {
		k++;
		value[k] = stack[j];
		count[k] = 1;
	    }
	}


	free(idhead[i]->stack);
	idhead[i]->stack = NULL;

	num = 0;
	for(j = 0; j < k+1; j++){if(count[j]>1)num ++;}
	idhead[i]->num = num;

	tmpp = (Pair_t **)calloc(num, sizeof(Pair_t *));
	if(tmpp == NULL){
	    my_error("fail to alloc mem for idhead[%d]->head, size = %u\n",i, num*sizeof(Pair_t *));
	}
	idhead[i]->head = tmpp;


	for(j = 0; j < num; j++){
	    tmp = calloc(1, sizeof(Pair_t));
	    if(tmp == NULL){
		my_error("fail to alloc mem for (idhead[%d]->head)[%d], size = %u\n",i,j, sizeof(Pair_t));
	    }
	    memSizeByThread[th_index] += sizeof(Pair_t) + sizeof(Pair_t *);
	    tmp -> index = -1;
	    tmp -> segnum = 0;
	    tmp -> s1 = 0;
	    tmp -> s2 = 0;
	    tmp -> ibdhead = NULL;
	    (idhead[i]->head)[j] = tmp;
	}
	jj = 0;
	for(j = 0; j < k+1; j++){
	    if(count[j]>1){
		tmp = (idhead[i]->head)[jj];
		tmp->index = value[j];
		jj ++;
	    }
	}
	free(value);
	free(count);
    }
    return NULL;
}




void cal_kinship(void)
{
    double a;
    long int b;

    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }
    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tassign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, cal_kinship_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n error code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    for ( i = 0; i < Nthreads; i++ )pthread_join(thread_id[i], NULL);
    if(checkTime == 1){ 
	Tcal_kinship += clock() - b;
	t_cal_kinship += my_wallclock() - a;
    }
    return ;
}

void *cal_kinship_by_thread(void *args)
{
    long int b = clock();
    int th_index = *((int *) args);
    IBD_t *ibdcur;
    int i, ii;
    int id1, id2;
    char h1, h2; 
    float g1, g2;
    for(i = 0; i < buffi; i++){

	ibdcur = IBDdat[i];
	if(ibdcur->pass != '1'){continue;}

	id1 = ibdcur->id1;
	id2 = ibdcur->id2;
	if(id1 % Nthreads != th_index)continue;

	g1 = ibdcur->g1;
	g2 = ibdcur->g2;

	h1 = ibdcur->h1;
	h2 = ibdcur->h2;

	//threadCount[th_index] ++;

	ii = binary_search_pair(idhead[id1], id2);
//fprintf(stderr, "%d %d\n", id2, ii);
	if(ii == -1){
	    if(g2 - g1 >= IBDcM){
		hapIBD_t *tmp = calloc(1, sizeof(hapIBD_t));
		tmp->h1 = h1;tmp->h2 = h2;
		tmp->g1 = g1;tmp->g2 = g2;
		tmp->next = NULL;
		output(tmp, id1, id2, i);
		free(tmp);}
	    else strbuff[i][0]='\0';
	}
	else {
	    memSizeByThread[th_index] += push_ibd(&((idhead[id1] -> head)[ii]->ibdhead), h1, h2, g1, g2);
	    strbuff[i][0]='\0';
	}
    }
    threadCount[th_index] += clock()/10000 - b;
    return NULL;
}



void output(hapIBD_t *head, int idi, int ids, int buf)
//kinship coefficient based on one segment
{	

    int segnum;
    float s0, s1, s2, kinship;
    segnum = 0; s0 = 0; s1 = 0; s2 = 0;
    segnum = calIBD12_pair(head, &s0, &s1, &s2);
    s1 = s1/totg;
    s0 = 1 - s1;
    kinship = s1/4;
    if(kinship >= kincut){
	sprintf(strbuff[buf], "%s %s %d %f %f 0 %f %d\n", idhead[idi]->id, idhead[ids]->id, segnum, s0, s1, kinship, cal_degree(kinship));
    }
    else strbuff[buf][0]='\0';
    return ;
}


void cal_kinship2(void)//for close relatives
{
    double a;
    long int b;

    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }
    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	rc = pthread_create(&thread_id[i], NULL, cal_kinship2_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n error code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    for ( i = 0; i < Nthreads; i++ ){
	pthread_join(thread_id[i], NULL);
    }
    if(checkTime == 1){
	Tcal_kinship += clock() - b;
	t_cal_kinship += my_wallclock() - a;
    }
    return ;
}

void *cal_kinship2_by_thread(void *args)
{
    int th_index = *((int *) args);
    Pair_t *prcur;
    int i, j;
    //int id1, id2;
    for(i = 0; i < idNum; i++){
	if(i % Nthreads != th_index)continue;
	for(j = 0; j < idhead[i]->num; j++){
	    prcur = (idhead[i]->head)[j];
	    output2(prcur);
	    //if(!(prcur->head))fprintf(stderr, "%d %d %d *%s* *%s*\n",i, id2, prcur->segnum, prcur->head, prcur->next);
	    memSizeByThread[th_index] += freehapIBDList(prcur->ibdhead);
	    prcur->ibdhead = NULL;
	}
    }
    return NULL;
}

void output2(Pair_t *head)
{

    int segnum;
    float s0, s1, s2;
    s0=0; s1 = 0; s2 = 0;
    segnum = calIBD12_pair(head->ibdhead, &s0, &s1, &s2);
    s1 = s1;
    s2 = s2;

    head-> segnum += segnum;
    head-> s1 += s1;
    head-> s2 += s2;

    return ;
}



void store_buff_ibd(void)
{
    double a;
    long int b;

    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }
    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;

    for ( i = 0; i < Nthreads; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tassign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, store_buff_ibd_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n error code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    for ( i = 0; i < Nthreads; i++ )pthread_join(thread_id[i], NULL);
    //fprintf(stderr, "done!\n");
    if(checkTime == 1){
	Tstore_buff_ibd += clock() - b;
	t_store_buff_ibd += my_wallclock() - a;
    }
    return ;
}

void *store_buff_ibd_by_thread(void *args)
{
    int b = clock()/10000;

    int th_index = *((int *) args);
    IBD_t *ibdcur;
    int i, chr;
    int id1, id2;
    float g1, g2;
    int num, L;
    for(i = 0; i < buffi; i++){

	ibdcur = IBDdat[i];
	if(ibdcur->pass != '1')continue;
	id1 = ibdcur->id1;
	id2 = ibdcur->id2;
	if(id1 % Nthreads != th_index)continue;

	//add length filters
	if(ibdcur->l < IBDcM)continue;//filters

	//update genetic ranges
	g1 = ibdcur->g1;
	g2 = ibdcur->g2;
	chr = ibdcur->chr;
	if(g1 < minpos[chr][th_index]){minpos[chr][th_index] = g1;}
	if(g2 > maxpos[chr][th_index]){maxpos[chr][th_index] = g2;}
	//fprintf(stderr, "%d %d %d %f %f %f %c\n", i, id1, id2, g1, g2, ibdcur->l, ibdcur->pass);


	//push pairs
	idhead[id1]->num += 1;
	num = idhead[id1]->num;
	L = (idhead[id1]->stack)[0];
	if(num > L - 1){
	    memSizeByThread[th_index] -= L * sizeof(int); 
	    while (num > L - 1) L = (int)(L * 1.25) + 5;
	    idhead[id1]->stack = realloc(idhead[id1]->stack, L * sizeof(int));
	    memSizeByThread[th_index] += L * sizeof(int); 
	    (idhead[id1]->stack)[0] = L;
	}
	(idhead[id1]->stack)[num] = id2;
    }
    threadCount[th_index] += clock()/10000 - b;
    return NULL;
}




void fill_buff_ibd_1st_pass(void)
{
    double a;
    long int b;

    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }


    pthread_t thread_id[Nthreads];
    int thread_args[Nthreads];
    int rc, i;
    int N;
    if (Nthreads > 1) N = Nthreads - 1;
    else N = Nthreads;


    for ( i = 0; i < N; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tAssign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, fill_buff_ibd_1st_pass_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n ERROR code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    read_buff();
    for ( i = 0; i < N; i++ ){
	rc = pthread_join(thread_id[i], NULL);
	if(rc > 0) printf("\n ERROR code from thread %d is %d \n", i, rc);

    }
    //fprintf(stderr, "done!\n");
    if(checkTime == 1){
	Tfill_buff_ibd1 += clock() - b;
	t_fill_buff_ibd1 += my_wallclock() - a;
    }
    return ;
}

void *fill_buff_ibd_1st_pass_by_thread(void *args)
{
    int th_index = *((int *) args);
    IBD_t *cur;
    int i;

    char str1[BUFF1], str2[BUFF1];
    int tmpid;
    int chr, id1, id2;
    float p1, p2, l;
    int N;


//fprintf(stderr, "thread-id%d:%d, fill-1, round %d\n",th_index, (unsigned int)pthread_self(), Round);
    if (Nthreads > 1) N = Nthreads-1;
    else N = Nthreads;

    for(i = 0; i < buffi; i++){
	if(i % N != th_index)continue;
	sscanf(strbuff[i], "%s %*s %s %*s %d %f %f %f", str1, str2, &chr, &p1, &p2, &l);
	if(l > IBDcM2) cal_coverage(th_index, chr, p1, p2);
	cur = IBDdat[i];

	if(Parts > 1){
	    if((hash_str(str1) + hash_str(str2)) % Parts == part - 1) cur->pass = '1';
	    else cur->pass = '0';
	}
	else cur->pass = '1'; //"1" means writing the pair to the memory, '0' means skipping this pair

	threadCount[th_index] ++;
	id1 = binary_search_ID_index(str1);
	id2 = binary_search_ID_index(str2);
	//fprintf(stderr, "%s: %d; %s: %d\n", str1, id1, str2, id2);

	if(cur->pass == '0')strbuff[i][0] = '\0';

	if(id1 == -1 || id2 == -1){
	    cur->pass = '0';
	    strbuff[i][0] = '\0';
	    cur->id1 = id1;
	    cur->id2 = id2;
	    continue;
	}

	if(id1 > id2){ //keep the order id1 < id2, for further ID pair assignment
	    tmpid = id1;
	    id1 = id2;
	    id2 = tmpid;
	}


	if(((id1 + id2) & 1) == 1) { 
	    /*permutate the searching order*/
	    /*after this step, id1 is the head index and id2 is the search index*/
	    tmpid = id1;
	    id1 = id2;
	    id2 = tmpid;
	}

	cur->id1 = id1;
	cur->id2 = id2;
	cur->g1 = p1;
	cur->g2 = p2;

	cur->chr = chr;
	cur->l = l;

	//fprintf(stderr, "%d %s %s %d %d\n",i, cur->str1, cur->str2, cur->id1, cur->id2);
    }
    return NULL;
}



void fill_buff_ibd_2nd_pass(void)
{
    double a;
    long int b;

    if(checkTime == 1){
	a = my_wallclock();
	b = clock();
    }
    pthread_t thread_id[Nthreads+1];
    int thread_args[Nthreads];
    int rc, i;
    int N;
if (Nthreads > 1) N = Nthreads - 1;
else N = Nthreads;

    for ( i = 0; i < N; i++ ){
	thread_args[i] = i;
	//fprintf(stderr, "\tAssign mission to thread-%d\n", i);
	rc = pthread_create(&thread_id[i], NULL, fill_buff_ibd_2nd_pass_by_thread, &thread_args[i]);
	if(rc){
	    printf("\n ERROR code from thread %d is %d \n", i, rc);
	    exit(1);
	}

    }
    read_buff();
    for ( i = 0; i < N; i++ ){
rc = pthread_join(thread_id[i], NULL);
if(rc > 0) printf("\n ERROR code from thread %d is %d \n", i, rc);
}

    //fprintf(stderr, "done!\n");
    if(checkTime == 1){
	Tfill_buff_ibd2 += clock() - b;
	t_fill_buff_ibd2 += my_wallclock() - a;
    }
    return ;
}

void *fill_buff_ibd_2nd_pass_by_thread(void *args)
{
    int th_index = *((int *) args);
    IBD_t *cur;
    int i;
    int N;

    int tmpid;
    int chr, id1, id2;
    float p1, p2, l;
    char h1, h2, tmph;
    char str1[BUFF1], str2[BUFF1];

//fprintf(stderr, "thread-id%d:%d, fill-2, round %d\n",th_index, (unsigned int)pthread_self(), Round);

    if (Nthreads > 1) N = Nthreads - 1;
    else N = Nthreads;

    for(i = 0; i < buffi; i++){
	if(i % N != th_index)continue;
	sscanf(strbuff[i], "%s %c %s %c %d %f %f", str1, &h1, str2, &h2, &chr, &p1, &p2);

	if(chr != CHR){
	    fprintf(stderr, "error! Detecting different chromosome numbers in one IBD file\n");
	    exit(-1);
	}


	cur = IBDdat[i];

	if(Parts > 1){
	    if((hash_str(str1) + hash_str(str2)) % Parts == part - 1) cur->pass = '1';
	    else cur->pass = '0';
	}
	else cur->pass = '1'; //"1" means pass to the memory, '0' means skip this pair

	threadCount[th_index] ++;
	id1 = binary_search_ID_index(str1);
	id2 = binary_search_ID_index(str2);

	if(cur->pass == '0')strbuff[i][0] = '\0';

	if(id1 == -1 || id2 == -1){
	    cur->pass = '0';
	    strbuff[i][0] = '\0';
	    cur->id1 = id1;
	    cur->id2 = id2;
	    continue;
	}

	if(id1 > id2){ //keep the order id1 < id2, for further ID pair assignment
	    tmpid = id1; tmph = h1;
	    id1 = id2; h1 = h2;
	    id2 = tmpid; h2 = tmph;
	}


	if(((id1 + id2) & 1) == 1) { 
	    //permutate the searching order
	    //after this step, id1 is the head index and id2 is the search index
	    tmpid = id1; tmph = h1;
	    id1 = id2; h1 = h2;
	    id2 = tmpid; h2 = tmph;
	}

	cur->id1 = id1;
	cur->id2 = id2;


	cur->h1 = h1;
	cur->h2 = h2;

	cur->chr = chr;
	cur->g1 = interpolate2cM(chr, p1);
	cur->g2 = interpolate2cM(chr, p2);


	l = cur->g2 - cur->g1 - apply_mask(chr, cur->g1, cur->g2);
	if(l < IBDcM2) {cur->pass = '0';strbuff[i][0] = '\0';}
	//fprintf(stderr, "%d %s %s %d %d\n",i, cur->str1, cur->str2, cur->id1, cur->id2);
    }
    return NULL;
}




int calIBD12_pair(hapIBD_t* head, float *S0, float *S1, float *S2)
{
    int L, segnum;
    int i, j, chr;
    int state;
    float from, to; 
    float a, b, mid, inc;
    int htag[2][2];
    char h1, h2;
    hapIBD_t *cur;
    float s0, s1,s2, chrL;
    s0 = 0; s1 = 0; s2 = 0; 


    cur = head;segnum = 0;
    chr = CHR;
    while(cur){
	segnum += 1;
	cur = cur-> next;
    }
    chrL = maxPos[chr] - minPos[chr] - apply_mask(chr, minPos[chr], maxPos[chr]);

    if(segnum == 0){
	s0 =  chrL;
	s1 = 0;
	s2 = 0;
    }
    else if(segnum == 1){
	s1 = head->g2 -head->g1 - apply_mask(chr, head->g1, head->g2);
	s0 = chrL - s1;
	s2 = 0;
    }
    else {
	cur = head;
	L = segnum * 2 + 2; //number of IBD segments
	float *segs; //store segments
	int *ibdstate; //store ibd0, 1, 2
	segs = (float *)calloc(L, sizeof(float));
	ibdstate = (int *)calloc(L, sizeof(int));
	segs[0] = minPos[chr];
	segs[L-1] = maxPos[chr];
	i = 1;
	cur = head;
	while(1){
	    segs[i] = cur->g1;
	    segs[i+1] = cur->g2;
	    i += 2;
	    if(cur->next == NULL)break;
	    cur = cur-> next;
	}

	//for(i = 0; i < L; i++)fprintf(stderr, "%f\n", segs[i]);
	//fprintf(stderr, " %d %d\n", L, segnum);
	qsort(segs, L, sizeof(float), float_cmp);

	//calculate ibd0, ibd1, ibd2
	for (i = 0; i < L-1; i++){
	    a = segs[i];b = segs[i+1];mid = (a + b)/2;
	    if(a > b){
		for(j = 0; j < L; j++)fprintf(stderr, " %f", segs[j]);
		fprintf(stderr, "%f %f\n", a, b);
	    }
	    htag[0][0] = 0;htag[0][1] = 0;
	    htag[1][0] = 0;htag[1][1] = 0;
	    cur = head;
	    while(cur){
		if(cur->g1 <=mid && cur->g2 >mid){
		    h1 = cur->h1;h2 = cur->h2;
		    if(h1 == haptag[0] && h2 == haptag[0])htag[0][0] = 1;
		    if(h1 == haptag[0] && h2 == haptag[1])htag[0][1] = 1;
		    if(h1 == haptag[1] && h2 == haptag[0])htag[1][0] = 1;
		    if(h1 == haptag[1] && h2 == haptag[1])htag[1][1] = 1;
		}
		cur = cur-> next;
	    }
	    if((htag[0][0] == 1 && htag[1][1] == 1) \
		    || (htag[0][1] == 1 && htag[1][0] == 1)) ibdstate[i] = 2;
	    else if(htag[0][0] == 0 && htag[1][1] == 0 \
		    && htag[0][1] == 0 && htag[1][0] == 0 ) ibdstate[i] = 0;
	    else ibdstate[i] = 1;

	}

//for(i = 0; i < L-1; i++){fprintf(stderr, "*%d %f %d\n",i, segs[i], ibdstate[i]);}

	//gapfilling
	for(i = 0; i < L-1; i++){
	    ////fill ibd0
	    if(ibdstate[i] == 0){
		j = i;while(j > 0 ){if(ibdstate[j] > 0)break;j--;}
		if( j == 0)from = segs[0] - gap1;
		else from = segs[j+1];

		j = i;while(j < L ){if(ibdstate[j] > 0)break;j++;}
		if(j == L)to=segs[j-1] + gap1;
		else to=segs[j];
		if((to - from) < gap1)ibdstate[i] = 1;
	    }

	    ////fill ibd2
	    
	    if(ibdstate[i] != 2){
		j = i;while(j > 0 ){if(ibdstate[j] == 2)break;j--;}
		if( j == 0)from = segs[0] - gap2;
		else from = segs[j+1];

		j = i;while(j < L ){if(ibdstate[j] == 2)break;j++;}
		if(j == L)to=segs[j-1] + gap2;
                else to=segs[j];
		if((to - from) < gap2)ibdstate[i] = 2;
//fprintf(stderr, "%d %f %f %f %d %d %d\n", i, from, to, segs[i], ibdstate[i], j, L);
	    }

	}
//for(i = 0; i < L-1; i++){fprintf(stderr, "**%d %f %d\n",i, segs[i], ibdstate[i]);}
//exit(-1);

	//calculte s0, s1, s2
	for (i = 0; i < L-1; i++){
	    a = segs[i];b = segs[i+1];state=ibdstate[i];
	    inc = b -a - apply_mask(chr, a, b);
	    if(state == 0) s0 +=inc;
	    if(state == 1) s1 +=inc;
	    if(state == 2) s2 +=inc;
	}
	free(segs);
	free(ibdstate);
    }
    *S0 += s0;
    *S1 += s1;
    *S2 += s2;
    return segnum;
}
