#include "head.h"
#include "tools.h"

int hash_str(char *str)
{
    int i;
    int l = strlen(str);
    int out;
    long int sum = 0;
    const int p = 7;
    for(i = 0; i < l; i++)sum += (int)str[i]*(int)pow(p,i);
    out = (int)abs(sum) % 991817;
    return out;
}

int ID_t_cmp(const void *a, const void *b)
{
    const ID_t *ia = *(const ID_t **)a;
    const ID_t *ib = *(const ID_t **)b;
    return(strcmp(ia->id, ib->id));
}

int int_cmp(const void *a, const void *b)
{
    const int *ia = (const int *)a; // casting pointer 
    const int *ib = (const int *)b;
    return (*ia  - *ib);
}


int float_cmp(const void *a, const void *b)
{
    const float *ia = (const float *)a; // casting pointer 
    const float *ib = (const float *)b; 
    if ((*ia  - *ib) >= 0) return 1;
    else return -1;
}



int double_cmp(const void *a, const void *b)
{
    const double *ia = (const double *)a; // casting pointer 
    const double *ib = (const double *)b; 
    if ((*ia  - *ib) >= 0) return 1;
    else return -1;
}


int binary_search_ID_index(char *id)
{
    if(strcmp(id, idhead[0]->id) < 0 || strcmp(id, idhead[idNum-1]->id) > 0)return -1;
    else {
	int tmp, i;
	int from, to, new;
	char *str;
	from=0; to= idNum-1;
        while(1){
            new = (int)(from+to)/2;
	    str = idhead[new]->id;
	    tmp = strcmp(id, str);
	    if(tmp == 0) return new;
	    else if(tmp > 0)from = new;
	    else to = new;
	    if(to - from < 3){
		for(i = from; i <= to; i++){
		    if(strcmp(id, idhead[i]->id)==0)return i;
		}
		break;
	    }
	}
	return -1;
    }
}


int binary_search_pair(ID_t *id, int index)
{
    int num = id->num;
    if(num == 0)return -1;
    Pair_t **cur =id->head;

    //fprintf(stderr, "target %d\n", index);
    //for(int j = 0; j < num; j++){fprintf(stderr, "search %d\n", cur[j]->index);}
    if(index < (cur[0]->index) || index > (cur[num-1]->index))return -1;
    else {
	int tmp, i;
	int from, to, new;
	from=0; to= num-1;
	while(1){
            new = (int)(from+to)/2;
            tmp = index - cur[new]->index;
            if(tmp == 0) return new;
            else if(tmp > 0)from = new;
            else to = new;
            if(to - from < 3){
                for(i = from; i <= to; i++){
                    if((index - cur[i]->index)==0)return i;
                }
                break;
            }
	}
	return -1;
    }

}


size_t push_ibd(hapIBD_t **head, char h1, char h2, float g1, float g2)
{
    hapIBD_t *cur = (hapIBD_t*)calloc(1, sizeof(hapIBD_t));
    cur->h1 = h1;
    cur->h2 = h2;
    cur->g1 = g1;
    cur->g2 = g2;
    cur->next = (*head);
    (*head) = cur;
    return sizeof(hapIBD_t);
}

int binary_search(float pos, float *map, int mapl)
{
    if(pos < map[0]) return -1;
    else if (pos >= map[mapl - 1]) return mapl-1;
    else {
	int tag = 0;
	int i_old1, i_old2, i_new;
	i_old1=0; i_old2= mapl-1;
	while(tag == 0){
	    i_new = (int)(i_old1+i_old2)/2;
	    if(pos >= map[i_new] && pos < map[i_new+1])return(i_new);
	    else if(pos < map[i_new])i_old2 = i_new;
	    else if(pos >= map[i_new+1])i_old1 = i_new;
	    if(i_old2 == i_old1)break;
	}
    }
    fprintf(stderr, "Error in genetic position interpolation, please check map file\n");
    exit(-1);
    return -1;
}

float interpolate2cM(int chr, float pos)
{
    /*interpolate genetic position*/
    int k, kk, mapl;
    float out, rate;
    float *mapp, *mapg;
    mapp = mapP[chr];
    mapg = mapG[chr];
    mapl = mapL[chr];
    k=binary_search(pos, mapp, mapl);
    if(k<0 || k>(mapl-2)){
	if(k < 0)kk = binary_search(mapp[0] + 5000000, mapp, mapl);
	if(k > (mapl-2))kk = binary_search(mapp[mapl-1] - 5000000, mapp, mapl);
	/*using chromosome average if the marker out of map coverage*/
	if(kk < 0 || kk > (mapl-1)){
	    rate = (mapg[mapl-1]-mapg[0])/(mapp[mapl-1]-mapp[0]);
	    out = mapg[0]+(pos-mapp[0])*rate;
	}
	/*using 5Mb at chromosome ends*/
	else if(k < 0){
	    rate = (mapg[kk] - mapg[0]) / (mapp[kk] - mapp[0]);
	    out = mapg[0]+(pos-mapp[0])*rate;
	}
	else if(k > (mapl-2)){
	    rate = (mapg[mapl-1] - mapg[kk]) / (mapp[mapl-1] - mapp[kk]);
	    out = mapg[mapl-1] + (pos - mapp[mapl-1]) * rate;
	}
    }
    else {
	rate = (mapg[k+1]-mapg[k])/(mapp[k+1]-mapp[k]);
	out=mapg[k]+(pos-mapp[k])*rate;
    }
    //fprintf(stderr, "chr= %d pos = %f out = %f\n", chr, pos, out);
    return out;
}

double my_wallclock(void)
{
    double T;
    struct timeval timecheck;

    gettimeofday(&timecheck, NULL);
    T = timecheck.tv_sec  + (double)timecheck.tv_usec / 1000000;
    return T;
}

void init_lock(void)
{
    int i;
    for(i = 0; i < idNum; i++)lock[i]='0';
    return ;
}



void cal_coverage(int th_index, int chr, float P1, float P2)
{
    double **tf;
    int i, j;
    int k;
    float x, y;
    float p1 = P1 / (BINSIZE * 1000);
    float p2 = P2 / (BINSIZE * 1000);

    x = y = 0.000001;

    i = (int)floor(p1) + 1;
    j = (int)floor(p2);

    tf = coverage_per_thread[th_index]->coverage;

    if(j == i -1)tf[chr][j] += p2 - p1;
    if(i <= j){
	x = (float)(i - p1);
	tf[chr][i - 1] += x;
	y = (float)(p2 - j);        
	tf[chr][j] += y;
	for(k = i; k < j; k ++)tf[chr][k] += 1;
    }


    if(x>1 || y >1){
	fprintf(stderr, "input %d %d %.0f %.0f\n",th_index, chr, P1, P2);
	fprintf(stderr, "out %f %f %f=%f %f=%f %d %d\n",p1, p2, x, i-p1, y, p2-j, i, j);
	exit(-1);
    }
    return ;
}



void cal_coverage_median(void)
{
    int i, j;
    int chr;
    unsigned int L;
    double **tf;
    L = 0;
    for(chr = 1; chr <= 22; chr++){
	for(i = 0; i < coverageL[chr]; i ++){
	    //if(chr == 21)fprintf(stderr, "%d", i);
	    for(j = 0; j < Nthreads; j++){
		tf = coverage_per_thread[j]->coverage;
		coverage[chr][i] += tf[chr][i];
		//if(chr==21)fprintf(stderr, " %lf", tf[chr][i]);
	    }
	    //if(chr==21)fprintf(stderr, "\n");
	    if(coverage[chr][i]>0)L++;
	}
    }
    double *segs; /*store segments*/
    segs = (double *)calloc(L, sizeof(double));
    j = 0;
    for(chr = 1; chr <= 22; chr++){
	for(i = 0; i < coverageL[chr]; i ++){
	    if(coverage[chr][i]>0){segs[j] = coverage[chr][i];j++;}
	}
    }

    qsort(segs, L, sizeof(double), double_cmp);
    if(L % 2 == 0)median = (segs[(int)(L/2)] + segs[(int)(L/2)-1])/2;
    else median = segs[(int)(L/2)];

    fprintf(stderr, "\nGenome-wide coverage median: %lf\n", median);

    free(segs);
    return ;
}

size_t push_mask(msk_t **ref, int p1, int p2, float g1, float g2)
{
    msk_t* cur = (msk_t *) calloc (1, sizeof(msk_t));
    cur->p1 = p1;
    cur->p2 = p2;
    cur->g1 = g1;
    cur->g2 = g2;

    cur->next = (*ref);
    (*ref) = cur;
    return sizeof(msk_t);

}

void delete_mask(msk_t *head)
{
    msk_t* cur;
    while(head){
	cur = head;
	head = head ->next;
	free(cur);
    }
    return ;
}

void create_mask(void)
{
    int i;
    int chr;
    int p1, p2;
    int newp1, newp2;
    float g1, g2;
    for(chr = 1; chr <= 22; chr++){
	p1 = -1 ; p2 = -1;
	for(i = coverageL[chr] - 1; i >= 0; i --){
	    if(coverage[chr][i] >0 && \
		    (coverage[chr][i]> 4*median || coverage[chr][i] < median/4)){
		newp1 = (int) (i*BINSIZE*1000);
		newp2 = (int) ((i+1)*BINSIZE*1000);
		if(p1 == -1) p1 = newp1;
		if(p2 == -1) p2 = newp2;

		if(newp2 == p1) p1 = newp1;
		else if(newp1 != p1){
		    g1 = interpolate2cM(chr, p1);
		    g2 = interpolate2cM(chr, p2);
		    memSize += push_mask(&mask[chr], p1, p2, g1, g2);
		    //printf("%d:%d %d\n", chr, p1, p2);
		    p1 = newp1;
		    p2 = newp2;
		}
	    }
	}
	if(p1 != -1){
	    g1 = interpolate2cM(chr, p1);
	    g2 = interpolate2cM(chr, p2);
	    memSize += push_mask(&mask[chr], p1, p2, g1, g2);
	    //printf("%d:%d %d\n", chr, p1, p2);
	}
    }
    return ;


}




float apply_mask(int chr, float g1, float g2)
{

    if(g1 > g2 ) {
	fprintf(stderr, "error in apply_mask(), g1 > g2\n");
	exit(-1);
    }
    msk_t *cur;
    cur = mask[chr];
    float out;
    float mfrom, mto;
    if(cur == NULL)return 0;
    out = 0;
    while(cur){
	mfrom = cur->g1;
	mto = cur->g2;

	if(g1 >= mfrom && g2 < mto)out += g2 - g1;
	else if(g1 >= mfrom && g1 < mto && g2 > mto)out += mto - g1;
	else if(g1 < mfrom && g2 > mto)out += mto - mfrom;
	else if(g1 < mfrom && g2 >= mfrom && g2 < mto)out += g2 - mfrom;
	else {}
	cur = cur->next;
    }
    return out;
}


int cal_degree(float kinship)
{
    int degree;
    float a;
    degree = 0;
    while(1){
	a = pow(2, -degree-1.5);
	if(a<=kinship) return degree;
	degree ++;
    }
}


int check_chr(void)
{
    IBD_t *ibdcur;
    register int i, chr_old, chr;
    chr_old = -1; chr = -1;
    for(i = 0; i < buffi; i++){
	ibdcur = IBDdat[i];
	if(ibdcur->pass != '1'){continue;}
	chr = ibdcur->chr;
	if(chr_old == -1) chr_old = chr;
	if(chr_old != chr){
	    fprintf(stderr, "error! Detecting different chromosome numbers in one IBD file\n");
	    exit(-1);
	}
    }
    return chr;
}

long int freehapIBDList(hapIBD_t *head)
{
    hapIBD_t* tmp;
    long int size = 0;

    while (head != NULL)
    {
	tmp = head;
	head = head->next;
	free(tmp);
	size -= sizeof(hapIBD_t);
    }
    return size;
}



void free_all(void)
{
    int i, j;
    free(headfile);
    for(i = 0; i < Nfile; i ++)free(ibdfile[i]);free(ibdfile);

    free(mapfile);
    free(maskfile);
    free(coveragefile);

    free(idfile);

    free(threadCount);
    free(memSizeByThread);
    free(pairNumByThread);
    free(segNumByThread);

    free(outfile);

    for(i = 0; i < 23; i ++){
	free(minpos[i]);
	free(maxpos[i]);
	free(coverage[i]);
	delete_mask(mask[i]);
    }
    free(minpos);
    free(maxpos);
    free(coverage);
    free(mask);

    for(i = 0; i < BUFF2 ;i ++){
	free(IBDdat[i]);
	free(strbuff[i]);
	free(strbuff2[i]);
    }
    free(IBDdat);
    free(strbuff);
    free(strbuff2);

    for (i = 0; i < Nthreads; i++){
	for(j = 0; j < 23; j ++)free((coverage_per_thread[i]->coverage)[j]);
	free(coverage_per_thread[i]->coverage);
	free(coverage_per_thread[i]);
	free(outbuff[i]);
    }
    free(coverage_per_thread);
    free(outbuff);
    free(outbuffi);

    free(lock);
    for(i = 0; i < idNum; i ++ ){
	if(tagKinship == 1){
	    for(j = 0; j < idhead[i]->num; j++){
		free((idhead[i]->head)[j]);
	    }
	    free(idhead[i]->head);
	}
	if(idhead[i]->stack != NULL)free(idhead[i]->stack);
	free(idhead[i]->id);
	free(idhead[i]);
    }
    free(idhead);

    for(i = 1; i<= 22; i++){
	free(mapP[i]);
	free(mapG[i]);
    }
    free(mapP);
    free(mapG);
}


void my_error(const char * message, ...)
{
  va_list args;
  va_start (args, message);
  vfprintf (stderr, message, args);
  va_end (args);
  fprintf (stderr, "\n");
  exit (EXIT_FAILURE);
}
