#include "head.h"
#include "tools.h"
#include "print.h"
#include "read.h"
#include "parallel.h"

void read_headfile(char *headfile)
{
    char str[BUFF_filename] = "0";
    int i;
    Nfile = 0;
    FILE *ifp = fopen(headfile, "r");
    while(fscanf(ifp, "%s",str)==1){
	fprintf(stderr, "  %s", str);
	if(access(str, R_OK) == -1) {fprintf(stderr, " does not exist!\n"); exit(-1);}
	fprintf(stderr, "\n");
	Nfile++;}
    fclose(ifp);
    ibdfile = (char **)calloc(Nfile, sizeof(char*));
    memSize += Nfile * sizeof(char*) + sizeof(char**);

    ifp = fopen(headfile, "r"); i = 0;
    while(fscanf(ifp, "%s",str)==1){
	ibdfile[i] = (char * )calloc(strlen(str)+10, sizeof(char));
	memSize += (strlen(str)+10) * sizeof(char);
	strcpy(ibdfile[i], str);

	i ++;
    }
    fclose(ifp);
}


int check_input(int arg, char**argv)
{
    fprintf(stderr, "############################\n");
    fprintf(stderr, "#Welcome to IBDkin V2.8.7.7#\n");
    fprintf(stderr, "############################\n");
    fprintf(stderr, "Copyright (C) 2020 Ying Zhou\n\n");
    //fprintf(stderr, "zlib verison: %s\n\n", ZLIB_VERSION);
    int i, fsize;
    IBDcM = 0;
    checkTime = 0; checkMem = 0;
    BINSIZE = 1000; FOLD = 4;
    IBDcM = 4; IBDcM2 = 2; gap1 = 5; gap2 = 20;
    kincut = 0.000001;
    Parts = 1; part = 1;
    tagKinship = 1;
    tagMask = 0;
    tagCoverage = 0;
    Nthreads = 2;
    outfile = (char *) calloc (16, sizeof(char));
    maskfile = (char *) calloc (16, sizeof(char));
    coveragefile = (char *) calloc (16, sizeof(char));
    memSize += 3* 16 * sizeof(char);
    strcpy(outfile, "kinship.gz");
    strcpy(maskfile, "mask.gz");
    strcpy(coveragefile, "coverage.gz");



    if(arg <=1)print_help();

    for(i = 1; i< arg; i++){
	if(strcmp(argv[i], "--ibdfile")==0){
	    fsize = strlen(argv[i+1])+10;
	    headfile = (char *) calloc (fsize, sizeof(char));
	    memSize += (fsize) * sizeof(char);
	    strcpy(headfile, argv[i+1]); i = i + 1;
	}
	else if(strcmp(argv[i], "--map")==0){
	    fsize = strlen(argv[i+1])+10;
	    mapfile = (char *) calloc (fsize, sizeof(char));
	    memSize += (fsize) * sizeof(char);
	    strcpy(mapfile, argv[i+1]); i = i + 1;
	}
        else if(strcmp(argv[i], "--range")==0){
            fsize = strlen(argv[i+1])+10;
            rangefile = (char *) calloc (fsize, sizeof(char));
            memSize += (fsize) * sizeof(char);
            strcpy(rangefile, argv[i+1]); i = i + 1;
        }
	else if(strcmp(argv[i], "--ind")==0){
	    fsize = strlen(argv[i+1])+10;
	    idfile = (char *) calloc (fsize, sizeof(char));
	    memSize += (fsize) * sizeof(char);
	    strcpy(idfile, argv[i+1]); 
	    i = i + 1;
	}
	else if(strcmp(argv[i], "--out")==0){
	    fsize = strlen(argv[i+1])+20;
            outfile = realloc (outfile, fsize*sizeof(char));
            maskfile = realloc (maskfile, fsize*sizeof(char));
            coveragefile = realloc (coveragefile, fsize*sizeof(char));
            memSize += 3* (fsize - 16) * sizeof(char);
            sprintf(outfile, "%s.kinship.gz", argv[i+1]); 
            sprintf(maskfile, "%s.mask.gz", argv[i+1]); 
            sprintf(coveragefile, "%s.coverage.gz", argv[i+1]); 
	    i = i + 1;
        }

	else if(strcmp(argv[i], "--cutcm")==0){
	    IBDcM = atof(argv[i+1]);
	    IBDcM2 = atof(argv[i+2]);
	    i = i + 2;
	}
	else if(strcmp(argv[i], "--binkb")==0){
            BINSIZE = atof(argv[i+1]);i = i + 1;
        }
	else if(strcmp(argv[i], "--fold")==0){
            FOLD = atof(argv[i+1]);i = i + 1;
        }
	else if(strcmp(argv[i], "--kinship")==0){
            kincut = atof(argv[i+1]);i = i + 1;
        }
	else if(strcmp(argv[i], "--nthreads")==0){
	    Nthreads = atoi(argv[i+1]);i = i + 1;
	}
	else if(strcmp(argv[i], "--checkTime")==0){
	    checkTime = 1;
	}
	else if(strcmp(argv[i], "--checkMem")==0){
            checkMem = 1;
        }
	else if(strcmp(argv[i], "--nokinship")==0){
	    tagKinship = 0;
	}
	else if(strcmp(argv[i], "--outmask")==0){
	    tagMask = 1;
	}
	else if(strcmp(argv[i], "--outcoverage")==0){
	    tagCoverage = 1;
	}
	else if(strcmp(argv[i], "--part")==0){
	    Parts = atoi(argv[i+1]);
	    part = atoi(argv[i+2]);
	    i = i + 2;
	}
	else if(strcmp(argv[i], "--merge")==0){
	    gap1 = atof(argv[i+1]);
	    gap2 = atof(argv[i+2]);
	    i = i + 2;
	}

	else {
	    print_help();
	}
    }
    if( tagKinship + tagMask + tagCoverage == 0 ){
	fprintf(stderr, "No output is specified!\n");
	print_help();
    }

    return 0;
}


void write_relate(void)
{
    double a = my_wallclock();
    long int b = clock();
    gzFile fp = gzopen(outfile, "a");
    gzbuffer(fp, BUFF_gz);
    Pair_t *prcur;
    int i, j, idi, ids, segnum;
    float kinship, s0, s1, s2;
    char str[BUFF_col]; //for storing each output row
    int tmp;
    char degree[8];
    for(i = 0; i < idNum; i++){
	for(j = 0; j < idhead[i]->num; j++){
	    prcur = (idhead[i]->head)[j];
	    idi = i;
	    ids = prcur->index;
	    segnum = prcur->segnum;;
	    s1 = prcur->s1 / totg;
	    s2 = prcur->s2 / totg;
	    s0 = 1 - s1 - s2;
	    kinship = s1/4+s2/2;
	    if(kinship >= kincut){
		tmp = cal_degree(kinship);
		if(tmp <= 3)sprintf(degree, "%d", tmp);
		else sprintf(degree, ">3");
		if(s2 <= 0.00001)sprintf(str, "%s\t%s\t%d\t%f\t%f\t0\t%f\t%s\n", idhead[idi]->id, idhead[ids]->id, segnum, s0, s1, kinship, degree);
		else sprintf(str, "%s\t%s\t%d\t%f\t%f\t%f\t%f\t%s\n", idhead[idi]->id, idhead[ids]->id, segnum, s0, s1, s2, kinship, degree);
		gzputs(fp, str);
	    }
	}
    }
    gzclose(fp);
    Tgzputs += clock() - b;
    t_gzputs += my_wallclock() - a;
    return ;
}

void read_ibd_1st_pass(void)
{
    double a;
    float maskL;
    Round = 0;
    if(checkTime == 1 )a = my_wallclock();


    int i, j;
    char str[BUFF_id];
    for ( i = 1; i < chrN; i++ ){
	minPos[i] = INT_MAX; //bp
	maxPos[i] = 0; //bp
    }
    for (i = 0; i < Nfile; i++){
	ifp = gzopen(ibdfile[i], "r");
	gzbuffer(ifp, BUFF_gz);
	fprintf(stderr, "\treading %s\n", ibdfile[i]);
	read_buff();

	sscanf(strbuff2[0], "%*s %*s %*s %*s %s", str);
	CHR = binary_search_string(chrs, chrN, str);

	if(CHR < 1 || CHR >= chrN){
	    fprintf(stderr, "unrecognized chromsome identifier: %s, skipping!\n", str);
	    continue;
	}

	while(1){
	    Round += 1;
//fprintf(stderr, "round = %d\n", Round);
	    if(buffi2 == 0)break;
	    buffi = buffi2;
//fprintf(stderr, "copy_buff()\n");
	    copy_buff();
//fprintf(stderr, "fill_buff_ibd_1st_pass()\n");
	    fill_buff_ibd_1st_pass();
	    //p_buff();
//fprintf(stderr, "store_buff_ibd()\n");
	    store_buff_ibd();
	}
	gzclose(ifp);
    }
    //p_datamap();
    //p_pair_depth();
    //cal_pair_num();

    cal_coverage_median();
    create_mask();

    fprintf(stderr, "\n  IBD coverages:\n");
    totg = 0; maskL = 0;
    for ( i = 1; i < chrN; i++ ){
	for (j = 0; j < Nthreads; j++){
	    if(minPos[i] > minpos[i][j]) minPos[i] = minpos[i][j];
	    if(maxPos[i] < maxpos[i][j]) maxPos[i] = maxpos[i][j];
	}
	if( minPos[i] < maxPos[i]){
	    minPos[i] = interpolate2cM(i, minPos[i]);
	    maxPos[i] = interpolate2cM(i, maxPos[i]);
	    //totg += maxPos[i] - minPos[i];
	    maskL = apply_mask(i, minPos[i], maxPos[i]);
	    totg += maxPos[i] - minPos[i] - maskL;
	    fprintf(stderr, "\t%s: %.2f -> %.2f cM, masked length = %.2fcM\n",chrs[i], minPos[i], maxPos[i], maskL);
	}

    }
    fprintf(stderr, "  Total chromosome length: %f cM\n", totg);
    if(checkTime ==1 ){
	p_std_time(my_wallclock() - a);
	fprintf(stderr,"\nread_ibd_1st_pass() Time = %s\n", timestr);}

    return ;
}


void read_ibd_2nd_pass(void)
{
    double a;
    if(checkTime)a = my_wallclock();
    int i;
    char str[BUFF_id];

    Round = 0;

    //ofp = fopen(outfile, "a");
    //gzFile ifp;
    for (i = 0; i < Nfile; i++){
	ifp = gzopen(ibdfile[i], "r");
	fprintf(stderr, "\treading %s\n", ibdfile[i]);
	read_buff();

        sscanf(strbuff2[0], "%*s %*s %*s %*s %s", str);
        CHR = binary_search_string(chrs, chrN, str);

        if(CHR < 1 || CHR >= chrN){
            fprintf(stderr, "unrecognized chromsome identifier: %s, skipping!\n", str);
            continue;
        }

	while(1){
	    Round += 1;
	    if(buffi2 == 0)break;
	    buffi = buffi2;
	    copy_buff();
	    fill_buff_ibd_2nd_pass();
	    //p_buff();
	    cal_kinship();//for distant relatives
	    write_buff();
	}
	cal_pair_num();
	p_mem();
	cal_kinship2();//for close relatives
	gzclose(ifp);
    }
    //fclose(ofp);
    write_relate();

    //p_datamap();
    if(checkTime){
	p_std_time(my_wallclock() - a);
	fprintf(stderr,"read_ibd_2nd_pass() Time = %s\n", timestr);
    }
}


void read_range(void)
{
    int i,j;
    chrs = (char **)calloc(100, sizeof(char *));
    for(i = 0; i < 100; i++){
        chrs[i] = (char *)calloc(32, sizeof(char));
        chrs[i][0] = '\0';
    }
    chrN = 1;
    gzFile fp;
    int chr, to;
    char str[BUFF_col];
    fp = gzopen(rangefile, "r");
    while(gzgets(fp, str, BUFF_col-1)!=NULL){
	sscanf(str, "%s %*s %d", chrs[99], &to);
	if(strcmp(chrs[99], chrs[chrN-1]) != 0){
	    strcpy(chrs[chrN], chrs[99]);
	    chrN += 1;
	}
    }
    gzclose(fp);
    if(chrN >= 100){fprintf(stderr, "Error, please contact the author!\n");exit(-1);}

    for(i = chrN; i < 100; i++)free(chrs[i]);
    chrs = realloc(chrs, chrN*sizeof(char *));

    qsort(chrs, chrN, sizeof(char *), string_cmp);

    //for(i=0;i<chrN;i++){fprintf(stderr, "%d %s %d\n", i, chrs[i], binary_search_string(chrs, chrN, chrs[i]));}
    //exit(-1);

    minpos = (float **)calloc(chrN, sizeof(float *));
    maxpos = (float **)calloc(chrN, sizeof(float *));
    for(i = 0; i < chrN; i++){
        minpos[i] = (float *)calloc(Nthreads, sizeof(float));
        maxpos[i] = (float *)calloc(Nthreads, sizeof(float));
        for(j = 0; j < Nthreads; j++){
            minpos[i][j] = INT_MAX;
            maxpos[i][j] = -1;
        }
    }
    memSize += 2 * chrN * (sizeof(float *) + Nthreads * sizeof(float));

    float bin = BINSIZE * 1000;
    char sstr[BUFF_id];
    fp = gzopen(rangefile, "r");
    coverageL[0] = 1;
    while(gzgets(fp, str, BUFF_col-1)!=NULL){
        sscanf(str, "%s %*s %d", sstr, &to);
	chr = binary_search_string(chrs, chrN, sstr);
	coverageL[chr] = (int) (to / bin) + 5;
    }
    gzclose(fp);

    /*the final coverage file*/
    coverage = (double **)calloc(chrN, sizeof(double*));
    memSize += sizeof(double **);
    for(chr = 0; chr < chrN; chr++){
	coverage[chr] = (double *)calloc(coverageL[chr], sizeof(double));
	memSize += sizeof(double *) + coverageL[chr] * sizeof(double);
	for(i = 0; i < coverageL[chr]; i ++){
	    coverage[chr][i] = -0.00000001;
	}
    }

    /*store coverage for each thread*/
    double **tf;
    coverage_per_thread = (cvrg_t **)calloc(Nthreads, sizeof(cvrg_t *));
    memSize += sizeof(cvrg_t **);
    for (i = 0; i < Nthreads; i++){
	coverage_per_thread[i] = (cvrg_t *)calloc(1, sizeof(cvrg_t));
	tf = (double **)calloc(chrN, sizeof(double*));
	memSize += sizeof(cvrg_t) + chrN * sizeof(double*);
	for(chr = 0; chr < chrN; chr++){
	    tf[chr] = (double *)calloc(coverageL[chr], sizeof(double));
	    memSize += coverageL[chr] * sizeof(double);
	    for(j = 0; j < coverageL[chr]; j ++){
		tf[chr][j] = -0.00000001;
	    }
	}
	coverage_per_thread[i]->coverage = tf;
    }

    mask = (msk_t **)calloc(chrN, sizeof(msk_t *));
    memSize += sizeof(msk_t **);
    for(chr = 0; chr < chrN; chr++){
	mask[chr] = (msk_t *)calloc(1, sizeof(msk_t));
	memSize += sizeof(msk_t);
	mask[chr]->p1 = -1;
	mask[chr]->p2 = -1;
	mask[chr]->g1 = -1;
	mask[chr]->g2 = -1;
	mask[chr]->next = NULL;
    }



    return ;
}

void read_map(char *plinkmap)
{
    double a = my_wallclock();
    gzFile fp;
    char str[BUFF_col], sstr[BUFF_id];
    int i, chr, j;//chr is the chromsome index

    float p, g;
    int mapi[chrN];
    for(i = 0; i < chrN; i++)mapL[i]=0;
    fp=gzopen(plinkmap, "r");
    while(gzgets(fp, str, BUFF_col-1)!=NULL){
	chr = -1;
	if(sscanf(str, "%s", sstr) == 1){
	    chr = binary_search_string(chrs, chrN, sstr);
	    if(chr > 0)mapL[chr] ++;
	}
    }
    gzclose(fp);

    mapP = (float **)calloc(chrN, sizeof(float*));
    mapG = (float **)calloc(chrN, sizeof(float*));
    memSize += 2*chrN*sizeof(float*);
    for(i = 0; i< chrN; i++){
	mapP[i] = (float *)calloc(mapL[i], sizeof(float));
	mapG[i] = (float *)calloc(mapL[i], sizeof(float));
	memSize += 2*mapL[i]*sizeof(float);
    }
    fp=gzopen(plinkmap, "r");
    i=0;
    for(chr = 1; chr < chrN; chr++)mapi[chr]=0;
    while(gzgets(fp, str, BUFF_col-1)!=NULL){
	if(sscanf(str, "%s %*s %f %f", sstr, &g, &p)==3){
	    chr = binary_search_string(chrs, chrN, sstr);
	    if(chr > 0 && chr < chrN){
		j = mapi[chr];
		mapG[chr][j] = g;
		mapP[chr][j] = p;
		mapi[chr] ++;
	    }
	}
    }
    gzclose(fp);
    /*check input map*/
    for (chr = 1; chr < chrN; chr ++){
	for(i=0;i<(mapL[chr]-1);i++){
//fprintf(stderr, "%d %s %.0f %f\n", chr, chrs[chr], mapP[chr][i], mapG[chr][i]);
	    if(mapP[chr][i+1]<mapP[chr][i]||mapG[chr][i+1]<mapG[chr][i]){
		fprintf(stderr, "Error: map position should be in increasing order!\n error in line %d %.0f %.0f\n",chr, mapP[chr][i], mapP[chr][i+1]);
		exit(-1);
	    }
	}
    }
    p_std_time(my_wallclock() - a);
    fprintf(stderr,"read_genetic_map() Time = %s\n", timestr);
    return ;
}

void read_ind(char *idfile)
{
    double a = my_wallclock();
    gzFile fp;
    int i;
    char str[BUFF_id];
    fp = gzopen(idfile, "r"); idNum = 0;
    while(gzgets(fp, str, BUFF_id-1)!=NULL)idNum ++;
    gzclose(fp);
    idhead = (ID_t **)calloc(idNum, sizeof(ID_t*));
    lock = (char *)calloc(idNum, sizeof(char));
    memSize += idNum * (sizeof(ID_t*) + sizeof(char)) + sizeof(ID_t**);
    fp = gzopen(idfile, "r"); i = 0;
    while(gzgets(fp, str, BUFF_id-1)!=NULL){
	idhead[i] = (ID_t *)calloc(1, sizeof(ID_t));
	idhead[i]->id = (char *)calloc(strlen(str)+1,sizeof(char));
	idhead[i]->head = NULL;
	idhead[i]->num = 0; 
	idhead[i]->stack = (int *)calloc(4, sizeof(int));
	(idhead[i]->stack)[0] = 4;
	memSize += sizeof(ID_t) + strlen(str)*sizeof(char) + 4 * sizeof(int);
	sscanf(str, "%s", idhead[i]->id);
	i++;
    }
    gzclose(fp);
    qsort(idhead, idNum, sizeof(ID_t *), ID_t_cmp);

    for(i = 0; i < (idNum-1); i++){
	if(strcmp(idhead[i]->id, idhead[i+1]->id)==0){
	    fprintf(stderr, "duplicated individuals: %s\n", idhead[i]->id);
	    exit(-1);
	}
    }

    fprintf(stderr, "%d individuals\n", idNum);
    p_std_time(my_wallclock() - a);
    fprintf(stderr,"read_ind() Time = %s\n", timestr);
    return ;
}


void init(void)
{
    int i;
    Tfill_buff_ibd1 =0;
    Tfill_buff_ibd2 =0;
    Tstore_buff_ibd = 0;
    Tcal_kinship = 0;
    Tgzread = 0;
    Tgzputs = 0;
    Tshrinke_id_pair = 0;
    Tcal_pair_num = 0;
    t_fill_buff_ibd1 =0;
    t_fill_buff_ibd2 =0;
    t_store_buff_ibd = 0;
    t_cal_kinship = 0;
    t_gzread = 0;
    t_gzputs = 0;
    t_shrinke_id_pair = 0;
    t_cal_pair_num = 0;

    haptag[0] = '1'; haptag[1] = '2';

    fprintf(stderr, "allocate IBD buffer\n");
    IBDdat = (IBD_t **)calloc(BUFF_row, sizeof(IBD_t*));
    strbuff = (char **)calloc(BUFF_row, sizeof(char*));
    strbuff2 = (char **)calloc(BUFF_row, sizeof(char*));
    for(i = 0; i < BUFF_row ;i ++){
	IBDdat[i] = (IBD_t *)calloc(1, sizeof(IBD_t));
	strbuff[i] = (char *)calloc(BUFF_col, sizeof(char));
	strbuff2[i] = (char *)calloc(BUFF_col, sizeof(char));
	memSize += sizeof(IBD_t) + 2 * BUFF_col * sizeof(char);
    }
    memSize += BUFF_row * (sizeof(IBD_t*) + 2 * sizeof(char*));
    return ;
}


