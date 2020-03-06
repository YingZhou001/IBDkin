#include "head.h"
#include "print.h"
#include "read.h"
#include "tools.h"

void print_help(void)
{
    fprintf(stderr, "Usage: IBDkin [options] parameters\n\n \
    (Required inputs:)\n \
    \t-ibdfile [file]\t#<string> a list of input IBD pathnames\n \
    \t-map [file]\t#<string> genetic map with cM distance in plink format\n \
    \t-ind [file]\t#<string> a list individuals to be analyzed\n\n \
    (Optional parameters:)\n \
    \t-out ./\t\t#<string> output prefix\n \
    \t-nthreads 2 \t#<int> number of threads\n \
    \t-degree 9 \t#<int> max relationship degree\n \
    \t-binkb 1000 \t#<float> bin size in kbp to calculate IBD coverage\n \
    \t-fold 4 \t#<float> max fold deviation to determine masked regions\n \
    \t-part 1 1 \t#<int> <int> total partitions and current partition\n \
    \t-cutcm 4 2 \t#<float> <float> minimum long and short IBD segment lengths in cM \n \
    \t-merge 5 20\t#<float> <float> max cM merge lengths for IBD1 and IBD2 regions \n\n \
    (Other flags:)\n \
    \t--nokinship\t#do not output kinship coefficients\n \
    \t--outmask\t#output genome mask\n \
    \t--outcoverage\t#output IBD coverage\n \
    \n");
    exit(-1);
}


void p_parameters(void)
{
    int i;
    fprintf(stderr, "###Parameters###\n");
    fprintf(stderr, "Input:\n");
    fprintf(stderr, " ibd: *%s*\n", headfile);
    read_headfile(headfile);
    fprintf(stderr, " map: *%s*\n", mapfile);
    fprintf(stderr, " ind: *%s*\n", idfile);

    fprintf(stderr, "\nOutputs:\n");
    if(tagKinship == 1){
	fprintf(stderr, "  *%s*\n", outfile);
	gzFile fp = gzopen(outfile, "w");
	gzbuffer(fp, BUFF3);
	gzprintf(fp, "ID1 ID2 segnum ibd0 ibd1 ibd2 kinship degree\n");
	gzclose(fp);
    }
    if(tagMask == 1)fprintf(stderr, "  *%s*\n", maskfile);
    if(tagCoverage == 1)fprintf(stderr, "  *%s*\n", coveragefile);

    fprintf(stderr, "\nOther settings:\n");
    fprintf(stderr, " *ibd cutoff-1 = %fcM*\n", IBDcM);
    fprintf(stderr, " *ibd cutoff-2 = %fcM*\n", IBDcM2);
    fprintf(stderr, " *degree = %d, kinship = %f*\n",degree, kincut);
    fprintf(stderr, " *fold = %f*\n", FOLD);
    fprintf(stderr, " *binsize = %.3fkb*\n", BINSIZE);
    fprintf(stderr, " *nthreads = %d*\n", Nthreads);

    threadCount = (long int *)calloc(Nthreads, sizeof(long int));
    memSizeByThread = (long int *)calloc(Nthreads, sizeof(long int));
    pairNumByThread = (long int *)calloc(Nthreads, sizeof(long int));
    segNumByThread = (long int *)calloc(Nthreads, sizeof(long int));
    outbuffi = (int *)calloc(Nthreads, sizeof(int));
    outbuff = (char **)calloc(Nthreads, sizeof(char *));

    memSize += 4 * Nthreads * sizeof(long int) + Nthreads * (sizeof(int) + sizeof(char *));

    for(i = 0; i< Nthreads; i++) {
	threadCount[i] = 0;
	memSizeByThread[i] = 0;
	pairNumByThread[i] = 0;
	segNumByThread[i] = 0;
	outbuffi[i] = 0;
	outbuff[i] = (char *)calloc(BUFF3, sizeof(char));
	memSize += BUFF3 * sizeof(char);
    }

    if(gap1 > 0){
	fprintf(stderr, " *merge threshold for IBD1 = %fcM*\n", gap1);
	fprintf(stderr, " *merge threshold for IBD2 = %fcM*\n", gap2);
    }
    else fprintf(stderr, " **no merge**\n");
    if(Parts > 1)fprintf(stderr, " **%d partitions in total, this is the %d's run**\n", Parts, part);
    else fprintf(stderr, " **no partitions**\n");


    fprintf(stderr, "#################\n\n");
}



void p_coverage(void)
{
    double a = my_wallclock();
    long int b = clock();
    gzFile fp = gzopen(coveragefile, "w");
    int chr, i;
    char str[1024];
    float p1, p2, g1, g2;
    float  bin = BINSIZE * 1000;
    for(chr = 1; chr < 23; chr ++){
	if( minPos[chr] < maxPos[chr]){
	    for(i = 0; i < coverageL[chr]; i++){
		p1 = i*bin; 
		p2 = p1 + bin;
		g1 = interpolate2cM(chr, p1);
		g2 = interpolate2cM(chr, p2);
		if(g2 >= minPos[chr] && g1 <= maxPos[chr]){
		    //sprintf(str, "%d %.0f %.0f %f %f %f\n", chr, p1, p2, g1, g2, coverage[chr][i]);
		    sprintf(str, "%d %.0f %.0f %f\n", chr, p1, p2, coverage[chr][i]);
		    gzputs(fp, str);
		}
	    }
	}
    }
    gzclose(fp);
    Tgzputs += clock() - b;
    t_gzputs += my_wallclock() - a;
    return ;
}

void p_mask(void)
{
    double a = my_wallclock();
    long int b = clock();
    gzFile fp = gzopen(maskfile, "w");
    char str[BUFF];
    int chr;
    msk_t *cur;
    for(chr = 1; chr < 23; chr ++){
	cur = mask[chr];
	while(cur){
	    if(cur->p1 > 0){
		//sprintf(str,"%d:%.0f %.0f %f %f\n", chr, cur->p1, cur->p2, cur->g1, cur->g2);
		sprintf(str,"%d %.0f %.0f\n", chr, cur->p1, cur->p2);
		gzputs(fp, str);
	    }
	    cur = cur -> next;
	}
    }
    gzclose(fp);
    Tgzputs += clock() - b;
    t_gzputs += my_wallclock() - a;
    return ;
}


void p_mem(void)
{
    if(checkMem != 1)return;
    p_maxrss();
    double mem;
    long int tmp;
    int i;
    tmp = memSize; 
    for(i = 0; i < Nthreads; i++ )tmp += memSizeByThread[i];
    mem = (double)tmp / 1073741824;
    if(mem > 1){fprintf(stderr, " *Memory in use: %.5lf Gigabytes\n\n", mem);return;}
    mem = (double)tmp / 1048576;
    if(mem > 1){fprintf(stderr, " *Memory in use: %.5lf Megabytes\n\n", mem);return;}
    mem = (double)tmp / 1024;
    if(mem > 1){fprintf(stderr, " *Memory in use: %.5lf Kilobytes\n\n", mem);return;}
    fprintf(stderr, " *Memory in use: %ld Bytes\n\n", tmp);

    return ;

}

void p_maxrss(void)
{
    struct rusage buf;
    getrusage(RUSAGE_SELF, &buf);
    double mem, tmp = buf.ru_maxrss;
    mem = (double)tmp / 1048576;
    if(mem > 1){fprintf(stderr, " *maxrss: %.5lf Gigabytes\n", mem);return;}
    mem = (double)tmp / 1024;
    if(mem > 1){fprintf(stderr, " *maxrss: %.5lf Megabytes\n", mem);return;}
    fprintf(stderr, " *maxrss: %.5lf Kilobytes\n", tmp);
    return ;
}


void p_time(void)
{
    fprintf(stderr,"\n function time:\n");
    p_std_time(t_gzread);
    fprintf(stderr," *gz_read() cpu=%.3fs, wall=%s, r=%.3f\n", Tgzread*1.0/CLOCKS_PER_SEC, timestr, Tgzread*1.0/CLOCKS_PER_SEC/t_gzread);
    p_std_time(t_fill_buff_ibd1);
    fprintf(stderr," *fill_buff_ibd1() cpu=%.3fs, wall=%s, r=%.3f\n", Tfill_buff_ibd1*1.0/CLOCKS_PER_SEC, timestr, Tfill_buff_ibd1*1.0/CLOCKS_PER_SEC/t_fill_buff_ibd1);
    if(tagKinship == 1){
	p_std_time(t_fill_buff_ibd2);
	fprintf(stderr," *fill_buff_ibd2() cpu=%.3fs, wall=%s, r=%.3f\n", Tfill_buff_ibd2*1.0/CLOCKS_PER_SEC, timestr, Tfill_buff_ibd2*1.0/CLOCKS_PER_SEC/t_fill_buff_ibd2);
    }
    if(tagKinship == 1){
	p_std_time(t_shrinke_id_pair);
	fprintf(stderr," *shrinke_id_pair() cpu=%.3fs, wall=%s, r=%.3f\n", Tshrinke_id_pair*1.0/CLOCKS_PER_SEC, timestr, Tshrinke_id_pair*1.0/CLOCKS_PER_SEC/t_shrinke_id_pair);
    }
    p_std_time(t_cal_pair_num);
    fprintf(stderr," *cal_pair_num() cpu=%.3fs, wall=%s, r=%.3f\n", Tcal_pair_num*1.0/CLOCKS_PER_SEC, timestr, Tcal_pair_num*1.0/CLOCKS_PER_SEC/t_cal_pair_num);

    p_std_time(t_store_buff_ibd);
    fprintf(stderr," *store_buff_ibd() cpu=%.3fs, wall=%s, r=%.3f\n", Tstore_buff_ibd*1.0/CLOCKS_PER_SEC, timestr, Tstore_buff_ibd*1.0/CLOCKS_PER_SEC/t_store_buff_ibd);
    if(tagKinship == 1){
	p_std_time(t_cal_kinship);
	fprintf(stderr," *cal_kinship() cpu=%.3fs, wall=%s, r=%.3f\n", Tcal_kinship*1.0/CLOCKS_PER_SEC, timestr, Tcal_kinship*1.0/CLOCKS_PER_SEC/t_cal_kinship);
    }
    p_std_time(t_gzputs);
    fprintf(stderr," *gz_write() cpu=%.3fs, wall=%s, r=%.3f\n", Tgzputs*1.0/CLOCKS_PER_SEC, timestr, Tgzputs*1.0/CLOCKS_PER_SEC/t_gzputs);
}

void p_std_time(double seconds)
{
    int hour, min;
    double sec;
    hour = (int) (seconds / 3600);
    min = (int) ((seconds - hour * 3600)/ 60);
    sec = seconds - hour * 3600 - min * 60;
    if(hour > 0)sprintf(timestr,"%dh%dm%.3fs", hour, min, sec);
    else if (min > 0)sprintf(timestr,"%dm%.3fs", min, sec);
    else sprintf(timestr,"%.3fs", sec);
    return ;
}
