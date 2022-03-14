#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
#include<limits.h>
#include<time.h>
#include<pthread.h>
#include<zlib.h>
#include<assert.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <errno.h> 
#include <stdarg.h>


#define BUFF_id 32 //for ID/Chr string length
#define BUFF_col 2048 //should larger than BUFF_id*2+1024
#define BUFF_row 1048576 //BUFF_col * BUFF_row * 2  determine mem in use 
#define BUFF_gz 1048576 //buffer size for gzip 
#define BUFF_filename 4096 //string that define the file name

long int Tfill_buff_ibd1, Tfill_buff_ibd2, Tstore_buff_ibd, Tcal_kinship, Tgzputs, Tgzread, Tshrinke_id_pair, Tcal_pair_num, Tcpbuff; //timing
double t_fill_buff_ibd1, t_fill_buff_ibd2, t_store_buff_ibd, t_cal_kinship, t_gzputs, t_gzread, t_shrinke_id_pair, t_cal_pair_num, t_cpbuff; //timing

int checkTime, checkMem;
int tagKinship, tagMask, tagCoverage;

int Parts, part;//data would be devide into ${Parts} and output the ${part} in this run
char *headfile, **ibdfile;
char *idfile, *mapfile, *outfile;
char *maskfile, *coveragefile, *rangefile;
char **chrs;
int Nfile;
int buffi, buffi2;
int KEYnumPAIR;
int KEYnumID;
int Nsample;
int idNum;
int chrN, CHR; //working chromosome number


long int *threadCount;
long int memSize; // in bytes
long int *memSizeByThread;
long int pairNum, segNum;
long int *pairNumByThread, *segNumByThread;

float minPos[100], maxPos[100]; //record minimum and maximum genetic positions for each chromosome. 22 for 22 automosomes
float **minpos, **maxpos; 
float BINSIZE; //kb
double **coverage, median;

int Nthreads;
int mapL[100];
int coverageL[100];
//int *outbuffi;
int Round;

float IBDcM, IBDcM2;
float **mapP;
float **mapG;
float totg;
float kincut;
float FOLD;
float gap1, gap2;

char haptag[2];
char *lock;
char **strbuff;
char **strbuff2;
//char **outbuff;
char timestr[BUFF_col];

gzFile ifp;
FILE *ofp;

typedef struct Pair
{
    int index;
    int segnum;
    float s1, s2;
    struct hapIBD *ibdhead;
} Pair_t;

/*save ID*/
typedef struct ID
{
    char *id;
    int num;//num of pairs stored
    int *stack;//stack[0] is the capacity, this array will be freed in the second pass
    struct Pair **head;
} ID_t;

ID_t **idhead;

typedef struct IBD //BUFFER
{
    int id1, id2;
    char chr;
    char pass;
    char h1, h2; //used to distinguish IBD1 and IBD2
    float g1, g2, l; //genetic positions
} IBD_t;

IBD_t **IBDdat;

typedef struct hapIBD //store ibs information
{
    char h1, h2; //used to distinguish IBD1 and IBD2
    float g1, g2; //genetic positions
    struct hapIBD *next;
} hapIBD_t;

typedef struct cvrg
{
   double **coverage;
} cvrg_t;

cvrg_t **coverage_per_thread;


typedef struct msk
{
    float p1, p2;
    float g1, g2;
    struct msk *next;
} msk_t;

msk_t **mask;
