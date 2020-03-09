#include "head.h"
#include "read.h"
#include "tools.h"
#include "print.h"
#include "parallel.h"


int main(int argc, char **argv)
{
    double a = my_wallclock();
    memSize = 0;

    //fprintf(stderr, "ID_t-%lu hapIBD_t-%lu IBD_t-%lu Pair_t-%lu void*-%lu long int-%lu int-%lu\n", sizeof(ID_t), sizeof(hapIBD_t), sizeof(IBD_t), sizeof(Pair_t), sizeof(void*), sizeof(long int), sizeof(int));

    check_input(argc, argv);
    p_parameters();
    init();    
    p_mem();

    read_ind(idfile);
    p_mem();

    read_map(mapfile);
    p_mem();

    fprintf(stdout, "\n\nFirst pass ...\n");
    read_ibd_1st_pass(); //first pass, index samples and save pairs
    cal_pair_num();
    p_mem();


    if(tagMask == 1)p_mask();
    if(tagCoverage == 1)p_coverage();

    if(tagKinship == 1){
	fprintf(stdout, "\n\nShrink ID pairs...\n");
	shrink_id_pair(); //first pass, index samples and save pairs
	cal_pair_num();
	p_mem();

	fprintf(stdout, "\n\nSecond pass ...\n");
	read_ibd_2nd_pass(); 
    }
    if(checkTime == 1)p_time();

    p_std_time(my_wallclock() - a);
    fprintf(stdout,"\n\nDone, Wallclock Time = %s\n", timestr);
    free_all();
    return 0;
}
