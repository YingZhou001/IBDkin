void fill_buff_ibd_1st_pass(void);
void *fill_buff_ibd_1st_pass_by_thread(void *args);

void fill_buff_ibd_2nd_pass(void);
void *fill_buff_ibd_2nd_pass_by_thread(void *args);

void store_buff_ibd(void);
void *store_buff_ibd_by_thread(void *args);

int calIBD12_pair(hapIBD_t* head, float *S0, float *S1, float *S2);


//for remote relatives
void cal_kinship(void);
void *cal_kinship_by_thread(void *args);
void output(hapIBD_t *head, int idi, int ids, int buf);

//for close relatives
void cal_kinship2(void);
void *cal_kinship2_by_thread(void *args);
void output2(Pair_t *head);

void shrink_id_pair(void);
void *shrink_id_pair_by_thread(void *args);

void cal_pair_num(void);
int cal_seg_num(hapIBD_t *head);
void *cal_pair_num_by_thread(void *args);

void write_buff(void);
void read_buff(void);

void *copy_buff_by_thread(void *args);
void copy_buff(void);

