int hash_str(char *id);

int ID_t_cmp(const void *a, const void *b);
int int_cmp(const void *a, const void *b);
int float_cmp(const void *a, const void *b);
int double_cmp(const void *a, const void *b);
int string_cmp(const void *p1, const void *p2);

int binary_search(float pos, float *map, int mapl);
int binary_search_string(char **S, int L, char *s);
float interpolate2cM(int chr, float p);

int binary_search_ID_index(char *id);
int binary_search_pair(ID_t *id, int index);

void init_lock(void);

void cal_coverage(int th_index, int chr, float p1, float p2);
void cal_coverage_median(void);
size_t push_mask(msk_t **ref, int p1, int p2, float g1, float g2);
void create_mask(void);
void delete_mask(msk_t *head);
float apply_mask(int chr, float g1, float g2);

int cal_degree(float kinship);
int check_chr(void);

void free_all(void);
long int freehapIBDList(hapIBD_t *head);

double my_wallclock(void);

size_t push_ibd(hapIBD_t **head, char h1, char h2, float g1, float g2);
void my_error(const char * message, ...);
