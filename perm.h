typedef enum { AUTO, LINEAR_INTERPOLATION, NEAREST_RANK } cp_mode;
typedef enum { CASE, CONTROL } group;

struct perm_data 
{
  int n_var;
  int n_cases;
  int n_controls;
  group *grp;
  float *data;
  float *pc_cases;
  float *pc_controls;
  float *pc_diff;
  float *fold_change;
  double *p_values;
};

int check_percentile(double, cp_mode, int, struct perm_data *);
