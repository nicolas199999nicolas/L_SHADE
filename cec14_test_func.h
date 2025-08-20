#ifndef CEC14_TEST_FUNC_H
#define CEC14_TEST_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif

void cec14_test_func(double *x, double *f, int nx, int mx, int func_num);

void sphere_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void ellips_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void bent_cigar_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void discus_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void dif_powers_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void rosenbrock_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void schaffer_F7_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void ackley_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void weierstrass_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void griewank_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void schwefel_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void katsuura_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void bi_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void grie_rosen_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void escaffer6_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void step_rastrigin_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void happycat_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);
void hgbat_func(double *x, double *f, int nx, double *Os, double *Mr, int s_flag, int r_flag);

void hf01(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);
void hf02(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);
void hf03(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);
void hf04(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);
void hf05(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);
void hf06(double *x, double *f, int nx, double *Os, double *Mr, int *S, int s_flag, int r_flag);

void cf01(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf02(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf03(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf04(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf05(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf06(double *x, double *f, int nx, double *Os, double *Mr, int r_flag);
void cf07(double *x, double *f, int nx, double *Os, double *Mr, int *SS, int r_flag);
void cf08(double *x, double *f, int nx, double *Os, double *Mr, int *SS, int r_flag);

#ifdef __cplusplus
}
#endif

#endif // CEC14_TEST_FUNC_H