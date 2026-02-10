/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/
#ifndef FIELD_FUNCTIONS_H
#define FIELD_FUNCTIONS_H

#include "field.h"

double var_to_std(double rvar, double missval);

enum Func
{
  Func_Name,
  Func_Code,
  Func_Param
};

enum FieldFunc
{
  FieldFunc_Min = 100,
  FieldFunc_Max,
  FieldFunc_Range,
  FieldFunc_Sum,
  FieldFunc_Avg,
  FieldFunc_Mean,
  FieldFunc_Var,
  FieldFunc_Var1,
  FieldFunc_Std,
  FieldFunc_Std1,
  FieldFunc_Skew,
  FieldFunc_Kurt,
  FieldFunc_Median,
  FieldFunc_Count,
  FieldFunc_Pctl,

  FieldFunc_Cor,
  FieldFunc_Covar,
  FieldFunc_Avgw,
  FieldFunc_Meanw,
  FieldFunc_Stdw,
  FieldFunc_Std1w,
  FieldFunc_Varw,
  FieldFunc_Var1w,
  FieldFunc_Minidx,
  FieldFunc_Maxidx,
  FieldFunc_Rmsd,

  FieldFunc_Add,
  FieldFunc_Sub,
  FieldFunc_Mul,
  FieldFunc_Div,
  FieldFunc_Mod,

  FieldFunc_EQ,
  FieldFunc_NE,
  FieldFunc_LE,
  FieldFunc_LT,
  FieldFunc_GE,
  FieldFunc_GT,

  FieldFunc_Atan2,
  FieldFunc_Setmiss,
};

// field_memory.cc
void field2D_init(FieldVector2D &field2D, VarList const &varList);
void field2D_init(FieldVector2D &field2D, VarList const &varList, int ptype);
void field2D_init(FieldVector2D &field2D, VarList const &varList, int ptype, double fillValue);
void field1Dvars_init(FieldVector &field1D, VarList const &varList);
void field1Dvars_init(FieldVector &field1D, VarList const &varList, int ptype);
void field1Dlevels_init(FieldVector &field1D, VarList const &varList);
void field1Dlevels_init(FieldVector &field1D, VarList const &varList, int ptype);

// field.cc
double field_function(Field const &field, int function);

double field_min(Field const &field);
double field_range(Field const &field);
double field_max(Field const &field);
double field_sum(Field const &field);
double field_mean(Field const &field);
double field_meanw(Field const &field);
double field_avg(Field const &field);
double field_avgw(Field const &field);
double field_std(Field const &field);
double field_std1(Field const &field);
double field_var(Field const &field);
double field_var1(Field const &field);
double field_stdw(Field const &field);
double field_std1w(Field const &field);
double field_varw(Field const &field);
double field_var1w(Field const &field);
double field_skew(Field const &field);
double field_kurt(Field const &field);
double field_median(Field const &field);
double field_count(Field const &field);

// ENS validation
double field_rank(Field &field);

double field_pctl(Field &field, double pn);

// field_zonal.cc
void zonal_function(Field const &field1, Field &field2, int function);
void zonal_min(Field const &field1, Field &field2);
void zonal_max(Field const &field1, Field &field2);
void zonal_range(Field const &field1, Field &field2);
void zonal_sum(Field const &field1, Field &field2);
void zonal_avg(Field const &field1, Field &field2);
void zonal_mean(Field const &field1, Field &field2);
void zonal_std(Field const &field1, Field &field2);
void zonal_std1(Field const &field1, Field &field2);
void zonal_var(Field const &field1, Field &field2);
void zonal_var1(Field const &field1, Field &field2);
void zonal_skew(Field const &field1, Field &field2);
void zonal_kurt(Field const &field1, Field &field2);
void zonal_median(Field const &field1, Field &field2);
void zonal_pctl(Field const &field1, Field &field2, double pn);

// field_meridional.cc
void meridional_function(Field const &field1, Field &field2, int function);
void meridional_pctl(Field const &field1, Field &field2, double pn);

void field_rms(Field const &field1, Field const &field2, Field &field3);

// fieldc.cc
void fieldc_function(Field &field, double rconst, int function);

void fieldc_mul(Field &field, double rconst);
void fieldc_div(Field &field, double rconst);
void fieldc_add(Field &field, double rconst);
void fieldc_sub(Field &field, double rconst);
void fieldc_min(Field &field, double rconst);
void fieldc_max(Field &field, double rconst);
void fieldc_mod(Field &field, double divisor);

// fieldc_complex.cc
void fieldc_function_complex(Field &field, const double rconstcplx[2], int function);

// field2.cc
void field2_function(Field &field1, Field const &field2, int function);

void field2_add(Field &field1, Field const &field2);
void field2_sum(Field &field1, Field const &field2);
void field2_sub(Field &field1, Field const &field2);
void field2_mul(Field &field1, Field const &field2);
void field2_div(Field &field1, Field const &field2);
void field2_min(Field &field1, Field const &field2);
void field2_max(Field &field1, Field const &field2);
void field2_atan2(Field &field1, Field const &field2);

void field2_sumq(Field &field1, Field const &field2);
void field2_sumw(Field &field1, Field const &field2, double w);
void field2_sumqw(Field &field1, Field const &field2, double w);
void field2_sumtr(Field &field1, Field const &field2, double refval);
void field2_vinit(Field &field1, Field const &field2);
void field2_vincr(Field &field1, Field const &field2);
void field2_vinit(Field &field1, Field const &field2, int vinit);
void field2_vincr(Field &field1, Field const &field2, int vincr);

void field2_sumsumq(Field &field1, Field &field2, Field const &field3);
void field2_maxmin(Field &field1, Field &field2, Field const &field3);
void field2_minidx(Field &field1, Field &field2, Field const &field3, int idx);
void field2_maxidx(Field &field1, Field &field2, Field const &field3, int idx);
void field2_var(Field &field1, Field const &field2, Field const &field3, int divisor);
void field2_std(Field &field1, Field const &field2, Field const &field3, int divisor);
void fieldc_var(Field &field1, Field const &field2, int numSets, int divisor);
void fieldc_std(Field &field1, Field const &field2, int numSets, int divisor);
void field2_moq(Field &field1, Field const &field2);
void field2_moqw(Field &field1, Field const &field2, double w);

void field2_count(Field &field1, Field const &field2);

// field2_complex.cc
void field2_function_complex(Field &field1, Field const &field2, int function);

#endif /* FIELD_FUNCTIONS_H */
