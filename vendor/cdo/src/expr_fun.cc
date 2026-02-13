/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include "cdo_output.h"
#include "cdo_zaxis.h"
#include "cdo_omp.h"
#include <mpim_grid.h>
#include "expr.h"
#include "expr_fun.h"
#include "expr_yacc.hh"

// clang-format off
auto expr_func_con_var = [](auto hasMV, auto n, auto mv, auto cVal, auto const &vIn, auto &vOut, auto binary_operator)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = fp_is_equal(vIn[i], mv) ? mv : binary_operator(cVal, vIn[i]); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(vIn[i], mv) ? mv : binary_operator(cVal, vIn[i]); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_operator(cVal, vIn[i]); }
    }
};

auto expr_mul_con_var = [](auto hasMV, auto n, auto mv, auto cVal, auto const &vIn, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(cVal, 0.0) ? 0.0 : fp_is_equal(vIn[i], mv) ? mv : binary_op_MUL(cVal, vIn[i]); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(cVal, 0.0) ? 0.0 : is_equal(vIn[i], mv) ? mv : binary_op_MUL(cVal, vIn[i]); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_op_MUL(cVal, vIn[i]); }
    }
};

auto expr_div_con_var = [](auto hasMV, auto n, auto mv, auto cVal, auto const &vIn, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (fp_is_equal(vIn[i], mv) || fp_is_equal(vIn[i], 0.0)) ? mv : binary_op_DIV(cVal, vIn[i]); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (is_equal(vIn[i], mv) || is_equal(vIn[i], 0.0)) ? mv : binary_op_DIV(cVal, vIn[i]); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(vIn[i], 0.0) ? mv : binary_op_DIV(cVal, vIn[i]); }
    }
};

auto expr_func_var_con = [](auto hasMV, auto n, auto mv, auto const &vIn, auto cVal, auto &vOut, auto binary_operator)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = fp_is_equal(vIn[i], mv) ? mv : binary_operator(vIn[i], cVal); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(vIn[i], mv) ? mv : binary_operator(vIn[i], cVal); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_operator(vIn[i], cVal); }
    }
};

auto expr_mul_var_con = [](auto hasMV, auto n, auto mv, auto const &vIn, auto cVal, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(cVal, 0.0) ? 0.0 : fp_is_equal(vIn[i], mv) ? mv : binary_op_MUL(vIn[i], cVal); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(cVal, 0.0) ? 0.0 : is_equal(vIn[i], mv) ? mv : binary_op_MUL(vIn[i], cVal); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_op_MUL(vIn[i], cVal); }
    }
};

auto expr_div_var_con = [](auto hasMV, auto n, auto mv, auto const &vIn, auto cVal, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (fp_is_equal(vIn[i], mv) || fp_is_equal(cVal, 0.0)) ? mv : binary_op_DIV(vIn[i], cVal); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (is_equal(vIn[i], mv) || is_equal(cVal, 0.0)) ? mv : binary_op_DIV(vIn[i], cVal); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(cVal, 0.0) ? mv : binary_op_DIV(vIn[i], cVal); }
    }
};

auto expr_func_var_var = [](auto hasMV, auto n, auto mv1, auto mv2, auto const &vIn1, auto const &vIn2, auto &vOut, auto binary_operator)
{
  if (hasMV)
    {
      if (std::isnan(mv1) || std::isnan(mv2))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (fp_is_equal(vIn1[i], mv1) || fp_is_equal(vIn2[i], mv2)) ? mv1 : binary_operator(vIn1[i], vIn2[i]); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (is_equal(vIn1[i], mv1) || is_equal(vIn2[i], mv2)) ? mv1 : binary_operator(vIn1[i], vIn2[i]); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_operator(vIn1[i], vIn2[i]); }
    }
};

auto expr_mul_var_var = [](auto hasMV, auto n, auto mv1, auto mv2, auto const &vIn1, auto const &vIn2, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv1) || std::isnan(mv2))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i)
          {
            vOut[i] = (fp_is_equal(vIn1[i], 0.0) || fp_is_equal(vIn2[i], 0.0)) ? 0.0 : (fp_is_equal(vIn1[i], mv1) || fp_is_equal(vIn2[i], mv2)) ? mv1 : binary_op_MUL(vIn1[i], vIn2[i]);
          }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i)
          {
            vOut[i] = (is_equal(vIn1[i], 0.0) || is_equal(vIn2[i], 0.0)) ? 0.0 : (is_equal(vIn1[i], mv1) || is_equal(vIn2[i], mv2)) ? mv1 : binary_op_MUL(vIn1[i], vIn2[i]);
          }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = binary_op_MUL(vIn1[i], vIn2[i]); }
    }
};

auto expr_div_var_var = [](auto hasMV, auto n, auto mv1, auto mv2, auto const &vIn1, auto const &vIn2, auto &vOut)
{
  if (hasMV)
    {
      if (std::isnan(mv1) || std::isnan(mv2))
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (fp_is_equal(vIn1[i], mv1) || fp_is_equal(vIn2[i], mv2) || fp_is_equal(vIn2[i], 0.0)) ? mv1 : binary_op_DIV(vIn1[i], vIn2[i]); }
      else
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
        for (size_t i = 0; i < n; ++i) { vOut[i] = (is_equal(vIn1[i], mv1) || is_equal(vIn2[i], mv2) || is_equal(vIn2[i], 0.0)) ? mv1 : binary_op_DIV(vIn1[i], vIn2[i]); }
    }
  else
    {
#ifdef _OPENMP
#pragma omp parallel for if (n > cdoMinLoopSize) default(shared)
#endif
      for (size_t i = 0; i < n; ++i) { vOut[i] = is_equal(vIn2[i], 0.0) ? mv1 : binary_op_DIV(vIn1[i], vIn2[i]); }
    }
};
// clang-format on

nodeType *
expr_con_con(int oper, const nodeType *p1, const nodeType *p2)
{
  auto p = new nodeType;

  p->type = NodeEnum::typeCon;
  p->isTmpObj = true;

  auto cval1 = p1->con().value;
  auto cval2 = p2->con().value;

  // clang-format off
  switch (oper)
    {
    case LT:   cval1 = static_cast<double>(binary_op_LT(cval1, cval2)); break;
    case GT:   cval1 = static_cast<double>(binary_op_GT(cval1, cval2)); break;
    case LE:   cval1 = static_cast<double>(binary_op_LE(cval1, cval2)); break;
    case GE:   cval1 = static_cast<double>(binary_op_GE(cval1, cval2)); break;
    case NE:   cval1 = static_cast<double>(binary_op_NE(cval1, cval2)); break;
    case EQ:   cval1 = static_cast<double>(binary_op_EQ(cval1, cval2)); break;
    case LEG:  cval1 = static_cast<double>(binary_op_LEG(cval1, cval2)); break;
    case AND:  cval1 = static_cast<double>(binary_op_AND(cval1, cval2)); break;
    case OR:   cval1 = static_cast<double>(binary_op_OR(cval1, cval2)); break;
    case '^':  cval1 = static_cast<double>(binary_op_POW(cval1, cval2)); break;
    case '+':  cval1 = static_cast<double>(binary_op_ADD(cval1, cval2)); break;
    case '-':  cval1 = static_cast<double>(binary_op_SUB(cval1, cval2)); break;
    case '*':  cval1 = static_cast<double>(binary_op_MUL(cval1, cval2)); break;
    case '/':  cval1 = static_cast<double>(binary_op_DIV(cval1, cval2)); break;
    default:   cdo_abort("%s: operator %d unsupported!", __func__, oper); break;
    }
  // clang-format on

  p->v = conNodeType(cval1);

  return p;
}

void
oper_expr_con_var(int oper, bool hasMV, size_t n, double mv, double *odat, double cval, const double *idat)
{
  // clang-format off
  switch (oper)
    {
    case LT:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_LT); break;
    case GT:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_GT); break;
    case LE:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_LE); break;
    case GE:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_GE); break;
    case NE:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_NE); break;
    case EQ:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_EQ); break;
    case LEG: expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_LEG); break;
    case AND: expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_AND); break;
    case OR:  expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_OR); break;
    case '^': expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_POW); break;
    case '+': expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_ADD); break;
    case '-': expr_func_con_var(hasMV, n, mv, cval, idat, odat, binary_op_SUB); break;
    case '*': expr_mul_con_var(hasMV, n, mv, cval, idat, odat); break;
    case '/': expr_div_con_var(hasMV, n, mv, cval, idat, odat); break;
    default: cdo_abort("%s: operator %d unsupported!", __func__, oper); break;
    }
  // clang-format on
}

void
oper_expr_var_con(int oper, bool hasMV, size_t n, double mv, double *odat, const double *idat, double cval)
{
  // clang-format off
  switch (oper)
    {
    case LT:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_LT); break;
    case GT:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_GT); break;
    case LE:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_LE); break;
    case GE:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_GE); break;
    case NE:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_NE); break;
    case EQ:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_EQ); break;
    case LEG: expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_LEG); break;
    case AND: expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_AND); break;
    case OR:  expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_OR); break;
    case '^': expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_POW); break;
    case '+': expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_ADD); break;
    case '-': expr_func_var_con(hasMV, n, mv, idat, cval, odat, binary_op_SUB); break;
    case '*': expr_mul_var_con(hasMV, n, mv, idat, cval, odat); break;
    case '/': expr_div_var_con(hasMV, n, mv, idat, cval, odat); break;
    default: cdo_abort("%s: operator %d unsupported!", __func__, oper); break;
    }
  // clang-format on
}

void
oper_expr_var_var(int oper, bool hasMV, size_t n, double mv1, double mv2, double *odat, const double *idat1, double *idat2)
{
  // clang-format off
  switch (oper)
    {
    case LT:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_LT); break;
    case GT:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_GT); break;
    case LE:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_LE); break;
    case GE:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_GE); break;
    case NE:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_NE); break;
    case EQ:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_EQ); break;
    case LEG: expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_LEG); break;
    case AND: expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_AND); break;
    case OR:  expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_OR); break;
    case '^': expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_POW); break;
    case '+': expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_ADD); break;
    case '-': expr_func_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat, binary_op_SUB); break;
    case '*': expr_mul_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat); break;
    case '/': expr_div_var_var(hasMV, n, mv1, mv2, idat1, idat2, odat); break;
    default: cdo_abort("%s: operator %d unsupported!", __func__, oper); break;
    }
  // clang-format on
}

void
fld_field_init(Field &field, size_t numMissVals, double missval, size_t ngp, double const *array, double const *w)
{
  field.size = ngp;
  field.numMissVals = numMissVals;
  field.missval = missval;
  if (array != nullptr) array_copy(ngp, array, field.vec_d.data());
  if (w != nullptr) array_copy(ngp, w, field.weightv.data());
}
/*
double *
fld_weights(int gridID, size_t ngp)
{
  static auto printWarning = true;
  double *weights = (double *) std::malloc(ngp * sizeof(double));
  for (size_t i = 0; i < ngp; ++i) weights[i] = 1;

  if (ngp > 1)
    {
      int wstatus = gridcell_weights(gridID, weights);
      if (wstatus != 0 && printWarning)
        {
          printWarning = false;
          cdo_warning("Grid cell bounds not available, using constant grid cell area weights!");
        }
    }

  return weights;
}
*/

void
vert_weights(int zaxisID, size_t nlev, Varray<double> &weights)
{
  weights.resize(nlev, 1.0);

  if (nlev > 1)
  {
    static auto printWarning = true;
    Varray<double> thickness(nlev);
    auto wstatus = get_layer_thickness(1, 0, 0, zaxisID, nlev, thickness, weights);
    if (wstatus == 0 && printWarning)
    {
      printWarning = false;
      cdo_warning("Layer bounds not available, using constant vertical weights!");
    }
  }
}
