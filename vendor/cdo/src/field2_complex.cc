/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cdi.h>

#include "arithmetic.h"
#include "field_functions.h"
#include "cdo_output.h"

static void
field2_add_complex(Field &field1, Field const &field2)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  if (field1.nwpv != 2) cdo_abort("Field1 is not complex!");
  if (field2.nwpv != 2) cdo_abort("Field2 is not complex!");

  auto gridsize = gridInqSize(field1.grid);
  if (gridsize != gridInqSize(field2.grid)) cdo_abort("Fields have different size (%s)", __func__);

  for (size_t i = 0; i < gridsize; ++i)
  {
    array1[2 * i] = ADDM(array1[2 * i], array2[2 * i]);
    array1[2 * i + 1] = ADDM(array1[2 * i + 1], array2[2 * i + 1]);
  }
}

static void
field2_sub_complex(Field &field1, Field const &field2)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  if (field1.nwpv != 2) cdo_abort("Field1 is not complex!");
  if (field2.nwpv != 2) cdo_abort("Field2 is not complex!");

  auto gridsize = gridInqSize(field1.grid);
  if (gridsize != gridInqSize(field2.grid)) cdo_abort("Fields have different size (%s)", __func__);

  for (size_t i = 0; i < gridsize; ++i)
  {
    array1[2 * i] = SUBM(array1[2 * i], array2[2 * i]);
    array1[2 * i + 1] = SUBM(array1[2 * i + 1], array2[2 * i + 1]);
  }
}

static void
field2_mul_complex(Field &field1, Field const &field2)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  if (field1.nwpv != 2) cdo_abort("Field1 is not complex!");
  if (field2.nwpv != 2) cdo_abort("Field2 is not complex!");

  auto gridsize = gridInqSize(field1.grid);
  if (gridsize != gridInqSize(field2.grid)) cdo_abort("Fields have different size (%s)", __func__);

  // z1 x z2 = (x1x2 - y1y2) + i(x1y2 + x2y1)
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto a1r = array1[2 * i];
    auto a1i = array1[2 * i + 1];
    array1[2 * i] = SUBM(MULM(a1r, array2[2 * i]), MULM(a1i, array2[2 * i + 1]));
    array1[2 * i + 1] = ADDM(MULM(a1r, array2[2 * i + 1]), MULM(a1i, array2[2 * i]));
  }
}

static void
field2_div_complex(Field &field1, Field const &field2)
{
  auto is_EQ = fp_is_equal;
  auto missval1 = field1.missval;
  auto missval2 = field2.missval;
  auto &array1 = field1.vec_d;
  auto const &array2 = field2.vec_d;

  if (field1.nwpv != 2) cdo_abort("Field1 is not complex!");
  if (field2.nwpv != 2) cdo_abort("Field2 is not complex!");

  auto gridsize = gridInqSize(field1.grid);
  if (gridsize != gridInqSize(field2.grid)) cdo_abort("Fields have different size (%s)", __func__);

  // z1 / z2 = (x1x2 + y1y2) / (x2x2 + y2y2) + i (y1x2 - x1y2) / (x2x2 + y2y2)
  for (size_t i = 0; i < gridsize; ++i)
  {
    auto a1r = array1[2 * i];
    auto a1i = array1[2 * i + 1];
    auto denominator = ADDM(MULM(array2[2 * i], array2[2 * i]), MULM(array2[2 * i + 1], array2[2 * i + 1]));
    array1[2 * i] = DIVM(ADDM(MULM(a1r, array2[2 * i]), MULM(a1i, array2[2 * i + 1])), denominator);
    array1[2 * i + 1] = DIVM(SUBM(MULM(a1i, array2[2 * i]), MULM(a1r, array2[2 * i + 1])), denominator);
  }
}

void
field2_function_complex(Field &field1, Field const &field2, int function)
{
  // clang-format off
  switch (function)
    {
    case FieldFunc_Add:     field2_add_complex(field1, field2);   break;
    case FieldFunc_Sub:     field2_sub_complex(field1, field2);   break;
    case FieldFunc_Mul:     field2_mul_complex(field1, field2);   break;
    case FieldFunc_Div:     field2_div_complex(field1, field2);   break;
    default: cdo_abort("%s: function %d not implemented!", __func__, function);
    }
  // clang-format on
}
