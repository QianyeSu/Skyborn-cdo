/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/
#ifndef ECA_H_
#define ECA_H_

#include "field.h"

enum ECA_EPILOG
{
  ECA_NONE,
  MEAN,
  PERCENT_OF_TIME,
  PERCENT_OF_TOTAL_AMOUNT
};

using ECA_FUNC_1 = void (*)(Field &, double);
using ECA_FUNC_2 = void (*)(Field &, Field const &);
using ECA_FUNC_3 = void (*)(Field &, Field const &, double);

/**
 * Structure describing a processing request of the form
 *
 * o = F3(F2(a * F1(i) + b))
 *
 * where i and o denote the input and output fields, and
 * F1, F2 and F3 are field operators.
 *
 * The structure contains the following elements:
 *
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f1arg     the argument of the 1st field operator
 * f2        the 2nd field operator
 * f3        the 3rd field operator
 * mulc      the multiplier a
 * addc      the addend b
 * epilog    the final operation carried out after processing
 */
struct ECA_MAJOR_REQUEST_ELEMENT_1
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  int refdate;
  ECA_FUNC_1 f1 = nullptr;
  double f1arg = 1.0;
  ECA_FUNC_2 f2 = nullptr;
  ECA_FUNC_2 f3 = nullptr;
  double mulc = 0.0;
  double addc = 0.0;
  ECA_EPILOG epilog = ECA_NONE;
};

/**
 * Structure describing a processing request of the form
 *
 * o = H3(H2(H1(i)))
 *
 * where i and o denote the input and output fields, and
 * H1, H2 and H3 are field operators.
 *
 * The structure contains the following elements:
 *
 * name      the name of the output variable
 * longname  the longname of the output variable
 * h1        the 1st field operator
 * h1arg     the argument of the 1st field operator
 * h2        the 2nd field operator
 * h3        the 3rd field operator
 */
struct ECA_MINOR_REQUEST_ELEMENT_1
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  ECA_FUNC_1 h1 = nullptr;
  double h1arg = 0.0;
  ECA_FUNC_2 h2 = nullptr;
  ECA_FUNC_2 h3 = nullptr;
};

struct ECA_REQUEST_1
{
  ECA_MAJOR_REQUEST_ELEMENT_1 var1;
  ECA_MINOR_REQUEST_ELEMENT_1 var2;
  int compare_type = -1;
};

/**
 * Structure describing a processing request of the form
 *
 * o = F5(F4(F3(F1(i1), F2(i2))))
 *
 * where i1, i2 and o denote the input and output fields,
 * and F1, F2, F3, F4 and F3 are field operators.
 *
 * The structure contains the following elements:
 *
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f1arg     the argument of the 1st field operator
 * f2        the 2nd field operator
 * f2arg     the argument of the 2nd field operator
 * f3        the 3rd field operator
 * f4        the 4th field operator
 * f5        the 5th field operator
 * f5arg     the argument of the 5th field operator
 * epilog    the final operation carried out after processing
 */
struct ECA_MAJOR_REQUEST_ELEMENT_2
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  int refdate = 0;
  ECA_FUNC_1 f1 = nullptr;
  double f1arg = 0.0;
  ECA_FUNC_1 f2 = nullptr;
  double f2arg = 0.0;
  ECA_FUNC_2 f3 = nullptr;
  ECA_FUNC_2 f4 = nullptr;
  ECA_FUNC_3 f5 = nullptr;
  double f5arg = 0.0;
  ECA_EPILOG epilog = ECA_NONE;
};

/**
 * Structure describing a processing request of the form
 *
 * o = H2(H1(i))
 *
 * where i and o denote the input and output fields, and
 * H1, and H2 are field operators.
 *
 * The structure contains the following elements:
 *
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * h1        the 1st field operator
 * h1arg     the argument of the 1st field operator
 * h2        the 2nd field operator
 */
struct ECA_MINOR_REQUEST_ELEMENT_2
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  ECA_FUNC_1 h1 = nullptr;
  double h1arg = 0.0;
  ECA_FUNC_2 h2 = nullptr;
};

struct ECA_REQUEST_2
{
  ECA_MAJOR_REQUEST_ELEMENT_2 var1;
  ECA_MINOR_REQUEST_ELEMENT_2 var2;
  int compare_type = -1;
};

/**
 * Structure describing a processing request of the form
 *
 * o = F3(F1(i1), F2(i2))
 *
 * where i1, i2 and o denote the input and output fields,
 * and F1, F2 and F3 are field operators.
 *
 * The structure contains the following elements:
 *
 * name      the name of the output variable
 * longname  the longname of the output variable
 * units     the units of the output variable
 * f1        the 1st field operator
 * f2        the 2nd field operator
 * f3        the 3rd field operator
 */
struct ECA_REQUEST_3
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  int refdate = 0;
  ECA_FUNC_2 f1 = nullptr;
  ECA_FUNC_2 f2 = nullptr;
  ECA_FUNC_2 f3 = nullptr;
  int compare_type = -1;
};

/**
 * Structure describing a GSL-like processing request. The structure
 * contains the following elements:
 *
 * name       the name of the 1st output variable
 * longname   the longname of the 1st output variable
 * units      the units of the 1st output variable
 * name2      the name of the 2nd output variable
 * longname2  the longname of the 2nd output variable
 * units2     the units of the 2nd output variable
 * name3      the name of the 3rd output variable
 * longname3  the longname of the 3rd output variable
 * units3     the units of the 3rd output variable
 * s1         the 1st field selector
 * s1arg      argument of the 1st field selector
 * s2         the 2nd field selector
 * s2arg      argument of the 2nd field selector
 * consecutiveDays  the number od concecutive days
 */
struct ECA_REQUEST_4
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  const char *name2 = nullptr;
  const char *longname2 = nullptr;
  const char *units2 = nullptr;
  ECA_FUNC_1 s1 = nullptr;
  double s1arg = 0.0;
  ECA_FUNC_1 s2 = nullptr;
  double s2arg = 0.0;
  ECA_FUNC_1 s3 = nullptr;
  double s3arg = 0.0;
  int consecutiveDays = 0;
  int compare_type = -1;
};

/**
 * Structure describing an ETCCDI index processing request. The structure
 * contains the following elements:
 *
 * name       the name of the output variable
 * longname   the longname of the output variable
 * units      the units of the output variable
 * pn         the percentile for the threshold
 * ndates     the window size
 * startboot  the begin year of the bootstrapping interval
 * endboot    the end year of the bootstrapping interval
 * func2      either FieldFunc_Sum or func_average
 */

struct ETCCDI_REQUEST
{
  const char *name = nullptr;
  const char *longname = nullptr;
  const char *units = nullptr;
  int pn = 0;
  int ndates = 0;
  int startboot = 0;
  int endboot = 0;
  int func2 = 0;
  int compare_type = -1;
};

/**
 * Function processing a request of type 1.
 *
 * @param request the processing request
 */
void eca1(const ECA_REQUEST_1 &request);

/**
 * Function processing a request of type 2.
 *
 * @param request the processing request
 */
void eca2(const ECA_REQUEST_2 &request);

/**
 * Function processing a request of type 3.
 *
 * @param request the processing request
 */
void eca3(const ECA_REQUEST_3 &request);

/**
 * Function processing a request of type 4.
 *
 * @param request the processing request
 */
void eca4(const ECA_REQUEST_4 &request);

/**
 * Function processing a etccdi request.
 *
 * @param request the processing request
 */

void etccdi_op(ETCCDI_REQUEST &request);

#endif /* ECA_H_ */
