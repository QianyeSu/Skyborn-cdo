/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author:

*/

#include <fstream>

#include <cdi.h>

#include "process_int.h"
#include "cdo_zaxis.h"
#include "cdi_lockedIO.h"

// NOTE: All operators in this module works only on GRIB edition 1 files!

extern "C" void streamGrbChangeParameterIdentification(int code, int ltype, int lev);

/*
Supported notations:
======================
Selection provided on commandline:
---------------------------------
cdo selmulti,'(33/34;105;10)'
 - or -
cdo selmulti,'(33/34;105;10);(11/6;109;55)'
cdo selmulti,'(33/34;105;10);(11/6;109;40/55)'
cdo selmulti,'(*;105;10);(11/6;109;40/55);(*;105;2)'
cdo selmulti,'{(33/34;105;10);(11/32,8;109;51/52/53/54/55)}'

NOTE: ' .. ' are mandatory !

Selection provided from a text file:
---------------------------------

cdo selmulti,selection_10m_wind.txt

(*A*) Compact general notation, selection file content:

(1; 103; 0)
(33,34; 105; 10)
(11,17; 105; 2)
(71,73,74,75,61,62,65,117,67,122,121,11,131,66,84,111,112; 105; 0)
# If nothing <'sel(' or 'del('>  is specified then
# the operator -selmulti or -delmulti decides if it will be selection of
extraction or delete # Explicite select or delete is also possible:
#(11; 109; *)
#(*; 105; *)
#del(*; 109; *)
#sel(*; 105; *)
#sel(*; 100; *)

# BUT simple array arithmetics should be also possible ("*" ~= mulc;  "+' ~=
addc) sel(33,34;105,1000,3000):math(*2;)        # not implemented yet
sel(11;105,500,1500,3000):math(+273.15;)  # not implemented yet

(*B*) HIP.X notation (KNMI specific), selection file content:

SELECT, PARAMETER=1, LEVTYPE=103, LEVEL=0
SELECT, PARAMETER=33/34, LEVTYPE=105, LEVEL=10
SELECT, PARAMETER=11/17, LEVTYPE=105, LEVEL=2
SELECT, PARAMETER=71/73/74/75/61/62/65/117/67/122/121/11/131/66/84/111/112,
LEVTYPE=105, LEVEL=0 # Explicite delete is also possible: #DELETE,
PARAMETER=128, LEVTYPE=109, LEVEL=*

# BUT simple array arithmetics should be also possible (SCALE ~= mulc;  OFFSET
~= addc)

# The following will convert Pressure from Pa into HPa; Temp from Kelvin to
Celsius: SELECT, PARAMETER=1, LEVTYPE= 103, LEVEL=0, SCALE=0.01 SELECT,
PARAMETER=11, LEVTYPE=105, LEVEL=2, OFFSET=273.15 SELECT, PARAMETER=33/34,
LEVTYPE=105, LEVEL=10

If SCALE and/or OFFSET are defined, then the data values are scaled as
SCALE*(VALUE-OFFSET). The default value for SCALE is 1.0; the default for OFFSET
is 0.0.

**** changemulti ***********

cdo
changemulti,'{(134;1;*|1;105;*);{(6;1;*|6;105;*)};{(246;*;*|76;*;*)};{(247;*;*|58;*;*)};{(248;*;*|71;*;*)}'
fileIN fileOUT


cdo changemulti,'{(134;1;*|1;105;*)}' fileIN fileOUT
# surface pressure has ECMWF code; change it into Hirlam notation ..
grib_set -w indicatorOfParameter=134,indicatorOfTypeOfLevel=1 -s
indicatorOfParameter=1,indicatorOfTypeOfLevel=105 ECMWF_H11_test0.grb
ECMWF_H11_test1.grb

cdo changemulti,'{(6;1;*|6;105;*)}' fileIN fileOUT
# orography has wrong level-type, should be 105
grib_set -w indicatorOfParameter=6,indicatorOfTypeOfLevel=1 -s
indicatorOfParameter=6,indicatorOfTypeOfLevel=105 ECMWF_H11_test1.grb
ECMWF_H11_test2.grb

cdo changemulti,'{(246;*;*|76;*;*)}' fileIN fileOUT
# change code for cloud_water
grib_set -w indicatorOfParameter=246 -s indicatorOfParameter=76
ECMWF_H11_test2.grb ECMWF_H11_test3.grb

cdo changemulti,'{(247;*;*|58;*;*)}' fileIN fileOUT
# change code for cloud_ice
grib_set -w indicatorOfParameter=247 -s indicatorOfParameter=58
ECMWF_H11_test3.grb ECMWF_H11_test4.grb

cdo changemulti,'{(248;*;*|71;*;*)}' fileIN fileOUT
# change code for total_cloud_cover
grib_set -w indicatorOfParameter=248 -s indicatorOfParameter=71
ECMWF_H11_test4.grb ECMWF_H11_test.grb


*/

namespace
{
struct TUPLEREC
{
  std::vector<int> codeLST;
  int ncodes;

  std::vector<int> levelTypeLST;
  int nlevelTypes;

  std::vector<int> levelLST;
  int nlevels;
  int sel_or_del_or_change;  // sel_or_del_or_change:  0:  operator decides,
                             // 1:select , 2:delete, 3:change
  int simpleMath;            // 1:  simple array arithmetics ( *,+), 0: do nothing
  float scale;
  float offset;

  int changedCode;  // used only changemulti mode
  int changedLevelType;
  int changedLevel;
};
}  // namespace

int checkListContainsInt(int value, std::vector<int> const &list, int arraylen);

#define MAX_TUPLES 1000
static int NUMTUPLES = 0;
static TUPLEREC *SelTUPLEREC[MAX_TUPLES];

TUPLEREC *TUPLERECNew();
void push_backSelTuple(TUPLEREC *tp);
TUPLEREC *getSelTuple(int index);

void printSelectionTuples();
int getNumberOfSelectionTuples();
int getNumberOfDeleteSelectionTuples();

int multiSelectionParser(const char *filenameOrString);

class Selmulti : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Selmulti",
    .operators = { { "selmulti", 0, 0, "filename/string with selection specification", SelmultiHelp },
                   { "delmulti", 0, 0, "filename/string with selection specification", SelmultiHelp },
                   { "changemulti", 0, 0, "filename/string with selection specification", SelmultiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Selmulti> registration = RegisterEntry<Selmulti>();

  int SELMULTI{}, DELMULTI{}, CHANGEMULTI{};
  CdoStreamID streamID1{};
  CdoStreamID streamID2{};

  int taxisID1{ CDI_UNDEFID };
  int taxisID2{ CDI_UNDEFID };

  int vlistID1{ CDI_UNDEFID };
  int vlistID2{ CDI_UNDEFID };

  int operatorID{};
  int numVars{};

  bool lcopy = false;

  VarList varList1{};

public:
  void
  init() override
  {
    SELMULTI = module.get_id("selmulti");
    DELMULTI = module.get_id("delmulti");
    CHANGEMULTI = module.get_id("changemulti");

    operatorID = cdo_operator_id();

    operator_input_arg(cdo_operator_enter(operatorID));

    // operator_check_argc(1);
    const char *filenameOrString = cdo_operator_argv(0).c_str();
    if (cdoDebugExt)
    {
      printf("Given operator arguments (nr=%d): \n", cdo_operator_argc());
      for (int i = 0; i < cdo_operator_argc(); ++i) printf("%s", cdo_operator_argv(i).c_str());
      printf("\n");
    }
    if (!multiSelectionParser(filenameOrString))
      cdo_warning("Error processing file with selection description!\n%s", filenameOrString);

    if (operatorID == SELMULTI)
      if (getNumberOfSelectionTuples() == 0)
        cdo_abort("Error! You must provide at lease ONE selection tuple!\n"
                  "Notations: 'SELECT,  .. or sel(/;;) or (/;;)'\nCheck the file: %s",
                  filenameOrString);

    if (operatorID == DELMULTI)
      if (getNumberOfDeleteSelectionTuples() == 0)
        cdo_abort("Error! You must provide at lease ONE selection tuple!\n"
                  "Notations: 'DELETE,  .. or del(/;;) or (/;;)'\nCheck the file: %s",
                  filenameOrString);

    if (operatorID == CHANGEMULTI)
      if (getNumberOfSelectionTuples() == 0)
        cdo_abort("Error! You must provide at lease ONE selection tuple!\n"
                  "Notations: 'CHANGE,  .. or (/;;|;;;)'\nCheck the file: %s",
                  filenameOrString);

    streamID1 = cdo_open_read(0);

    vlistID1 = cdo_stream_inq_vlist(streamID1);

    varList1 = VarList(vlistID1);

    vlistClearFlag(vlistID1);
    numVars = varList1.numVars();

    Debug(cdoDebugExt, " Total number of variables: %d", numVars);

    for (int varID = 0; varID < numVars; ++varID)
    {
      auto const &var = varList1.vars[varID];
      auto ltype = zaxis_to_ltype(var.zaxisID);

      for (int levelID = 0; levelID < var.nlevels; ++levelID)
      {
        auto level = zaxisInqLevel(var.zaxisID, levelID);

        if (operatorID == DELMULTI) vlistDefFlag(vlistID1, varID, levelID, true);  // set initially, override bellow if in selection
        if (operatorID == CHANGEMULTI)
        {
          vlistDefFlag(vlistID1, varID, levelID, true);  // change operation copies all fields
          continue;
        }

        for (int ii = 0; ii < NUMTUPLES; ++ii)
        {
          TUPLEREC *tuplerec = getSelTuple(ii);
          // if ( cdoDebugExt ) cdo_print(" Processing: (code %d,
          // ltype %d, level %d);  nvars=%d, varID=%d", code, ltype,
          // (int)level, nvars, varID);
          // Note: When the list is Empty then function
          // checkListContainsInt() also returns true !
          auto selcode = checkListContainsInt(var.code, tuplerec->codeLST, tuplerec->ncodes);
          auto selltype = checkListContainsInt(ltype, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
          auto sellevel = checkListContainsInt((int) level, tuplerec->levelLST, tuplerec->nlevels);
          if (selcode && selltype && sellevel)
          {
            if (operatorID == SELMULTI)
            {
              switch (tuplerec->sel_or_del_or_change)
              {
                case 0:  // operator decides ...
                  vlistDefFlag(vlistID1, varID, levelID, true);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code, ltype,
                                (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; "
                                "OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
                case 1:
                  vlistDefFlag(vlistID1, varID, levelID, true);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code, ltype,
                                (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; "
                                "OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
                case 2:
                  vlistDefFlag(vlistID1, varID, levelID, false);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code,
                                ltype, (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   "
                                "[varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
              }
            }
            else if (operatorID == DELMULTI)
            {
              switch (tuplerec->sel_or_del_or_change)
              {
                case 0:  // operator decides ...
                  vlistDefFlag(vlistID1, varID, levelID, false);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code,
                                ltype, (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   "
                                "[varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
                case 1:
                  vlistDefFlag(vlistID1, varID, levelID, true);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code, ltype,
                                (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting : (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; SCALE=%f; "
                                "OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
                case 2:
                  vlistDefFlag(vlistID1, varID, levelID, false);
                  if (cdoDebugExt)
                  {
                    if (!tuplerec->simpleMath)
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", var.code,
                                ltype, (int) (level), varID, levelID);
                    else
                      cdo_print(" Selecting for removal: (code %3i, ltype %3i, level %3i)   "
                                "[varID(%d),levelID(%d)]; SCALE=%f; OFFSET=%f",
                                var.code, ltype, (int) (level), varID, levelID, tuplerec->scale, tuplerec->offset);
                  }
                  break;
              }
            }
            break;
          }
        }  // end for ( .. NUMTUPLES
      }    // end for ( levelID
    }      // end for ( varID

    Debug(cdoDebugExt, " Writing the selected fields ...");

    vlistID2 = vlistCreate();
    cdo_vlist_copy_flag(vlistID2, vlistID1);
    vlistDefNtsteps(vlistID2, vlistNtsteps(vlistID1));

    {
      numVars = vlistNvars(vlistID2);
      int varID = 0;
      for (; varID < numVars; ++varID)
        if (vlistInqVarTimetype(vlistID2, varID) != TIME_CONSTANT) break;
      if (varID == numVars) vlistDefNtsteps(vlistID2, 0);
    }

    taxisID1 = vlistInqTaxis(vlistID1);
    taxisID2 = taxisDuplicate(taxisID1);
    vlistDefTaxis(vlistID2, taxisID2);

    streamID2 = cdo_open_write(1);
    cdo_def_vlist(streamID2, vlistID2);
  }

  void
  run() override
  {
    auto gridsizeMax = varList1.gridsizeMax();
    if (vlistNumber(vlistID1) != CDI_REAL) gridsizeMax *= 2;
    Varray<double> array(gridsizeMax);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;

      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);

      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);
        auto const &var = varList1.vars[varID];

        if (vlistInqFlag(vlistID1, varID, levelID) == true)
        {
          int simpleMath = 0;  // 1:  simple array arithmetics ( *,+), 0: do nothing
          double scale = 1.0;
          double offset = 0.0;
          auto code = var.code;
          auto level = zaxisInqLevel(var.zaxisID, levelID);
          auto ltype = zaxis_to_ltype(var.zaxisID);
          for (int ii = 0; ii < NUMTUPLES; ++ii)
          {
            TUPLEREC *tuplerec = getSelTuple(ii);
            // Note: When the list is Empty then function
            // checkListContainsInt() also returns true !
            auto selcode = checkListContainsInt(code, tuplerec->codeLST, tuplerec->ncodes);
            auto selltype = checkListContainsInt(ltype, tuplerec->levelTypeLST, tuplerec->nlevelTypes);
            auto sellevel = checkListContainsInt((int) level, tuplerec->levelLST, tuplerec->nlevels);
            lcopy = true;
            if (selcode && selltype && sellevel)
            {
              if (operatorID == CHANGEMULTI)
              {
                if (cdoDebugExt)
                  cdo_print(" Processing: (code %d, ltype %d, level %d);  nvars=%d, varID=%d => (selcode %d, selltype "
                            "%d, sellevel %d) => change (%d,%d,%d)",
                            code, ltype, (int) level, numVars, varID, selcode, selltype, sellevel, tuplerec->changedCode,
                            tuplerec->changedLevelType, tuplerec->changedLevel);
                if ((tuplerec->changedCode == -1) && (tuplerec->changedLevelType == -1) && (tuplerec->changedLevel == -1))
                  cdo_print(" WARNING: Cannot CHANGE identification!");
                else
                  streamGrbChangeParameterIdentification(tuplerec->changedCode, tuplerec->changedLevelType, tuplerec->changedLevel);
                // Calling PROXY function
                // streamGrbChangeParameterIdentification() which
                // results in later calling func.
                // gribapiChangeParameterIdentification(); see
                // stream_gribapi.c The change happens during the last
                // step of writing a grib-record into a file.
              }
              else
              {
                if (cdoDebugExt)
                  cdo_print(" Processing: (code %d, ltype %d, level %d);  nvars=%d, varID=%d => (selcode %d, "
                            "selltype %d, sellevel %d)",
                            code, ltype, (int) level, numVars, varID, selcode, selltype, sellevel);
              }
              simpleMath = tuplerec->simpleMath;  // 1:  simple array arithmetics ( *,+),
                                                  // 0: do nothing
              if (simpleMath)
              {
                scale = tuplerec->scale;
                offset = tuplerec->offset;
                lcopy = false;
              }
              break;  // get out of this for loop
            }
          }  // end of for (ii=0; ii<NUMTUPLES ..

          auto varID2 = vlistFindVar(vlistID2, varID);
          auto levelID2 = vlistFindLevel(vlistID2, varID, levelID);

          // tijdelijk PATCH M.K.
          if ((varID2 == -1) || (levelID2 == -1))
          {
            cdo_print(" Warning: Missing varID or levelID with (code %3i, "
                      "ltype %3i, level %3i)   [varID(%d),levelID(%d)] .. #2[varID(%d),levelID(%d)]",
                      code, ltype, (int) (level), varID, levelID, varID2, levelID2);
            continue;
          }
          cdo_def_field(streamID2, varID2, levelID2);

          if (lcopy)
          {
            if (cdoDebugExt)
              cdo_print(" Copying record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", fieldID, code,
                        ltype, (int) (level), varID, levelID);
            cdo_copy_field(streamID1, streamID2);
          }
          else
          {
            size_t numMissVals;
            cdo_read_field(streamID1, array.data(), &numMissVals);

            if (!simpleMath)
            {
              if (cdoDebugExt)
                cdo_print(" Writing record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]", fieldID, code,
                          ltype, (int) (level), varID, levelID);
            }
            else  // 1:  simple array arithmetics ( *,+)
            {
              if (cdoDebugExt)
                cdo_print(" Writing record [%4d] with (code %3i, ltype %3i, level %3i)   [varID(%d),levelID(%d)]; "
                          "SCALE=%f; OFFSET=%f",
                          fieldID, code, ltype, (int) (level), varID, levelID, scale, offset);
              for (size_t li = 0; li < gridsizeMax; ++li)
                if (!fp_is_equal(array[li], var.missval))
                {
                  array[li] = scale * (array[li] - offset);
                  // If SCALE and/or OFFSET are defined, then the data
                  // values are scaled as SCALE*(VALUE-OFFSET).
                }
            }

            cdo_write_field(streamID2, array.data(), numMissVals);
          }  // end of else ( lcopy )
        }    // end if ( vlistInqFlag(vlistID1, varID, levelID) == true )
      }      // end for ( fieldID ..

      tsID++;
    }  // end while
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
    cdo_stream_close(streamID2);

    vlistDestroy(vlistID2);

    cdoDebugExt = 0;
  }
};

TUPLEREC *
TUPLERECNew()
{
  TUPLEREC *tpl = (TUPLEREC *) std::malloc(sizeof(TUPLEREC));

  tpl->ncodes = 0;

  tpl->nlevelTypes = 0;

  tpl->nlevels = 0;
  tpl->sel_or_del_or_change = 0;  // sel_or_del_or_change:  0:  operator
                                  // decides, 1:select , 2:delete, 3:change

  tpl->simpleMath = 0;  // 1:  simple array arithmetics ( *,+), 0: do nothing
  tpl->scale = 1.0;
  tpl->offset = 0.0;

  tpl->changedCode = -1;  // used only changemulti mode
  tpl->changedLevelType = -1;
  tpl->changedLevel = -1;

  return tpl;
}

void
push_backSelTuple(TUPLEREC *tp)
{
  if (NUMTUPLES < MAX_TUPLES)
  {
    SelTUPLEREC[NUMTUPLES] = tp;
    NUMTUPLES++;
  }
}

TUPLEREC *
getSelTuple(int index)
{
  if (index < NUMTUPLES) { return SelTUPLEREC[index]; }
  return nullptr;
}

int
checkListContainsInt(int value, std::vector<int> const &list, int arraylen)
// Note: When the list is Empty it also returns true !
{
  if (arraylen == 0)
  {
    // if ( cdoDebugExt ) cdo_print(" value=%d: found. (list empty)");
    return 1;
  }
  for (int i = 0; i < arraylen; ++i)
  {
    if (list[i] == -1)
    {  // this is for '*' selection; can be any code, any level or any
       // level-type
      // if ( cdoDebugExt ) cdo_print(" value=%d: found.");
      return 1;
    }

    if (list[i] == value)
    {
      // if ( cdoDebugExt ) cdo_print(" value=%d: found.");
      return 1;
    }
  }
  // if ( cdoDebugExt ) cdo_print(" value=%d: NOT found.");
  return 0;
}

#define MAX_LINE_LEN 65536

static char *
removeSpaces(char *pline)
{
  if (pline == nullptr) return nullptr;
  while (std::isspace((int) *pline)) pline++;
  return pline;
}

static char *
skipSeparator(char *pline)
{
  if (pline == nullptr) return nullptr;
  while (std::isspace((int) *pline)) pline++;
  if (*pline == '=' || *pline == ':' || *pline == '/' || *pline == ',') pline++;
  while (std::isspace((int) *pline)) pline++;
  return pline;
}

static char *
goToNextSeparator(char *pline)
{
  if (pline == nullptr) return nullptr;
  int separatorFound = 0;
  while (std::isspace((int) *pline) || !separatorFound)
  {
    if (*pline == '\0') return nullptr;
    pline++;
    if (*pline == '|') return nullptr;
    if (*pline == '=' || *pline == ':' || *pline == '/' || *pline == ',')
      separatorFound = 1;
    else if (*pline == ';')
    {
      pline++;
      pline = removeSpaces(pline);
      return (pline);
    }
  }
  if (separatorFound) pline++;
  Debug(cdoDebugExt >= 100, "goToNextSeparator():  pline= ('%s') ", pline);
  // while ( std::isspace((int) *pline) ) pline++;
  pline = removeSpaces(pline);
  return pline;
}

static char *
strContains(char *str, const char *substr)
{
  if (str == nullptr) return nullptr;
  if (substr == nullptr) return nullptr;

  str = removeSpaces(str);

  size_t lensub = std::strlen(substr);
  size_t lenstr = std::strlen(str);

  if (lensub > lenstr)
  {
    if (cdoDebugExt >= 100)
      cdo_print("strContains():  substr('%s') NOT found in str('%s');  lensub(%zu)>lenstr(%zu) ", substr, str, lensub, lenstr);
    return nullptr;
  }
  char *rv = std::strstr(str, substr);
  if (rv)
  {
    Debug(cdoDebugExt >= 100, "strContains():  substr('%s') FOUND in str('%s')", substr, str);
    return (rv + lensub);  // points after subStr ..
  }
  else
  {
    Debug(cdoDebugExt >= 100, "strContains():  substr('%s') NOT found in str('%s')", substr, str);
    return rv;
  }
}

static char *
findParamEnd(char *str)
{
  char *ptr = str;
  char *ptrEnding = nullptr;

  if (str == nullptr) return nullptr;
  // supported endings are: ", " or ";"
  if (ptrEnding == nullptr) ptrEnding = strContains(ptr, ", ");  // HIP notation
  if (ptrEnding == nullptr) ptrEnding = strContains(ptr, ";");   // compact notation
  if (ptrEnding != nullptr)
  {
    ptrEnding = removeSpaces(ptrEnding);
    Debug(cdoDebugExt >= 100, " ptrEnding='%s'", ptrEnding);
    return ptrEnding;
  }
  Debug(cdoDebugExt >= 100, " ptrEnding=end-of-string");
  size_t lenstr = std::strlen(str);
  ptrEnding = str + lenstr;
  ptrEnding = removeSpaces(ptrEnding);
  return ptrEnding;
}

static char *
findTupleEnd(char *str)
{
  char *ptr = str;
  char *ptrEnding;

  if (str == nullptr) return nullptr;
  // supported endings are: ")" or ");"
  ptrEnding = strContains(ptr, ")");
  if (ptrEnding == nullptr) ptrEnding = strContains(ptr, ");");
  if (ptrEnding != nullptr)
  {
    ptrEnding = removeSpaces(ptrEnding);
    Debug(cdoDebugExt >= 100, " findTupleEnd='%s'", ptrEnding);
    return ptrEnding;
  }
  Debug(cdoDebugExt >= 100, " findTupleEnd=end-of-string");
  size_t lenstr = std::strlen(str);
  ptrEnding = str + lenstr;
  ptrEnding = removeSpaces(ptrEnding);
  return ptrEnding;
}

static char *
readlineForParsing(std::ifstream &gfile, char *strToParsePtr, char *line)
{
  if (gfile.is_open())  // file is open => we parse text from a file
  {
    std::string strLine;
    if (std::getline(gfile, strLine))
    {
      strcpy(line, strLine.c_str());
      return line;
    }
    else { return nullptr; }
  }
  else if (strToParsePtr != nullptr)  // we parse a given string
  {
    if (cdoDebugExt >= 30) cdo_print("%s(): Parsing selection string:  %s", __func__, strToParsePtr);
    char *tpEnd = nullptr;
    if (std::strlen(strToParsePtr) > 0) tpEnd = findTupleEnd(strToParsePtr);
    if (tpEnd == nullptr)
    {
      Debug(cdoDebugExt >= 100, "%s(): End of selection string reached.", __func__);
      return nullptr;
    }
    else
    {
      tpEnd[0] = 0;
      if (std::strlen(strToParsePtr) <= MAX_LINE_LEN) std::strcpy(line, strToParsePtr);
      Debug(cdoDebugExt >= 100, "%s(): Current selection line=%s", __func__, line);
      strToParsePtr = tpEnd + 1;
      return strToParsePtr;
    }
  }
  else
  {
    cdo_abort(" Cannot parse selection string:  %s", strToParsePtr);
    return nullptr;
  }
}

int
multiSelectionParser(const char *filenameOrString)
{
  std::ifstream gfile;
  char *strpos;
  char *parEnd;
  int val;
  float floatval;
  TUPLEREC *tuplerec;
  int selectionRec;
  char strToParse[MAX_LINE_LEN];
  char *strToParsePtr = nullptr;
  strToParse[0] = 0;
  int convertLevelsHPa2Pa = 0;  // this is HIP - KNMI selection specific
  char first3chars[4];
  std::strncpy(first3chars, filenameOrString, 3);

  if ((filenameOrString[0] == '{') || (filenameOrString[0] == '(') || strContains(first3chars, "del")
      || strContains(first3chars, "sel"))
  {
    // cdo selmulti,'(33/34;105;10)'
    // - or -
    // cdo selmulti,'{(33/34;105;10);(11/32,8;109;51/52/53/54/55)}'
    // - or -
    // cdo changemulti,'{(134;1;0|1;105;0);(6;1;0|6;105;0)}'
    std::strncpy(strToParse, filenameOrString, MAX_LINE_LEN - 1);
    strToParsePtr = &strToParse[0];
    if (strToParsePtr[0] == '{') strToParsePtr++;
    int strLn = std::strlen(strToParsePtr);
    if (strToParsePtr[strLn - 1] == '}') strToParsePtr[strLn - 1] = 0;
    Debug(cdoDebugExt, " Parsing selection string:  %s", strToParsePtr);
  }
  else
  {
    gfile.open(filenameOrString);
    Debug(cdoDebugExt, " Parsing file:  %s", filenameOrString);
    if (!gfile.is_open())
    {
      cdo_abort(" Missing file:  %s", filenameOrString);
      return 0;
    }
  }

  char line[MAX_LINE_LEN];
  while ((strToParsePtr = readlineForParsing(gfile, strToParsePtr, line)))
  {
    if (line[0] == '#') continue;
    if (line[0] == '\0') continue;
    char *pline = line;
    if (cdoDebugExt >= 30) cdo_print(": Line: %s", pline);
    while (std::isspace((int) *pline)) pline++;
    if (pline[0] == '\0') continue;
    strpos = strContains(pline, "SELECT, ");
    selectionRec = 0;  // default is 0;  sel_or_del_or_change:  0:  operator
                       // decides, 1:select , 2:delete, 3:change
    if (strpos != nullptr)
      selectionRec = 1;
    else
    {
      strpos = strContains(pline, "DELETE, ");
      if (strpos != nullptr)
        selectionRec = 2;
      else
      {
        strpos = strContains(pline, "REMOVE, ");
        if (strpos != nullptr) selectionRec = 2;
      }
    }
    if (strpos != nullptr)  // we have SELECT ..
    {
      Debug(cdoDebugExt, " Parsing notation SELECT: %s", strpos);
      pline = strpos;
      tuplerec = TUPLERECNew();
      tuplerec->sel_or_del_or_change = selectionRec;
      push_backSelTuple(tuplerec);
      // SELECT, PARAMETER=11/17, LEVTYPE=105, LEVEL=2
      while ((pline != nullptr) && (std::strlen(pline) != 0))
      {
        Debug(cdoDebugExt >= 100, "pline='%s'", pline);
        strpos = strContains(pline, "PARAMETER=");
        if ((strpos) != nullptr)
        {
          Debug(cdoDebugExt >= 100, ": PARAMETER=%s", strpos);
          pline = strpos;
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            if (pline[0] == '*')
              val = -1;
            else
              val = atoi(pline);
            Debug(cdoDebugExt >= 100, "code=%d", val);
            tuplerec->codeLST.push_back(val);
            tuplerec->ncodes = tuplerec->codeLST.size();
            strpos = goToNextSeparator(pline);
            pline = strpos;
            if (!strpos) break;
            pline = skipSeparator(strpos);
          }
        }
        Debug(cdoDebugExt >= 100, "pline='%s'", pline);
        strpos = strContains(pline, "LEVTYPE=");
        if ((strpos) != nullptr)
        {
          Debug(cdoDebugExt >= 100, ": LEVTYPE=%s", strpos);
          pline = strpos;
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            if (pline[0] == '*')
              val = -1;
            else
              val = atoi(pline);
            if (val == 100)
            {
              convertLevelsHPa2Pa = 1;
              if (cdoDebugExt >= 1)
                cdo_print("Detected levelType=100 ! HIP selection specifies HPa but CDO has levels in "
                          "Pascals.\nSelection HPa's will be converted into Pa's.");
            }
            else
              convertLevelsHPa2Pa = 0;
            Debug(cdoDebugExt >= 100, "levelType=%d", val);
            tuplerec->levelTypeLST.push_back(val);
            tuplerec->nlevelTypes = tuplerec->levelTypeLST.size();
            strpos = goToNextSeparator(pline);
            pline = strpos;
            if (!strpos) break;
            pline = skipSeparator(strpos);
          }
        }
        Debug(cdoDebugExt >= 100, "pline='%s'", pline);
        strpos = strContains(pline, "LEVEL=");
        if ((strpos) != nullptr)
        {
          Debug(cdoDebugExt >= 100, ": LEVEL=%s", strpos);
          pline = strpos;
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            if ((convertLevelsHPa2Pa) && (val != -1)) val *= 100;
            Debug(cdoDebugExt >= 100, "level=%d", val);
            tuplerec->levelLST.push_back(val);
            tuplerec->nlevels = tuplerec->levelLST.size();
            strpos = goToNextSeparator(pline);
            pline = strpos;
            if (!strpos) break;
            pline = skipSeparator(strpos);
          }
        }
        Debug(cdoDebugExt >= 100, "pline='%s' (check SCALE=...)", pline);
        strpos = strContains(pline, "SCALE=");
        if ((strpos) != nullptr)
        {
          Debug(cdoDebugExt >= 100, ": SCALE= %s", strpos);
          pline = strpos;
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            floatval = atof(pline);
            Debug(cdoDebugExt >= 100, "scale=%f", floatval);
            tuplerec->simpleMath = 1;    // 1:  simple array arithmetics
                                         // ( *,+), 0: do nothing
            tuplerec->scale = floatval;  // tuplerec->offset = 0.0;
            strpos = goToNextSeparator(pline);
            pline = strpos;
            if (!strpos) break;
            pline = skipSeparator(strpos);
          }
        }
        Debug(cdoDebugExt >= 100, "pline='%s' (check OFFSET=...)", pline);
        strpos = strContains(pline, "OFFSET=");
        if ((strpos) != nullptr)
        {
          Debug(cdoDebugExt >= 100, ": OFFSET= %s", strpos);
          pline = strpos;
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            floatval = atof(pline);
            Debug(cdoDebugExt >= 100, "offset=%f", floatval);
            tuplerec->simpleMath = 1;     // 1:  simple array arithmetics
                                          // ( *,+), 0: do nothing
            tuplerec->offset = floatval;  // tuplerec->scale = 1.0;
            strpos = goToNextSeparator(pline);
            pline = strpos;
            if (!strpos) break;
            pline = skipSeparator(strpos);
          }
        }
      }  // end while pline
      continue;
    }  // end if SELECT,

    // Here comes the short notation
    selectionRec = 0;  // default is 0;  sel_or_del_or_change:  0:  operator
                       // decides, 1:select , 2:delete, 3:change
    strpos = strContains(pline, "sel(");
    if (strpos != nullptr)
      selectionRec = 1;
    else
    {
      strpos = strContains(pline, "del(");
      if (strpos != nullptr)
        selectionRec = 2;
      else
      {
        strpos = strContains(pline, "|");
        if (strpos != nullptr)
          selectionRec = 3;
        else
        {
          strpos = strContains(pline, "(");
          if (strpos != nullptr) selectionRec = 0;
        }
      }
    }
    if (strpos != nullptr)
    {
      if (cdoDebugExt)
        cdo_print(" Parsing notation (code(s),..; levelType(s),..; level(s),..) : %s; [selectionRec =%d]", strpos, selectionRec);
      if (selectionRec == 3) strpos = strContains(pline, "(");
      pline = strpos;
      tuplerec = TUPLERECNew();
      tuplerec->sel_or_del_or_change = selectionRec;
      push_backSelTuple(tuplerec);
      // (33/34; 105; 10)
      while ((pline != nullptr) && (std::strlen(pline) != 0) && (pline[0] != ')'))
      {
        Debug(cdoDebugExt >= 100, "[1]: pline='%s'", pline);
        // 1st is code
        {
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            if (pline[0] == '*')
              val = -1;
            else
              val = atoi(pline);
            Debug(cdoDebugExt >= 100, "code=%d", val);
            tuplerec->codeLST.push_back(val);
            tuplerec->ncodes = tuplerec->codeLST.size();
            strpos = goToNextSeparator(pline);
            if (!strpos)
            {
              pline = parEnd;
              break;
            }
            else
              pline = strpos;
            pline = skipSeparator(strpos);
          }
        }
        // 2nd is level type
        Debug(cdoDebugExt >= 100, "[2]: pline='%s'", pline);
        {
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            if (val == 100)
            {
              convertLevelsHPa2Pa = 1;
              if (cdoDebugExt >= 1)
                cdo_print("Detected levelType=100 ! Selection specifies HPa but CDO has levels in "
                          "Pascals.\nSelection HPa's will be converted into Pa's.");
            }
            else
              convertLevelsHPa2Pa = 0;

            Debug(cdoDebugExt >= 100, "levelType=%d", val);
            tuplerec->levelTypeLST.push_back(val);
            tuplerec->nlevelTypes = tuplerec->levelTypeLST.size();
            strpos = goToNextSeparator(pline);
            if (!strpos)
            {
              pline = parEnd;
              break;
            }
            else
              pline = strpos;
            pline = skipSeparator(strpos);
          }
        }
        // 3rd is level
        Debug(cdoDebugExt >= 100, "[3]: pline='%s'", pline);
        {
          pline = removeSpaces(pline);
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            if ((convertLevelsHPa2Pa) && (val != -1)) val *= 100;
            Debug(cdoDebugExt >= 100, "level=%d", val);
            tuplerec->levelLST.push_back(val);
            tuplerec->nlevels = tuplerec->levelLST.size();
            Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            strpos = goToNextSeparator(pline);
            Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            if (!strpos)
            {
              strpos = strContains(pline, "|");  // compact notation for  changemulti
              if (strpos)
                pline = strpos - 1;  // strContains returns character after...
              else
                pline = parEnd;
              Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
              break;
            }
            else
              pline = strpos;
            pline = skipSeparator(strpos);
          }
        }

        // OPTIONAL:
        // cdo
        // changemulti,'{(134;1;*|1;105;*);{(6;1;*|6;105;*)};{(246;*;*|76;*;*)};{(247;*;*|58;*;*)};{(248;*;*|71;*;*)}'
        // fileIN fileOUT
        Debug(cdoDebugExt >= 100, "[OPT]: pline='%s'", pline);
        {
          pline = removeSpaces(pline);
          // pline points to: "=1;105;*);....."
          if (pline[0] == '|')
          {  // changemulti specification
            tuplerec->sel_or_del_or_change = 3;
            pline = &pline[1];
            // Get changedCode:
            parEnd = findParamEnd(pline);
            if ((!parEnd) || (pline[0] == 0)) cdo_abort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            tuplerec->changedCode = val;
            Debug(cdoDebugExt >= 100, "changedCode=%d", val);
            strpos = goToNextSeparator(pline);
            if (!strpos) cdo_abort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            pline = skipSeparator(strpos);
            // Get changedLevelType:
            parEnd = findParamEnd(pline);
            if ((!parEnd) || (pline[0] == 0)) cdo_abort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            tuplerec->changedLevelType = val;
            Debug(cdoDebugExt >= 100, "changedLevelType=%d", val);
            strpos = goToNextSeparator(pline);
            if (!strpos) cdo_abort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            pline = skipSeparator(strpos);
            // Get changedLevel:
            parEnd = findParamEnd(pline);
            if ((!parEnd) || (pline[0] == 0)) cdo_abort("Channot parse: strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            tuplerec->changedLevel = val;
            Debug(cdoDebugExt >= 100, "changedLevel=%d", val);
            pline = parEnd;
          }  // changemulti specification
          parEnd = findParamEnd(pline);
          if ((!parEnd) && (!pline)) Debug(cdoDebugExt >= 100, "strpos=%s; parEnd=%s ... pline=%s", strpos, parEnd, pline);
          while ((pline != parEnd) && (std::strlen(pline) > 0))
          {
            pline = removeSpaces(pline);
            val = (pline[0] == '*') ? -1 : atoi(pline);
            if ((convertLevelsHPa2Pa) && (val != -1)) val *= 100;
            Debug(cdoDebugExt >= 100, "level=%d", val);
            tuplerec->levelLST.push_back(val);
            tuplerec->nlevels = tuplerec->levelLST.size();
            strpos = goToNextSeparator(pline);
            if (!strpos)
            {
              pline = parEnd;
              break;
            }
            else
              pline = strpos;
            pline = skipSeparator(strpos);
          }
        }

      }  // end while pline
    }    // end if "("
  }      // end while ( readlineForParsing() )

  if (gfile.is_open()) gfile.close();

  printSelectionTuples();
  return 1;
}

void
printSelectionTuples()
{
  Debug(cdoDebugExt, " Printing selection tuples:");

  int ii, ri;
  for (ii = 0; ii < NUMTUPLES; ++ii)
  {
    TUPLEREC *tuplerec = getSelTuple(ii);
    char strval[1000];
    std::strcpy(strval, "(");
    char bff[200];
    bff[0] = '\0';
    if (cdoDebugExt)
      cdo_print(" Selection tuple [%d]: ncodes=%d, nlevelTypes=%d, nlevels=%d", ii, tuplerec->ncodes, tuplerec->nlevelTypes,
                tuplerec->nlevels);
    for (ri = 0; ri < tuplerec->ncodes; ri++)
    {
      std::snprintf(bff, sizeof(bff), "%d", tuplerec->codeLST[ri]);
      std::strcat(strval, bff);
      if ((ri + 1) < tuplerec->ncodes)
        std::strcat(strval, "/");
      else
        std::strcat(strval, ";");
    }
    for (ri = 0; ri < tuplerec->nlevelTypes; ri++)
    {
      std::snprintf(bff, sizeof(bff), "%d", tuplerec->levelTypeLST[ri]);
      std::strcat(strval, bff);
      if ((ri + 1) < tuplerec->nlevelTypes)
        std::strcat(strval, "/");
      else
        std::strcat(strval, ";");
    }
    for (ri = 0; ri < tuplerec->nlevels; ri++)
    {
      std::snprintf(bff, sizeof(bff), "%d", tuplerec->levelLST[ri]);
      std::strcat(strval, bff);
      if ((ri + 1) < tuplerec->nlevels)
        std::strcat(strval, "/");
      else
        std::strcat(strval, ")");
    }

    if (tuplerec->simpleMath)
    {
      std::snprintf(bff, sizeof(bff), " {scale = %f; offset = %f}", tuplerec->scale, tuplerec->offset);
      std::strcat(strval, bff);
    }

    // sel_or_del_or_change:  0:  operator decides, 1:select , 2:delete,
    // 3:change
    if (tuplerec->sel_or_del_or_change == 1) { Debug(cdoDebugExt, " Selection tuple [%d] = %s (select)", ii, strval); }
    else if (tuplerec->sel_or_del_or_change == 2) { Debug(cdoDebugExt, " Selection tuple [%d] = %s (delete)", ii, strval); }
    if (tuplerec->sel_or_del_or_change == 3)
    {
      if (cdoDebugExt)
        cdo_print(" Selection tuple [%d] = %s (change) => (%d; %d;%d;)", ii, strval, tuplerec->changedCode,
                  tuplerec->changedLevelType, tuplerec->changedLevel);
    }
    else
      cdo_print(" Selection tuple [%d] = %s (select/delete)", ii, strval);
  }
}

int
getNumberOfSelectionTuples()
{
  int nn = 0;
  for (int ii = 0; ii < NUMTUPLES; ++ii)
  {
    TUPLEREC *tuplerec = getSelTuple(ii);
    if (tuplerec->sel_or_del_or_change != 2) nn++;
  }
  return nn;
}

int
getNumberOfDeleteSelectionTuples()
{
  int nn = 0;
  cdo_print("NUMTUPLES = %d", NUMTUPLES);
  printSelectionTuples();
  for (int ii = 0; ii < NUMTUPLES; ++ii)
  {
    TUPLEREC *tuplerec = getSelTuple(ii);
    cdo_print("getNumberOfDeleteSelectionTuples() [%d]=%d\n", ii, tuplerec->sel_or_del_or_change);
    if (tuplerec->sel_or_del_or_change != 1) nn++;
  }
  return nn;
}
