/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Copyright (C) 2006 Brockmann Consult

  Author: Ralf Quast

*/

// clang-format off
/*
      MODULE      OPERATOR     INDEX    DESCRIPTION

      EcaCfd      eca_cfd      CFD      maximum number of consecutive frost days
      EcaCsu      eca_csu      CSU      maximum number of consecutive summer days
      EcaCwdi     eca_cwdi     CWDI     cold wave duration index
      EcaCwfi     eca_cwfi     CWFI     number of cold-spell days
      EcaEtr      eca_etr      ETR      intra-period extreme temperature range
      EcaFd       eca_fd       FD       number of frost days
      EcaGsl      eca_gsl      GSL      growing season length
      EcaHd       eca_hd       HD       heating degree days
      EcaHwdi     eca_hwdi     HWDI     heat wave duration index
      EcaHwfi     eca_hwfi     HWFI     number of warm-spell days
      EcaId       eca_id       ID       number of ice days
      EcaSu       eca_su       SU       number of summer days
      EcaTg10p    eca_tg10p    TG10p    percent of time TX < 10th percentile of daily mean temperature
      EcaTg90p    eca_tg90p    TG90p    percent of time TX > 90th percentile of daily mean temperature
      EcaTn10p    eca_tn10p    TN10p    percent of time TX < 10th percentile of daily minimum temperature
      EcaTn90p    eca_tn90p    TN90p    percent of time TX > 90th percentile of daily minimum temperature
      EcaTr       eca_tr       TR       number of tropical nights
      EcaTx10p    eca_tx10p    TX10p    percent of time TX < 10th percentile of daily maximum temperature
      EcaTx90p    eca_tx90p    TX90p    percent of time TX > 90th percentile of daily maximum temperature

      EcaCdd      eca_cdd      CDD      maximum number of consecutive dry days
      EcaCwd      eca_cwd      CWD      maximum number of consecutive wet days
      EcaR10mm    eca_r10mm    R10mm    number of days with precipitation >= 10 mm
      EcaR20mm    eca_r20mm    R20mm    number of days with precipitation >= 20 mm
      EcaR75p     eca_r75p     R75p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR75ptot  eca_r75ptot  R75pTOT  Percentage of annual total precipitation due to events with RR > 75th percentile of daily precipitation amount
      EcaR90p     eca_r90p     R90p     Percent of time RR > 90th percentile of daily precipitation amount
      EcaR90ptot  eca_r90ptot  R90pTOT  Percentage of annual total precipitation due to events with RR > 90th percentile of daily precipitation amount
      EcaR95p     eca_r95p     R95p     Percent of time RR > 95th percentile of daily precipitation amount
      EcaR95ptot  eca_r95ptot  R95pTOT  Percentage of annual total precipitation due to events with RR > 95th percentile of daily precipitation amount
      EcaR99p     eca_r99p     R99p     Percent of time RR > 75th percentile of daily precipitation amount
      EcaR99ptot  eca_r99ptot  R99pTOT  Percentage of annual total precipitation due to events with RR > 99th percentile of daily precipitation amount
      EcaRr1      eca_rr1      RR1      number of wet days
      EcaSdii     eca_sdii     SDII     simple daily intensity index

      Fdns        fdns                  frost days without surface snow

      Strwin      strwin                number of strong-wind days
      Strbre      strbre                number of strong-breeze days
      Strgal      strgal                number of strong-gale days
      Hurr        hurr                  number of hurricane days
*/
// clang-format on

#include "process_int.h"
#include "cdo_options.h"
#include "param_conversion.h"
#include "ecacore.h"
#include "ecautil.h"
#include "util_date.h"
#include "pmlist.h"
#include "field_functions.h"

#define TO_DEG_CELSIUS(x) ((x) -273.15)
#define TO_KELVIN(x) ((x) + 273.15)

constexpr int ECA_refdate = 19550101;
constexpr int ETC_refdate = 18500101;

// clang-format off

static const char CFD_NAME[]         = "consecutive_frost_days_index_per_time_period";
static const char CFD_LONGNAME[]     = "Consecutive frost days index is the greatest number of consecutive frost days in a given time period. Frost days is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
//static const char CFD_UNITS[]        = "No.";
static const char CFD_NAME2[]        = "number_of_cfd_periods_with_more_than_%ddays_per_time_period";
static const char CFD_LONGNAME2[]    = "Number of cfd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CFD_UNITS2[]       = "No.";

static const char CSU_NAME[]         = "consecutive_summer_days_index_per_time_period";
static const char CSU_LONGNAME[]     = "Consecutive summer days index is the greatest number of consecutive summer days in a given time period. Summer days is the number of days where maximum of temperature is above 25 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
//static const char CSU_UNITS[]        = "No.";
static const char CSU_NAME2[]        = "number_of_csu_periods_with_more_than_%ddays_per_time_period";
static const char CSU_LONGNAME2[]    = "Number of csu periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CSU_UNITS2[]       = "No.";

static const char CWDI_NAME[]        = "cold_wave_duration_index_wrt_mean_of_reference_period";
static const char CWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily minimum temperature is more than %1.0f degrees below a reference value. The reference value is calculated  as the mean of minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS[]       = "No.";
static const char CWDI_NAME2[]       = "cold_waves_per_time_period";
static const char CWDI_LONGNAME2[]   = "Number of cold waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWDI_UNITS2[]      = "No.";

static const char CWFI_NAME[]        = "cold_spell_days_index_wrt_10th_percentile_of_reference_period";
static const char CWFI_NAME_ET[]     = "csdiETCCDI";
static const char CWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is below a reference value. The reference value is calculated  as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_LONGNAME_ET[] = "Cold Spell Duration Index";
static const char CWFI_UNITS[]       = "No.";
static const char CWFI_UNITS_ET[]    = "days";
static const char CWFI_NAME2[]       = "cold_spell_periods_per_time_period";
static const char CWFI_LONGNAME2[]   = "Number of cold spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char CWFI_UNITS2[]      = "No.";

static const char ETR_NAME[]         = "intra_period_extreme_temperature_range";
static const char ETR_LONGNAME[]     = "Difference between the absolute extreme temperatures in observation period. The time period should be defined by the bounds of the time coordinate.";
//static const char ETR_UNITS[]        = "K";

static const char FD_NAME[]          = "frost_days_index_per_time_period";
static const char FD_NAME_ET[]       = "fdETCCDI";
static const char FD_LONGNAME[]      = "Frost days index is the number of days where minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char FD_LONGNAME_ET[]   = "Number of Frost Days";
//static const char FD_UNITS[]         = "No.";
static const char FD_UNITS_ET[]      = "days";

static const char GSL_NAME[]         = "thermal_growing_season_length";
static const char GSL_LONGNAME[]     = "Counted are the number of days per calendar year between the first occurrence of at least %d consecutive days where the daily mean temperature is above %1.0f degree Celsius and the first occurrence of at least %d consecutive days after 1st of July where the daily mean temperature is below %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS[]        = "No.";
static const char GSL_NAME2[]        = "day_of_year_of_growing_season_start";
static const char GSL_LONGNAME2[]    = "Day of year of growing season start. The time period should be defined by the bounds of the time coordinate.";
static const char GSL_UNITS2[]       = "No.";

static const char HD_NAME[]          = "heating_degree_days_per_time_period";
static const char HD_LONGNAME[]      = "Heating degree days relates the outside temperature with the room temperature during the heating period. It is the sum of the difference between room temperature X and daily mean temperature Y on days where Y is below a given constant A. X is 20 degree Celsius and A is 15 degree Celsius according to VDI guidelines. According to ECAD both X and A are 17 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char HD_UNITS[]         = "No.";

static const char HWDI_NAME[]        = "heat_wave_duration_index_wrt_mean_of_reference_period";
static const char HWDI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily maximum temperature is more than %1.0f degrees above a reference value. The reference value is calculated  as the mean of maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS[]       = "No.";
static const char HWDI_NAME2[]       = "heat_waves_per_time_period";
static const char HWDI_LONGNAME2[]   = "Number of heat waves per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWDI_UNITS2[]      = "No.";

static const char HWFI_NAME[]        = "warm_spell_days_index_wrt_90th_percentile_of_reference_period";
static const char HWFI_NAME_ET[]     = "wsdiETCCDI";
static const char HWFI_LONGNAME[]    = "This is the number of days per time period where in intervals of at least %d consecutive days the daily mean temperature is above a reference value. The reference value is calculated  as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_LONGNAME_ET[] = "Warm Spell Duration Index";
static const char HWFI_UNITS[]       = "No.";
static const char HWFI_UNITS_ET[]    = "days";
static const char HWFI_NAME2[]       = "warm_spell_periods_per_time_period";
static const char HWFI_LONGNAME2[]   = "Number of warm spell periods per time period. The time period should be defined by the bounds of the time coordinate.";
static const char HWFI_UNITS2[]      = "No.";

static const char ID_NAME[]          = "ice_days_index_per_time_period";
static const char ID_NAME_ET[]       = "idETCCDI";
static const char ID_LONGNAME[]      = "Ice days index is the number of days where maximum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char ID_LONGNAME_ET[]   = "Number of Icing Days";
static const char ID_UNITS[]         = "No.";
static const char ID_UNITS_ET[]      = "days";

static const char SU_NAME[]          = "summer_days_index_per_time_period";
static const char SU_NAME_ET[]       = "suETCCDI";
static const char SU_LONGNAME[]      = "Summer days index is the number of days where maximum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char SU_LONGNAME_ET[]   = "Number of Summer Days";
//static const char SU_UNITS[]         = "No.";
static const char SU_UNITS_ET[]      = "days";

static const char TG10P_NAME[]       = "cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TG10P_LONGNAME[]   = "This is the percent of time per time period where daily mean temperature is below a reference value. The reference value is calculated as the 10th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG10P_UNITS[]      = "Percent";

static const char TG90P_NAME[]       = "warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TG90P_LONGNAME[]   = "This is the percent of time per time period where daily mean temperature is above a reference value. The reference value is calculated as the 90th percentile of daily mean temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TG90P_UNITS[]      = "Percent";

static const char TN10P_NAME[]       = "cold_nights_percent_wrt_10th_percentile_of_reference_period";
static const char TN10P_LONGNAME[]   = "This is the percent of time per time period where daily minimum temperature is below a reference value. The reference value is calculated as the 10th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN10P_UNITS[]      = "Percent";

static const char TN90P_NAME[]       = "warm_nights_percent_wrt_90th_percentile_of_reference_period";
static const char TN90P_LONGNAME[]   = "This is the percent of time per time period where daily minimum temperature is above a reference value. The reference value is calculated as the 90th percentile of daily minimum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TN90P_UNITS[]      = "Percent";

static const char TR_NAME[]          = "tropical_nights_index_per_time_period";
static const char TR_NAME_ET[]       = "trETCCDI";
static const char TR_LONGNAME[]      = "Tropical nights index is the number of days where minimum of temperature is above %1.0f degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char TR_LONGNAME_ET[]   = "Number of Tropical Nights";
static const char TR_UNITS[]         = "No.";
static const char TR_UNITS_ET[]      = "days";

static const char TX10P_NAME[]       = "very_cold_days_percent_wrt_10th_percentile_of_reference_period";
static const char TX10P_LONGNAME[]   = "This is the percent of time per time period where daily maximum temperature is below a reference value. The reference value is calculated as the 10th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX10P_UNITS[]      = "Percent";

static const char TX90P_NAME[]       = "very_warm_days_percent_wrt_90th_percentile_of_reference_period";
static const char TX90P_LONGNAME[]   = "This is the percent of time per time period where daily maximum temperature is above a reference value. The reference value is calculated as the 90th percentile of daily maximum temperatures of a five day window centred on each calendar day of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char TX90P_UNITS[]      = "Percent";

static const char CDD_NAME[]         = "consecutive_dry_days_index_per_time_period";
static const char CDD_NAME_ET[]      = "cddETCCDI";
static const char CDD_LONGNAME[]     = "Consecutive dry days is the greatest number of consecutive days per time period with daily precipitation amount below %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_LONGNAME_ET[]  = "Maximum Number of Consecutive Days with Less Than 1mm of Precipitation [days]";
static const char CDD_UNITS[]        = "No.";
static const char CDD_UNITS_ET[]      = "days";
static const char CDD_NAME2[]        = "number_of_cdd_periods_with_more_than_%ddays_per_time_period";
static const char CDD_LONGNAME2[]    = "Number of cdd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CDD_UNITS2[]       = "No.";

static const char CWD_NAME[]         = "consecutive_wet_days_index_per_time_period";
static const char CWD_NAME_ET[]      = "cwdETCCDI";
static const char CWD_LONGNAME[]     = "Consecutive wet days is the greatest number of consecutive days per time period with daily precipitation above %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_LONGNAME_ET[]  = "Maximum Number of Consecutive Days with At Least 1mm of Precipitation";
static const char CWD_UNITS[]        = "No.";
static const char CWD_UNITS_ET[]     = "days";
static const char CWD_NAME2[]        = "number_of_cwd_periods_with_more_than_%ddays_per_time_period";
static const char CWD_LONGNAME2[]    = "Number of cwd periods in given time period with more than %d days. The time period should be defined by the bounds of the time coordinate.";
static const char CWD_UNITS2[]       = "No.";

static const char PD_NAME[]          = "precipitation_days_index_per_time_period";
static const char PD_NAME_ET[]       = "r1mmETCCDI";
static const char PD_LONGNAME[]      = "precipitation days is the number of days per time period with daily precipitation sum exceeding %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char PD_LONGNAME_ET[]   = "Count of Days with At Least 1mm of Precipitation";
static const char PD_UNITS[]         = "No.";
static const char PD_UNITS_ET[]      = "days";

static const char R10MM_NAME[]       = "heavy_precipitation_days_index_per_time_period";
static const char R10MM_NAME_ET[]    = "r10mmETCCDI";
static const char R10MM_LONGNAME[]   = "Heavy precipitation days is the number of days per time period with daily precipitation sum exceeding 10 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R10MM_LONGNAME_ET[]= "Count of Days with At Least 10mm of Precipitation";
static const char R10MM_UNITS[]      = "No.";
static const char R10MM_UNITS_ET[]   = "days";

static const char R20MM_NAME[]       = "very_heavy_precipitation_days_index_per_time_period";
static const char R20MM_NAME_ET[]    = "r20mmETCCDI";
static const char R20MM_LONGNAME[]   = "Very heavy precipitation days is the number of days with daily precipitation sum exceeding 20 mm. The time period should be defined by the bounds of the time coordinate.";
static const char R20MM_LONGNAME_ET[]= "Count of Days with At Least 20mm of Precipitation";
static const char R20MM_UNITS[]      = "No.";
static const char R20MM_UNITS_ET[]   = "days";

static const char R75P_NAME[]        = "moderate_wet_days_wrt_75th_percentile_of_reference_period";
static const char R75P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 75th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R75P_UNITS[]       = "Percent";

static const char R75PTOT_NAME[]     = "precipitation_percent_due_to_R75p_days";
static const char R75PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due to moderate_wet_days_wrt_75th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R75PTOT_UNITS[]    = "Percent";

static const char R90P_NAME[]        = "wet_days_wrt_90th_percentile_of_reference_period";
static const char R90P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 90th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R90P_UNITS[]       = "Percent";

static const char R90PTOT_NAME[]     = "precipitation_percent_due_to_R90p_days";
static const char R90PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due towet_days_wrt_90th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R90PTOT_UNITS[]    = "Percent";

static const char R95P_NAME[]        = "very_wet_days_wrt_95th_percentile_of_reference_period";
static const char R95P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 95th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R95P_UNITS[]       = "Percent";

static const char R95PTOT_NAME[]     = "precipitation_percent_due_to_R95p_days";
static const char R95PTOT_LONGNAME[] = "Percentage of total precipitation amount per time period due to very_wet_days_wrt_95th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
static const char R95PTOT_UNITS[]    = "Percent";

static const char R99P_NAME[]        = "extremely_wet_days_wrt_99th_percentile_of_reference_period";
static const char R99P_LONGNAME[]    = "This is the percent of time per time period of wet days (daily sum at least 1 mm / day) where daily precipitation amount of a wet day is above a reference value. The reference value is calculated as the 99th percentile of all wet days of a given 30 year climate reference period. The time period should be defined by the bounds of the time coordinate.";
static const char R99P_UNITS[]       = "Percent";

static const char R99PTOT_NAME[]     = "precipitation_percent_due_to_R99p_days";
static const char R99PTOT_LONGNAME[] = "percentage of total  precipitation amount per time period due to extremely_wet_days_wrt_99th_percentile_of_reference_period. The time period should be defined by the bounds of the time coordinate.";
//static const char R99PTOT_UNITS[]    = "Percent";

static const char RR1_NAME[]         = "wet_days_index_per_time_period";
static const char RR1_LONGNAME[]     = "Wet days index is the number of days per time period with daily precipitation of at least %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char RR1_UNITS[]        = "No.";

static const char RX1DAY_NAME[]      = "highest_one_day_precipitation_amount_per_time_period";
static const char RX1DAY_NAME_ET[]   = "rx1dayETCCDI";
static const char RX1DAY_LONGNAME[]  = "Highest one day precipitation is the maximum of one day precipitation amount in a given time period. The time period should be defined by the bounds of the time coordinate.";
static const char RX1DAY_LONGNAME_ET[]= "Maximum 1-day Precipitation";
static const char RX1DAY_UNITS[]     = "mm per day";
static const char RX1DAY_UNITS_ET[]  = "mm";

static const char RX5DAY_NAME[]      = "highest_five_day_precipitation_amount_per_time_period";
static const char RX5DAY_NAME_ET[]   = "rx5dayETCCDI";
static const char RX5DAY_LONGNAME[]  = "Highest precipitation amount for five day interval (including the calendar day as the last day). The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_LONGNAME_ET[]= "Maximum Consecutive 5-day Precipitation";
static const char RX5DAY_UNITS[]     = "mm per 5 day";
static const char RX5DAY_UNITS_ET[]  = "mm";
static const char RX5DAY_NAME2[]     = "number_of_5day_heavy_precipitation_periods_per_time_period";
static const char RX5DAY_LONGNAME2[] = "Number of 5day periods in given time period with precipitation amount exceeding %1.0f mm / 5 days. The time period should be defined by the bounds of the time coordinate.";
static const char RX5DAY_UNITS2[]    = "No.";

static const char SDII_NAME[]        = "simple_daily_intensity_index_per_time_period";
static const char SDII_NAME_ET[]     = "sdiiETCCDI";
static const char SDII_LONGNAME[]    = "Simple daily intensity index is the mean of precipitation amount on wet days. A wet day is a day with precipitation sum of at least %g mm. The time period should be defined by the bounds of the time coordinate.";
static const char SDII_LONGNAME_ET[] = "Simple Precipitation Intensity Index";
static const char SDII_UNITS[]       = "mm";
static const char SDII_UNITS_ET[]    = "mm d-1";

static const char FDNS_NAME[]        = "frost_days_where_no_snow_index_per_time_period";
static const char FDNS_LONGNAME[]    = "Frost days where no snow index is the number of days without snowcover and where the minimum of temperature is below 0 degree Celsius. The time period should be defined by the bounds of the time coordinate.";
static const char FDNS_UNITS[]       = "No.";

static const char STRWIN_NAME[]      = "strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME[]  = "Strong wind days index is the number of days per time period where maximum wind speed is above %1.0f m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS[]     = "No.";
static const char STRWIN_NAME2[]     = "consecutive_strong_wind_days_index_per_time_period";
static const char STRWIN_LONGNAME2[] = "Greatest number of consecutive strong wind days per time period. The time period should be defined by the bounds of the time coordinate.";
static const char STRWIN_UNITS2[]    = "No.";

static const char STRBRE_NAME[]      = "strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME[]  = "Strong breeze days index is the number of days per time period where maximum wind speed is above 10.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char STRBRE_NAME2[]     = "consecutive_strong_breeze_days_index_per_time_period";
static const char STRBRE_LONGNAME2[] = "Greatest number of consecutive strong breeze days per time period. The time period should be defined by the bounds of the time coordinate.";

//static const char STRGAL_NAME[]      = "strong_gale_days_index_per_time_period";
//static const char STRGAL_LONGNAME[]  = "Strong gale days index is the number of days per time period where maximum wind speed is above 20.5 m/s. The time period should be defined by the bounds of the time coordinate.";
//static const char STRGAL_NAME2[]     = "consecutive_strong_gale_days_index_per_time_period";
//static const char STRGAL_LONGNAME2[] = "Greatest number of consecutive strong gale days per time period. The time period should be defined by the bounds of the time coordinate.";

static const char HURR_NAME[]        = "hurricane_days_index_per_time_period";
static const char HURR_LONGNAME[]    = "Hurricane days index is the number of days per time period where maximum wind speed is above 32.5 m/s. The time period should be defined by the bounds of the time coordinate.";
static const char HURR_NAME2[]       = "consecutive_hurricane_days_index_per_time_period";
static const char HURR_LONGNAME2[]   = "Greatest number of consecutive hurricane days per time period. The time period should be defined by the bounds of the time coordinate.";

// clang-format on

/* ECA temperature indices */

static void
set_default_compare_type(int &compare_type)
{
  compare_type = cdo_operator_f2(cdo_operator_id());
}

static void
set_compare_type_from_params(int &compare_type, std::vector<std::string> const &params)
{
  KVList kvlist;
  if (kvlist.parse_arguments(params) != 0) cdo_abort("Argument parse error!");
  auto kv = kvlist.search("freq");
  if (kv && kv->nvalues > 0)
  {
    // clang-format off
      if      (kv->values[0] == "month") compare_type = CMP_MONTH;
      else if (kv->values[0] == "year")  compare_type  = CMP_YEAR;
      else cdo_abort("Frequency '%s' unknown.", kv->values[0]);
    // clang-format on
  }
}

#include <functional>

template <typename Request>
class EcaIndices : public Process
{
protected:
  EcaIndices(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module,
             std::function<void(Request p_request)> p_eca_func)
      : Process(p_ID, p_operatorName, p_arguments, p_module), eca_func(p_eca_func)
  {
  }
  std::function<void(Request p_request)> eca_func;
  Request request;
  // virtual void init() = 0;

public:
  void
  run() override
  {
    assert(request.compare_type != -1);
    eca_func(request);
  }

  void
  close() override
  {
  }
};

class EcaIndices1 : public EcaIndices<ECA_REQUEST_1>
{
public:
  EcaIndices1(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module)
      : EcaIndices<ECA_REQUEST_1>(p_ID, p_operatorName, p_arguments, p_module, eca1)
  {
  }
};

class EcaIndices2 : public EcaIndices<ECA_REQUEST_2>
{
public:
  EcaIndices2(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module)
      : EcaIndices<ECA_REQUEST_2>(p_ID, p_operatorName, p_arguments, p_module, eca2)
  {
  }
};

class EcaIndices3 : public EcaIndices<ECA_REQUEST_3>
{
public:
  EcaIndices3(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module)
      : EcaIndices<ECA_REQUEST_3>(p_ID, p_operatorName, p_arguments, p_module, eca3)
  {
  }
};

class EcaIndices4 : public EcaIndices<ECA_REQUEST_4>
{
public:
  EcaIndices4(int p_ID, std::string const &p_operatorName, std::vector<std::string> const &p_arguments, const CdoModule &p_module)
      : EcaIndices<ECA_REQUEST_4>(p_ID, p_operatorName, p_arguments, p_module, eca4)
  {
  }
};

class EcaCfd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaCfd",
    .operators = { { "eca_cfd", 0, CMP_DATE, Eca_cfdHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCfd> registration = RegisterEntry<EcaCfd>();
  int ndays = 5;
  char cfd_name2[1024];
  char cfd_longname2[1024];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else
    {
      if (cdo_operator_argc() > 0) ndays = parameter_to_int(cdo_operator_argv(0));
    }

    std::snprintf(cfd_name2, sizeof(cfd_name2), CFD_NAME2, ndays);
    std::snprintf(cfd_longname2, sizeof(cfd_longname2), CFD_LONGNAME2, ndays);

    request.var1.name = CFD_NAME;
    request.var1.longname = CFD_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.f1 = vfarselltc;
    request.var1.f1arg = TO_KELVIN(0.0);
    request.var1.f2 = vfarnum2;
    request.var1.f3 = field2_max;
    request.var2.name = cfd_name2;
    request.var2.longname = cfd_longname2;
    request.var2.units = CFD_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = ndays + 1;
    request.var2.h3 = vfarnum;
  }
};

class EcaCsu : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaCsu",
    .operators = { { "eca_csu", 0, CMP_DATE, Eca_csuHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCsu> registration = RegisterEntry<EcaCsu>();
  double argT = 25.0;
  int ndays = 5;
  char csu_name2[1024];
  char csu_longname2[1024];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 3) cdo_abort("Too many arguments!");
    if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else if (cdo_operator_argc() > 0)
    {
      argT = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
    }

    std::snprintf(csu_name2, sizeof(csu_name2), CSU_NAME2, ndays);
    std::snprintf(csu_longname2, sizeof(csu_longname2), CSU_LONGNAME2, ndays);

    request.var1.name = CSU_NAME;
    request.var1.longname = CSU_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.f1 = vfarselgtc;
    request.var1.f1arg = TO_KELVIN(argT);
    request.var1.f2 = vfarnum2;
    request.var1.f3 = field2_max;
    request.var2.name = csu_name2;
    request.var2.longname = csu_longname2;
    request.var2.units = CSU_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = ndays + 1;
    request.var2.h3 = vfarnum;
  }
};

class EcaCwdi : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaCwdi",
    .operators = { { "eca_cwdi", 0, CMP_DATE, Eca_cwdiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCwdi> registration = RegisterEntry<EcaCwdi>();
  int argN = 6;
  double argT = 5.0;
  char longname[sizeof(CWDI_LONGNAME) + 80];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      set_compare_type_from_params(request.compare_type, params);
      argT = parameter_to_double(cdo_operator_argv(1));
      argN = parameter_to_int(cdo_operator_argv(0));
    }
    else
    {
      if (cdo_operator_argc() > 1)
        argT = parameter_to_double(cdo_operator_argv(1));
      else if (cdo_operator_argc() > 0)
        argN = parameter_to_int(cdo_operator_argv(0));
    }

    std::snprintf(longname, sizeof(longname), CWDI_LONGNAME, argN, argT);

    request.var1.name = CWDI_NAME;
    request.var1.longname = longname;
    request.var1.refdate = ECA_refdate;
    request.var1.units = CWDI_UNITS;
    request.var1.f2 = fieldc_sub;
    request.var1.f2arg = argT;
    request.var1.f3 = vfarsellt;
    request.var1.f4 = vfarnum2;
    request.var1.f5 = vfarnum3;
    request.var1.f5arg = argN;
    request.var2.name = CWDI_NAME2;
    request.var2.longname = CWDI_LONGNAME2;
    request.var2.units = CWDI_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = argN;
    request.var2.h2 = vfarnum;
  }
};

class EcaCwfi : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaCwfi",
    .operators = { { "eca_cwfi", 0, CMP_DATE, Eca_cwfiHelp }, { "etccdi_csdi", 0, CMP_YEAR, Eca_cwfiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCwfi> registration = RegisterEntry<EcaCwfi>();
  int ECA_CWFI, ETCCDI_CSDI;
  int argN = 6;
  char longname[sizeof(CWFI_LONGNAME) + 40];

public:
  void
  init() override
  {
    ECA_CWFI = module.get_id("eca_cwfi");
    ETCCDI_CSDI = module.get_id("etccdi_csdi");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
      argN = parameter_to_int(cdo_operator_argv(0));
    }
    else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
    }

    if (ECA_CWFI == cdo_operator_id())
    {
      std::snprintf(longname, sizeof(longname), CWFI_LONGNAME, argN);
      request.var1.name = CWFI_NAME;
      request.var1.longname = longname;
      request.var1.units = CWFI_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_CSDI == cdo_operator_id())
    {
      request.var1.name = CWFI_NAME_ET;
      request.var1.longname = CWFI_LONGNAME_ET;
      request.var1.units = CWFI_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f3 = vfarsellt;
    request.var1.f4 = vfarnum2;
    request.var1.f5 = vfarnum3;
    request.var1.f5arg = argN;
    request.var2.name = CWFI_NAME2;
    request.var2.longname = CWFI_LONGNAME2;
    request.var2.units = CWFI_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = argN;
    request.var2.h2 = vfarnum;
  }
};

class EcaEtr : public EcaIndices3
{
public:
  using EcaIndices3::EcaIndices3;
  inline static CdoModule module = {
    .name = "EcaEtr",
    .operators = { { "eca_etr", 0, CMP_DATE, Eca_etrHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaEtr> registration = RegisterEntry<EcaEtr>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.name = ETR_NAME;
    request.longname = ETR_LONGNAME;
    request.refdate = ECA_refdate;
    request.f1 = field2_max;
    request.f2 = field2_min;
    request.f3 = field2_sub;
  }
};

class EcaFd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaFd",
    .operators = { { "eca_fd", 0, CMP_DATE, Eca_fdHelp }, { "etccdi_fd", 0, CMP_YEAR, Eca_fdHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaFd> registration = RegisterEntry<EcaFd>();
  int ECA_FD, ETCCDI_FD;

public:
  void
  init() override
  {
    ECA_FD = module.get_id("eca_fd");
    ETCCDI_FD = module.get_id("etccdi_fd");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }
    else {}

    if (ECA_FD == cdo_operator_id())
    {
      request.var1.name = FD_NAME;
      request.var1.longname = FD_LONGNAME;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_FD == cdo_operator_id())
    {
      request.var1.name = FD_NAME_ET;
      request.var1.longname = FD_LONGNAME_ET;
      request.var1.units = FD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselltc;
    request.var1.f1arg = TO_KELVIN(0.0);
    request.var1.f2 = vfarnum;
  }
};

/*
 * Definition of GSL: (Thermal) Growing Season Length start at the first span
 * of at least 6 (argN) days with T > 5.0°C (argT) in first half of the year
 * and ends at the first span of ar least 6 (argN) days with T < 5.0°C (argT)
 * in the second half.
 * ATTENTION: Year of the northern hemisphere starts in january to
 * december, whereas for the southern hemisphere is goes from july to june!
 * Hence, at least 18 Month of data is needed for computing the gsl of the whole earth.
 */
class EcaGsl : public EcaIndices4
{
public:
  using EcaIndices4::EcaIndices4;
  inline static CdoModule module = {
    .name = "EcaGsl",
    .operators = { { "eca_gsl", 0, CMP_YEAR, Eca_gslHelp }, { "etccdi_gsl", Eca_gslHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaGsl> registration = RegisterEntry<EcaGsl>();
  int argN = 6;
  double argT = 5.0;
  double minLandFraction = 0.5;
  char longname[sizeof(GSL_LONGNAME) + 160];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
    if (cdo_operator_argc() > 1) argT = parameter_to_double(cdo_operator_argv(1));
    if (cdo_operator_argc() > 2) minLandFraction = parameter_to_double(cdo_operator_argv(2));

    std::snprintf(longname, sizeof(longname), GSL_LONGNAME, argN, argT, argN, argT);

    request.name = GSL_NAME;
    request.longname = longname;
    request.units = GSL_UNITS;
    request.name2 = GSL_NAME2;
    request.longname2 = GSL_LONGNAME2;
    request.units2 = GSL_UNITS2;
    request.s1 = vfarselgtc;
    request.s1arg = TO_KELVIN(argT);
    request.s2 = vfarselltc;
    request.s2arg = TO_KELVIN(argT);
    request.s3 = vfarselgec;
    request.s3arg = minLandFraction;
    request.consecutiveDays = argN;
  }
};

class EcaHd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaHd",
    .operators = { { "eca_hd", 0, CMP_DATE, Eca_hdHelp }, { "etccdi_hd", Eca_hdHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaHd> registration = RegisterEntry<EcaHd>();
  double argX = 17.0;
  double argA = 17.0;

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      argX = parameter_to_double(cdo_operator_argv(0));
      argA = argX;
    }
    if (cdo_operator_argc() > 1) argA = parameter_to_double(cdo_operator_argv(1));

    request.var1.name = HD_NAME;
    request.var1.longname = HD_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = HD_UNITS;
    request.var1.f1 = vfarselltc;
    request.var1.f1arg = TO_KELVIN(argA);
    request.var1.f2 = field2_sum;
    request.var1.mulc = -1.0;
    request.var1.addc = TO_KELVIN(argX);
  }
};

class EcaHwdi : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaHwdi",
    .operators = { { "eca_hwdi", 0, CMP_DATE, Eca_hwdiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaHwdi> registration = RegisterEntry<EcaHwdi>();
  int argN = 6;
  double argT = 5.0;
  char longname[sizeof(HWDI_LONGNAME) + 80];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      set_compare_type_from_params(request.compare_type, params);
      argN = parameter_to_int(cdo_operator_argv(0));
      argT = parameter_to_double(cdo_operator_argv(1));
    }
    else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
      if (cdo_operator_argc() > 1) argT = parameter_to_double(cdo_operator_argv(1));
    }

    std::snprintf(longname, sizeof(longname), HWDI_LONGNAME, argN, argT);

    request.var1.name = HWDI_NAME;
    request.var1.longname = longname;
    request.var1.refdate = ECA_refdate;
    request.var1.units = HWDI_UNITS;
    request.var1.f2 = fieldc_add;
    request.var1.f2arg = argT;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum2;
    request.var1.f5 = vfarnum3;
    request.var1.f5arg = argN;
    request.var2.name = HWDI_NAME2;
    request.var2.longname = HWDI_LONGNAME2;
    request.var2.units = HWDI_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = argN;
    request.var2.h2 = vfarnum;
  }
};

class EcaHwfi : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaHwfi",
    .operators = { { "eca_hwfi", 0, CMP_DATE, Eca_hwfiHelp }, { "etccdi_wsdi", 0, CMP_YEAR, Eca_hwfiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaHwfi> registration = RegisterEntry<EcaHwfi>();
  int ECA_HWFI, ETCCDI_WSDI;
  int argN = 6;
  char longname[sizeof(HWFI_LONGNAME) + 40];

public:
  void
  init() override
  {
    ECA_HWFI = module.get_id("eca_hwfi");
    ETCCDI_WSDI = module.get_id("etccdi_wsdi");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
      argN = parameter_to_int(cdo_operator_argv(0));
    }
    else
    {
      if (cdo_operator_argc() > 0) argN = parameter_to_int(cdo_operator_argv(0));
    }

    if (ECA_HWFI == cdo_operator_id())
    {
      std::snprintf(longname, sizeof(longname), HWFI_LONGNAME, argN);
      request.var1.name = HWFI_NAME;
      request.var1.longname = longname;
      request.var1.units = HWFI_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_WSDI == cdo_operator_id())
    {
      request.var1.name = HWFI_NAME_ET;
      request.var1.longname = HWFI_LONGNAME_ET;
      request.var1.units = HWFI_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum2;
    request.var1.f5 = vfarnum3;
    request.var1.f5arg = argN;
    request.var2.name = HWFI_NAME2;
    request.var2.longname = HWFI_LONGNAME2;
    request.var2.units = HWFI_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = argN;
    request.var2.h2 = vfarnum;
  }
};

class EcaId : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaId",
    .operators = { { "eca_id", 0, CMP_DATE, Eca_idHelp }, { "etccdi_id", 0, CMP_YEAR, Eca_idHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaId> registration = RegisterEntry<EcaId>();
  int ECA_ID, ETCCDI_ID;

public:
  void
  init() override
  {

    ECA_ID = module.get_id("eca_id");
    ETCCDI_ID = module.get_id("etccdi_id");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }
    else {}

    if (ECA_ID == cdo_operator_id())
    {
      request.var1.name = ID_NAME;
      request.var1.longname = ID_LONGNAME;
      request.var1.units = ID_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_ID == cdo_operator_id())
    {
      request.var1.name = ID_NAME_ET;
      request.var1.longname = ID_LONGNAME_ET;
      request.var1.units = ID_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselltc;
    request.var1.f1arg = TO_KELVIN(0.0);
    request.var1.f2 = vfarnum;
  }
};

class EcaSu : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaSu",
    .operators = { { "eca_su", 0, CMP_DATE, Eca_suHelp }, { "etccdi_su", 0, CMP_YEAR, Eca_suHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaSu> registration = RegisterEntry<EcaSu>();
  int ECA_SU, ETCCDI_SU;
  double argT = 25.0;
  char longname[sizeof(SU_LONGNAME) + 40];

public:
  void
  init() override
  {

    ECA_SU = module.get_id("eca_su");
    ETCCDI_SU = module.get_id("etccdi_su");
    set_default_compare_type(request.compare_type);
    if (cdo_operator_argc() > 0) argT = parameter_to_double(cdo_operator_argv(0));
    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }

    if (ECA_SU == cdo_operator_id())
    {
      std::snprintf(longname, sizeof(longname), SU_LONGNAME, argT);
      request.var1.name = SU_NAME;
      request.var1.longname = longname;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_SU == cdo_operator_id())
    {
      request.var1.name = SU_NAME_ET;
      request.var1.longname = SU_LONGNAME_ET;
      request.var1.units = SU_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselgtc;
    request.var1.f1arg = TO_KELVIN(argT);
    request.var1.f2 = vfarnum;
  }
};

class EcaTg10p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTg10p",
    .operators = { { "eca_tg10p", 0, CMP_DATE, Eca_tg10pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTg10p> registration = RegisterEntry<EcaTg10p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = TG10P_NAME;
    request.var1.longname = TG10P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TG10P_UNITS;
    request.var1.f3 = vfarsellt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaTg90p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTg90p",
    .operators = { { "eca_tg90p", 0, CMP_DATE, Eca_tg90pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTg90p> registration = RegisterEntry<EcaTg90p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = TG90P_NAME;
    request.var1.longname = TG90P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TG90P_UNITS;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaTn10p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTn10p",
    .operators = { { "eca_tn10p", 0, CMP_DATE, Eca_tn10pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTn10p> registration = RegisterEntry<EcaTn10p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = TN10P_NAME;
    request.var1.longname = TN10P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TN10P_UNITS;
    request.var1.f3 = vfarsellt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaTn90p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTn90p",
    .operators = { { "eca_tn90p", 0, CMP_DATE, Eca_tn90pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTn90p> registration = RegisterEntry<EcaTn90p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = TN90P_NAME;
    request.var1.longname = TN90P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TN90P_UNITS;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaTr : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaTr",
    .operators = { { "eca_tr", 0, CMP_DATE, Eca_trHelp }, { "etccdi_tr", 0, CMP_YEAR, Eca_trHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTr> registration = RegisterEntry<EcaTr>();
  int ECA_TR, ETCCDI_TR;
  double argT = 20.0;
  char longname[1024];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    ECA_TR = module.get_id("eca_tr");
    ETCCDI_TR = module.get_id("etccdi_tr");
    if (cdo_operator_argc() > 0) argT = parameter_to_double(cdo_operator_argv(0));
    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }

    if (ECA_TR == cdo_operator_id())
    {
      std::snprintf(longname, sizeof(longname), TR_LONGNAME, argT);
      request.var1.name = TR_NAME;
      request.var1.longname = longname;
      request.var1.units = TR_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_TR == cdo_operator_id())
    {
      request.var1.name = TR_NAME_ET;
      request.var1.longname = TR_LONGNAME_ET;
      request.var1.units = TR_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselgtc;
    request.var1.f1arg = TO_KELVIN(argT);
    request.var1.f2 = vfarnum;
  }
};

class EcaTx10p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTx10p",
    .operators = { { "eca_tx10p", 0, CMP_DATE, Eca_tx10pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTx10p> registration = RegisterEntry<EcaTx10p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = TX10P_NAME;
    request.var1.longname = TX10P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TX10P_UNITS;
    request.var1.f3 = vfarsellt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaTx90p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaTx90p",
    .operators = { { "eca_tx90p", 0, CMP_DATE, Eca_tx90pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaTx90p> registration = RegisterEntry<EcaTx90p>();

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      if ('m' == cdo_operator_argv(0)[0])
        request.compare_type = CMP_MONTH;
      else
        cdo_warning("Parameter value '%s' is invalid. The only valid value is "
                    "'m' indicating monthly mode. Operating in yearly mode now.",
                    cdo_operator_argv(0));
    }

    request.var1.name = TX90P_NAME;
    request.var1.longname = TX90P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = TX90P_UNITS;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

// ECA precipitation indices

class EcaCdd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaCdd",
    .operators = { { "eca_cdd", 0, CMP_DATE, Eca_cddHelp }, { "etccdi_cdd", 0, CMP_YEAR, Eca_cddHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCdd> registration = RegisterEntry<EcaCdd>();
  int ECA_CDD, ETCCDI_CDD;
  double threshold = 1;
  int ndays = 5;
  char cdd_longname[1024];
  char cdd_longname2[1024];
  char cdd_name2[1024];

public:
  void
  init() override
  {

    ECA_CDD = module.get_id("eca_cdd");
    ETCCDI_CDD = module.get_id("etccdi_cdd");

    set_default_compare_type(request.compare_type);
    if (cdo_operator_argc() > 3)
      cdo_abort("Too many arguments!");
    else if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else if (cdo_operator_argc() > 0)
    {
      threshold = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
    }

    std::snprintf(cdd_longname, sizeof(cdd_longname), CDD_LONGNAME, threshold);
    std::snprintf(cdd_longname2, sizeof(cdd_longname2), CDD_LONGNAME2, ndays);
    std::snprintf(cdd_name2, sizeof(cdd_name2), CDD_NAME2, ndays);

    if (ECA_CDD == cdo_operator_id())
    {
      request.var1.name = CDD_NAME;
      request.var1.longname = cdd_longname;
      request.var1.units = CDD_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_CDD == cdo_operator_id())
    {
      request.var1.name = CDD_NAME_ET;
      request.var1.longname = CDD_LONGNAME_ET;
      request.var1.units = CDD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselltc;
    request.var1.f1arg = threshold;
    request.var1.f2 = vfarnum2;
    request.var1.f3 = field2_max;
    request.var2.name = cdd_name2;
    request.var2.longname = cdd_longname2;
    request.var2.units = CDD_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = ndays + 1;
    request.var2.h3 = vfarnum;
  }
};

class EcaCwd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaCwd",
    .operators = { { "eca_cwd", 0, CMP_DATE, Eca_cwdHelp }, { "etccdi_cwd", 0, CMP_YEAR, Eca_cwdHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaCwd> registration = RegisterEntry<EcaCwd>();

  int ECA_CWD, ETCCDI_CWD;
  double threshold = 1;
  int ndays = 5;
  char cwd_longname[1024];
  char cwd_longname2[1024];
  char cwd_name2[1024];

public:
  void
  init() override
  {
    ECA_CWD = module.get_id("eca_cwd");
    ETCCDI_CWD = module.get_id("etccdi_cwd");

    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 3)
      cdo_abort("Too many arguments!");
    else if (cdo_operator_argc() > 2)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 2, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else if (cdo_operator_argc() > 0)
    {
      threshold = parameter_to_double(cdo_operator_argv(0));
      if (cdo_operator_argc() == 2) ndays = parameter_to_int(cdo_operator_argv(1));
    }

    std::snprintf(cwd_longname, sizeof(cwd_longname), CWD_LONGNAME, threshold);
    std::snprintf(cwd_longname2, sizeof(cwd_longname2), CWD_LONGNAME2, ndays);
    std::snprintf(cwd_name2, sizeof(cwd_name2), CWD_NAME2, ndays);

    if (ECA_CWD == cdo_operator_id())
    {
      request.var1.name = CWD_NAME;
      request.var1.longname = cwd_longname;
      request.var1.units = CWD_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_CWD == cdo_operator_id())
    {
      request.var1.name = CWD_NAME_ET;
      request.var1.longname = CWD_LONGNAME_ET;
      request.var1.units = CWD_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = threshold;
    request.var1.f2 = vfarnum2;
    request.var1.f3 = field2_max;
    request.var2.name = cwd_name2;
    request.var2.longname = cwd_longname2;
    request.var2.units = CWD_UNITS2;
    request.var2.h1 = vfarseleqc;
    request.var2.h1arg = ndays + 1;
    request.var2.h3 = vfarnum;
  }
};

class EcaPd : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaPd",
    .operators = { { "eca_pd", 0, CMP_DATE, Eca_pdHelp },
                   { "eca_r10mm", 0, CMP_DATE, Eca_pdHelp },
                   { "eca_r20mm", 0, CMP_DATE, Eca_pdHelp },
                   { "etccdi_r1mm", 0, CMP_DATE, Eca_pdHelp },
                   { "etccdi_r10mm", 0, CMP_DATE, Eca_pdHelp },
                   { "etccdi_r20mm", 0, CMP_DATE, Eca_pdHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaPd> registration = RegisterEntry<EcaPd>();
  int ECA_PD, ETCCDI_R1MM, ECA_R10MM, ETCCDI_R10MM, ECA_R20MM, ETCCDI_R20MM;
  char lnamebuffer[1024];
  double threshold = 0;

public:
  void
  init() override
  {
    ECA_PD = module.get_id("eca_pd");
    ETCCDI_R1MM = module.get_id("etccdi_r1mm");
    ECA_R10MM = module.get_id("eca_r10mm");
    ETCCDI_R10MM = module.get_id("etccdi_r10mm");
    ECA_R20MM = module.get_id("eca_r20mm");
    ETCCDI_R20MM = module.get_id("etccdi_r20mm");

    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto params = cdo_get_oper_argv();
      KVList kvlist;
      if (std::strstr(cdo_operator_argv(0).c_str(), "=") || cdo_operator_argc() > 1)
      {
        if (cdo_operator_argc() > 1) params = std::vector<std::string>(params.begin() + 1, params.end());
        if (kvlist.parse_arguments(params) != 0) cdo_abort("Argument parse error!");
        auto kv = kvlist.search("freq");
        if (kv && kv->nvalues > 0)
        {
          // clang-format off
                if      (kv->values[0] == "month") request.compare_type = CMP_MONTH;
                else if (kv->values[0] == "year")  request.compare_type = CMP_YEAR;
          // clang-format on
        }
      }
    }

    auto operatorID = cdo_operator_id();

    if (operatorID == ECA_PD || operatorID == ETCCDI_R1MM)
    {
      if (operatorID == ECA_PD)
      {
        operator_input_arg("daily precipitation amount threshold in [mm]");

        if (cdo_operator_argc() < 1) cdo_abort("Too few arguments!");
        if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");
        threshold = parameter_to_double(cdo_operator_argv(0));
        std::snprintf(lnamebuffer, sizeof(lnamebuffer), PD_LONGNAME, threshold);
        request.var1.name = PD_NAME;
        request.var1.longname = lnamebuffer;
        request.var1.units = PD_UNITS;
      }
      else
      {
        threshold = 1;
        request.var1.name = PD_NAME_ET;
        request.var1.longname = PD_LONGNAME_ET;
        request.var1.units = PD_UNITS_ET;
      }
      if (threshold < 0) cdo_abort("Parameter out of range: threshold = %g", threshold);
    }
    else if (operatorID == ECA_R10MM || operatorID == ETCCDI_R10MM)
    {
      threshold = 10;
      if (operatorID == ECA_R10MM)
      {
        request.var1.name = R10MM_NAME;
        request.var1.longname = R10MM_LONGNAME;
        request.var1.units = R10MM_UNITS;
        request.var1.refdate = ECA_refdate;
      }
      else
      {
        request.var1.name = R10MM_NAME_ET;
        request.var1.longname = R10MM_LONGNAME_ET;
        request.var1.units = R10MM_UNITS_ET;
        request.var1.refdate = ETC_refdate;
      }
    }
    else if (operatorID == ECA_R20MM || operatorID == ETCCDI_R20MM)
    {
      threshold = 20;
      if (operatorID == ECA_R20MM)
      {
        request.var1.name = R20MM_NAME;
        request.var1.longname = R20MM_LONGNAME;
        request.var1.units = R20MM_UNITS;
        request.var1.refdate = ECA_refdate;
      }
      else
      {
        request.var1.name = R20MM_NAME_ET;
        request.var1.longname = R20MM_LONGNAME_ET;
        request.var1.units = R20MM_UNITS_ET;
        request.var1.refdate = ETC_refdate;
      }
    }

    if (Options::cdoVerbose) cdo_print("threshold = %g", threshold);

    request.var1.f1 = vfarselgec;
    request.var1.f1arg = threshold;
    request.var1.f2 = vfarnum;
  }
};

class EcaR75p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR75p",
    .operators = { { "eca_r75p", 0, CMP_DATE, Eca_r75pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR75p> registration = RegisterEntry<EcaR75p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R75P_NAME;
    request.var1.longname = R75P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R75P_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaR75ptot : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR75ptot",
    .operators = { { "eca_r75ptot", 0, CMP_DATE, Eca_r75ptotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR75ptot> registration = RegisterEntry<EcaR75ptot>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R75PTOT_NAME;
    request.var1.longname = R75PTOT_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R75PTOT_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = field2_sum;
    request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;
  }
};

class EcaR90p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR90p",
    .operators = { { "eca_r90p", 0, CMP_DATE, Eca_r90pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR90p> registration = RegisterEntry<EcaR90p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R90P_NAME;
    request.var1.longname = R90P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R90P_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaR95p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR95p",
    .operators = { { "eca_r95p", 0, CMP_DATE, Eca_r95pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR95p> registration = RegisterEntry<EcaR95p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R95P_NAME;
    request.var1.longname = R95P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R95P_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaR90ptot : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR90ptot",
    .operators = { { "eca_r90ptot", 0, CMP_DATE, Eca_r90ptotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR90ptot> registration = RegisterEntry<EcaR90ptot>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R90PTOT_NAME;
    request.var1.longname = R90PTOT_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R90PTOT_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = field2_sum;
    request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;
  }
};

// class ModuleEcaR96p : public EcaIndices2
//{
// public:
//   void
//   init() override
//   {
// //
//     cdo_operator_add("eca_r95p", 0, CMP_DATE, nullptr);
//     set_default_compare_type(request.compare_type);
//
//     request.var1.name = R95P_NAME;
//     request.var1.longname = R95P_LONGNAME;
//     request.var1.refdate = ECA_refdate;
//     request.var1.units = R95P_UNITS;
//     request.var1.f1 = vfarselgec;
//     request.var1.f3 = vfarselgt;
//     request.var1.f4 = vfarnum;
//     request.var1.epilog = PERCENT_OF_TIME;
//   }
// };
//
// void *
// EcaR96p(void *process)
//{
//   ModuleEcaR96p ecaR96p;
//   ecaR96p.init(process);
//   ecaR96p.run();
//   ecaR96p.close();
//   return nullptr;
// }

class EcaR95ptot : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR95ptot",
    .operators = { { "eca_r95ptot", 0, CMP_DATE, Eca_r95ptotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR95ptot> registration = RegisterEntry<EcaR95ptot>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R95PTOT_NAME;
    request.var1.longname = R95PTOT_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R95PTOT_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = field2_sum;
    request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;
  }
};

class EcaR99p : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR99p",
    .operators = { { "eca_r99p", 0, CMP_DATE, Eca_r99pHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR99p> registration = RegisterEntry<EcaR99p>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R99P_NAME;
    request.var1.longname = R99P_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = R99P_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = vfarnum;
    request.var1.epilog = PERCENT_OF_TIME;
  }
};

class EcaR99ptot : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "EcaR99ptot",
    .operators = { { "eca_r99ptot", 0, CMP_DATE, Eca_r99ptotHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaR99ptot> registration = RegisterEntry<EcaR99ptot>();

public:
  void
  init() override
  {

    set_default_compare_type(request.compare_type);

    request.var1.name = R99PTOT_NAME;
    request.var1.longname = R99PTOT_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.f1 = vfarselgec;
    request.var1.f3 = vfarselgt;
    request.var1.f4 = field2_sum;
    request.var1.epilog = PERCENT_OF_TOTAL_AMOUNT;
  }
};

class EcaRr1 : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaRr1",
    .operators = { { "eca_rr1", 0, CMP_DATE, Eca_rr1Help } },
    .aliases = { { "eca_r1mm", "eca_rr1" } },
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaRr1> registration = RegisterEntry<EcaRr1>();
  char longname[1024];

public:
  void
  init() override
  {
    double threshold = 1;

    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2)
      cdo_abort("Too many arguments!");
    else if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else
    {
      if (cdo_operator_argc() == 1) threshold = parameter_to_double(cdo_operator_argv(0));
    }

    std::snprintf(longname, sizeof(longname), RR1_LONGNAME, threshold);

    request.var1.name = RR1_NAME;
    request.var1.longname = longname;
    request.var1.units = RR1_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = threshold;
    request.var1.f2 = vfarnum;
  }
};

class EcaRx1day : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaRx1day",
    .operators = { { "eca_rx1day", 0, CMP_DATE, Eca_rx1dayHelp },
                   { "etccdi_rx1day", 0, CMP_YEAR, Eca_rx1dayHelp },
                   { "etccdi_rx1daymon", Eca_rx1dayHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaRx1day> registration = RegisterEntry<EcaRx1day>();
  int ECA_RX1DAY, ETCCDI_RX1DAY;

public:
  void
  init() override
  {

    ECA_RX1DAY = module.get_id("eca_rx1day");
    ETCCDI_RX1DAY = module.get_id("etccdi_rx1day");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }

    if (ECA_RX1DAY == cdo_operator_id())
    {
      request.var1.name = RX1DAY_NAME;
      request.var1.longname = RX1DAY_LONGNAME;
      request.var1.units = RX1DAY_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_RX1DAY == cdo_operator_id())
    {
      request.var1.name = RX1DAY_NAME_ET;
      request.var1.longname = RX1DAY_LONGNAME_ET;
      request.var1.units = RX1DAY_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
    request.var1.f2 = field2_max;
  }
};

class EcaRx5day : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaRx5day",
    .operators = { { "eca_rx5day", 0, CMP_DATE, Eca_rx5dayHelp },
                   { "etccdi_rx5day", 0, CMP_YEAR, Eca_rx5dayHelp },
                   { "etccdi_rx5daymon", Eca_rx5dayHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaRx5day> registration = RegisterEntry<EcaRx5day>();
  int ECA_RX5DAY, ETCCDI_RX5DAY;
  char longname2[sizeof(RX5DAY_LONGNAME2) + 40];

public:
  void
  init() override
  {
    double argX = 50.0;

    ECA_RX5DAY = module.get_id("eca_rx5day");
    ETCCDI_RX5DAY = module.get_id("etccdi_rx5day");
    set_default_compare_type(request.compare_type);
    if (cdo_operator_argc() > 1)
    {
      argX = parameter_to_double(cdo_operator_argv(0));
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }

    std::snprintf(longname2, sizeof(longname2), RX5DAY_LONGNAME2, argX);

    if (ECA_RX5DAY == cdo_operator_id())
    {
      request.var1.name = RX5DAY_NAME;
      request.var1.longname = RX5DAY_LONGNAME;
      request.var1.units = RX5DAY_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_RX5DAY == cdo_operator_id())
    {
      request.var1.name = RX5DAY_NAME_ET;
      request.var1.longname = RX5DAY_LONGNAME_ET;
      request.var1.units = RX5DAY_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }
    request.var1.f2 = field2_max;
    request.var2.name = RX5DAY_NAME2;
    request.var2.longname = longname2;
    request.var2.units = RX5DAY_UNITS2;
    request.var2.h1 = vfarselgec;
    request.var2.h1arg = argX;
    request.var2.h2 = vfarnum;
  }
};

class EcaSdii : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "EcaSdii",
    .operators = { { "eca_sdii", 0, CMP_DATE, Eca_sdiiHelp }, { "etccdi_sdii", 0, CMP_DATE, Eca_sdiiHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<EcaSdii> registration = RegisterEntry<EcaSdii>();
  int ECA_SDII, ETCCDI_SDII;
  char longname[1024];
  double threshold = 1;

public:
  void
  init() override
  {
    ECA_SDII = module.get_id("eca_sdii");
    ETCCDI_SDII = module.get_id("etccdi_sdii");
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2) cdo_abort("Too many arguments!");

    if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      if (ETCCDI_SDII == cdo_operator_id()) { request.compare_type = CMP_YEAR; }
      set_compare_type_from_params(request.compare_type, params);
    }
    else
    {
      if (cdo_operator_argc() == 1) threshold = parameter_to_double(cdo_operator_argv(0));
    }

    if (ECA_SDII == cdo_operator_id())
    {
      std::snprintf(longname, sizeof(longname), SDII_LONGNAME, threshold);
      request.var1.name = SDII_NAME;
      request.var1.longname = longname;
      request.var1.units = SDII_UNITS;
      request.var1.refdate = ECA_refdate;
    }
    else if (ETCCDI_SDII == cdo_operator_id())
    {
      request.var1.name = SDII_NAME_ET;
      request.var1.longname = SDII_LONGNAME_ET;
      request.var1.units = SDII_UNITS_ET;
      request.var1.refdate = ETC_refdate;
    }

    request.var1.f1 = vfarselgec;
    request.var1.f1arg = threshold;
    request.var1.f2 = field2_sum;
    request.var1.epilog = MEAN;
  }
};

class Fdns : public EcaIndices2
{
public:
  using EcaIndices2::EcaIndices2;
  inline static CdoModule module = {
    .name = "Fdns",
    .operators = { { "fdns", 0, CMP_DATE, FdnsHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 2, 1, NoRestriction },
  };
  inline static RegisterEntry<Fdns> registration = RegisterEntry<Fdns>();

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    request.var1.name = FDNS_NAME;
    request.var1.longname = FDNS_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = FDNS_UNITS;
    request.var1.f1 = vfarsellec;
    request.var1.f1arg = TO_KELVIN(0.0);
    request.var1.f2 = vfarsellec;
    request.var1.f2arg = 0.01;
    request.var1.f3 = field2_add;  // any f with f(a, b) = miss, if a = miss or b = miss will do here
    request.var1.f4 = vfarnum;
  }
};

class Strwin : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "Strwin",
    .operators = { { "strwin", 0, CMP_DATE, StrwinHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Strwin> registration = RegisterEntry<Strwin>();
  double maxWind = 10.5;
  char longname[sizeof(STRWIN_LONGNAME) + 40];

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 2)
      cdo_abort("Too many arguments!");
    else if (cdo_operator_argc() > 1)
    {
      auto params = cdo_get_oper_argv();
      params = std::vector<std::string>(params.begin() + 1, params.end());
      set_compare_type_from_params(request.compare_type, params);
    }
    else
    {
      if (cdo_operator_argc() > 0) maxWind = parameter_to_double(cdo_operator_argv(0));
    }

    std::snprintf(longname, sizeof(longname), STRWIN_LONGNAME, maxWind);

    request.var1.name = STRWIN_NAME;
    request.var1.longname = longname;
    request.var1.refdate = ECA_refdate;
    request.var1.units = STRWIN_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = maxWind;
    request.var1.f2 = vfarnum;
    request.var2.name = STRWIN_NAME2;
    request.var2.longname = STRWIN_LONGNAME2;
    request.var2.units = STRWIN_UNITS2;
    request.var2.h1 = vfarselgec;
    request.var2.h1arg = maxWind;
    request.var2.h2 = vfarnum2;
    request.var2.h3 = field2_max;
  }
};

class Strbre : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "Strbre",
    .operators = { { "strbre", 0, CMP_DATE, StrbreHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Strbre> registration = RegisterEntry<Strbre>();
  double maxWind = 10.5;

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }

    request.var1.name = STRBRE_NAME;
    request.var1.longname = STRBRE_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = STRWIN_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = maxWind;
    request.var1.f2 = vfarnum;
    request.var2.name = STRBRE_NAME2;
    request.var2.longname = STRBRE_LONGNAME2;
    request.var2.units = STRWIN_UNITS2;
    request.var2.h1 = vfarselgec;
    request.var2.h1arg = maxWind;
    request.var2.h2 = vfarnum2;
    request.var2.h3 = field2_max;
  }
};

class Strgal : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "Strgal",
    .operators = { { "strgal", 0, CMP_DATE, StrgalHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Strgal> registration = RegisterEntry<Strgal>();
  double maxWind = 20.5;

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }

    request.var1.name = STRBRE_NAME;
    request.var1.longname = STRBRE_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = STRWIN_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = maxWind;
    request.var1.f2 = vfarnum;
    request.var2.name = STRBRE_NAME2;
    request.var2.longname = STRBRE_LONGNAME2;
    request.var2.units = STRWIN_UNITS2;
    request.var2.h1 = vfarselgec;
    request.var2.h1arg = maxWind;
    request.var2.h2 = vfarnum2;
    request.var2.h3 = field2_max;
  }
};

class Hurr : public EcaIndices1
{
public:
  using EcaIndices1::EcaIndices1;
  inline static CdoModule module = {
    .name = "Hurr",
    .operators = { { "hurr", 0, CMP_DATE, HurrHelp } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Hurr> registration = RegisterEntry<Hurr>();
  double maxWind = 32.5;

public:
  void
  init() override
  {
    set_default_compare_type(request.compare_type);

    if (cdo_operator_argc() > 0)
    {
      auto const &params = cdo_get_oper_argv();
      set_compare_type_from_params(request.compare_type, params);
    }

    request.var1.name = HURR_NAME;
    request.var1.longname = HURR_LONGNAME;
    request.var1.refdate = ECA_refdate;
    request.var1.units = STRWIN_UNITS;
    request.var1.f1 = vfarselgec;
    request.var1.f1arg = maxWind;
    request.var1.f2 = vfarnum;
    request.var2.name = HURR_NAME2;
    request.var2.longname = HURR_LONGNAME2;
    request.var2.units = STRWIN_UNITS2;
    request.var2.h1 = vfarselgec;
    request.var2.h1arg = maxWind;
    request.var2.h2 = vfarnum2;
    request.var2.h3 = field2_max;
  }
};
