# skyborn-cdo

Pre-compiled [CDO (Climate Data Operators)](https://code.mpimet.mpg.de/projects/cdo) distributed as a pip-installable Python package.

This is a backend module for [Skyborn](https://github.com/QianyeSu/Skyborn) — an atmospheric science research toolkit.

### About CDO 2.6.0

This package bundles **CDO 2.6.0** with all required dependencies. For detailed CDO documentation and release information, see:

- **CDO Official Project**: https://code.mpimet.mpg.de/projects/cdo
- **CDO 2.6.0 Release Notes**: https://code.mpimet.mpg.de/projects/cdo/news

## Installation

```bash
pip install skyborn-cdo
```

This installs a pre-compiled CDO binary along with all required libraries (NetCDF, HDF5, ecCodes, FFTW3, PROJ, UDUNITS2). **No system-level CDO installation needed.**

### Optional dependencies

```bash
# xarray support (returnXArray)
pip install skyborn-cdo[xarray]

# Test dependencies
pip install skyborn-cdo[test]
```

### Supported platforms

| Platform | Architecture | Status |
|----------|-------------|--------|
| Linux    | x86_64      | ✅ Supported |
| macOS    | arm64 (Apple Silicon) | ✅ Supported |
| Windows  | x86_64      | ✅ Supported |

## Quick Start

```python
from skyborn_cdo import Cdo

cdo = Cdo()
print(cdo.version())  # CDO 2.6.0
```

## Usage

skyborn-cdo provides two API styles — **command-line style** and **method-call style** — both covering all 900+ CDO operators.

### 1. Command-line Style — `cdo()`

Pass a full CDO command string, just like typing in a terminal. The leading `cdo` prefix is optional.

```python
from skyborn_cdo import Cdo

cdo = Cdo()

# With 'cdo' prefix (copy/paste from terminal)
cdo("cdo -O mergetime in1.nc in2.nc out.nc")

# Without 'cdo' prefix (also valid)
cdo("mergetime in1.nc in2.nc out.nc")

# Complex command with options
cdo("-O -f nc4 sellonlatbox,0,30,0,30 input.nc output.nc")

# Chained operators (CDO pipe syntax)
cdo("-O -f nc4 -fldmean -sellonlatbox,70,140,10,55 input.nc output.nc")
cdo("-O -f nc4 -sp2gpl -setgridtype,regular input.nc output.nc")
```

### 2. Method-call Style — `cdo.operator()`

Each of CDO's 900+ operators is available as a Python method via dynamic dispatch.

```python
from skyborn_cdo import Cdo

cdo = Cdo()

# Basic: operator(parameters, input=..., output=...)
cdo.sellonlatbox("0,30,0,30", input="input.nc", output="output.nc")

# No parameters needed
cdo.copy(input="input.nc", output="output.nc")
cdo.topo(output="topo.nc")

# Multiple input files (space-separated string or list)
cdo.mergetime(input="in1.nc in2.nc in3.nc", output="merged.nc")
cdo.mergetime(input=["in1.nc", "in2.nc", "in3.nc"], output="merged.nc")

# Wildcard / glob patterns (automatically expanded)
cdo.mergetime(input="data_2020*.nc", output="merged.nc")
cdo.mergetime(input="/path/to/data_20200?.nc", output="merged.nc")
cdo.ensmean(input="ensemble_*.nc", output="ensmean.nc")
```

> **Wildcards**: Both `cdo.operator(input="*.nc")` and `cdo("mergetime *.nc out.nc")` support glob patterns (`*`, `?`, `[...]`). Files are sorted alphabetically before being passed to CDO.

### 3. Options

CDO options like `-O` (overwrite), `-s` (silent), `-f nc4` (output format) can be set globally or per-call.

```python
# Global options — applied to every command
cdo = Cdo(options="-O -s")
cdo.copy(input="in.nc", output="out.nc")

# Per-call options — merged with global options
cdo = Cdo()
cdo.copy(input="in.nc", output="out.nc", options="-O -f nc4")

# Common options:
#   -O          Overwrite existing output files
#   -s          Silent mode (suppress informational messages)
#   -f nc4      Output in NetCDF4 format
#   -f nc4c     Output in NetCDF4 classic format
#   -f grb2     Output in GRIB2 format
#   -P 4        Use 4 threads for parallel processing
```

> **⚠️ Important**: CDO **does not overwrite** existing output files by default.
> If the output file already exists, CDO will raise an error:
> `Outputfile out.nc already exists!`
> Use `-O` to enable overwriting: `cdo = Cdo(options="-O")` or pass `options="-O"` per call.

### 4. Info Operators (return text output)

Operators like `info`, `sinfo`, `griddes`, `showname`, etc. return their text output as a string instead of writing to a file.

```python
cdo = Cdo()

# File information
info = cdo.sinfo(input="input.nc")
print(info)

# Grid description
grid = cdo.griddes(input="input.nc")
print(grid)

# Show variable names
names = cdo.showname(input="input.nc")
print(names)

# Show timestamps
dates = cdo.showdate(input="input.nc")
print(dates)

# Number of time steps / levels / variables
print(cdo.ntime(input="input.nc"))
print(cdo.nlevel(input="input.nc"))
print(cdo.nvar(input="input.nc"))
```

### 5. Return as xarray / netCDF4 / numpy

```python
# Return as xarray.Dataset (requires: pip install skyborn-cdo[xarray])
ds = cdo.sellonlatbox("0,30,0,30", input="input.nc", returnXArray=True)
print(ds)

# Return as netCDF4.Dataset (requires: pip install netCDF4)
nc = cdo.copy(input="input.nc", returnCdf=True)
print(nc.variables.keys())

# Return as numpy array (first variable)
arr = cdo.copy(input="input.nc", returnArray=True)
print(arr.shape)

# Return as masked numpy array
marr = cdo.copy(input="input.nc", returnMaArray=True)
```

### 6. Chained Operators (Pipeline)

CDO supports nesting operators in a single command. This is the most powerful feature for building complex processing pipelines.

#### Command-line style (recommended for complex chains)

```python
cdo = Cdo(options="-O")

# Remapping → Spatial selection → Field mean
cdo("-f nc4 -fldmean -sellonlatbox,70,140,10,55 -remapbil,r180x90 input.nc output.nc")

# Spectral mode conversion pipeline
cdo("-f nc4 -sp2gpl -setgridtype,regular input_spectral.nc output_regular.nc")

# Arithmetic pipeline
cdo("-addc,273.15 -mulc,0.01 input.nc output.nc")

# Conditional masking
cdo("-ifthen -gtc,0 topo.nc topo.nc positive_topo.nc")
```

#### Method-call style (using input as sub-command)

```python
cdo = Cdo(options="-O")

# Use CDO operator syntax in the input parameter for chaining
cdo.remapbil("r360x180", input="-mergetime in1.nc in2.nc in3.nc", output="out.nc")
cdo.fldmean(input="-sellonlatbox,70,140,10,55 input.nc", output="mean.nc")
cdo.timmean(input="-sellonlatbox,0,360,-30,30 -remapbil,r180x90 input.nc", output="out.nc")
```

### 7. Common Operator Examples

#### Spatial Operations

```python
cdo = Cdo(options="-O")

# Select region by longitude/latitude box
cdo.sellonlatbox("70,140,10,55", input="global.nc", output="china.nc")

# Select by grid index box
cdo.selindexbox("1,100,1,80", input="input.nc", output="subset.nc")

# Remap to regular grid
cdo.remapbil("r360x180", input="input.nc", output="1deg.nc")     # Bilinear
cdo.remapcon("r360x180", input="input.nc", output="1deg.nc")     # Conservative
cdo.remapnn("r360x180", input="input.nc", output="1deg.nc")      # Nearest neighbor
cdo.remapdis("r360x180", input="input.nc", output="1deg.nc")     # Distance weighted

# Zonal / Meridional mean
cdo.zonmean(input="input.nc", output="zonmean.nc")
cdo.mermean(input="input.nc", output="mermean.nc")

# Invert latitude direction
cdo.invertlat(input="input.nc", output="flipped.nc")
```

#### Time Operations

```python
cdo = Cdo(options="-O")

# Time averaging
cdo.timmean(input="input.nc", output="time_avg.nc")
cdo.timstd(input="input.nc", output="time_std.nc")
cdo.timmin(input="input.nc", output="time_min.nc")
cdo.timmax(input="input.nc", output="time_max.nc")

# Monthly / Seasonal / Yearly statistics
cdo.monmean(input="input.nc", output="monthly_mean.nc")
cdo.seasmean(input="input.nc", output="seasonal_mean.nc")
cdo.yearmonmean(input="input.nc", output="yearly_mean.nc")

# Select time steps
cdo.selyear("2020", input="input.nc", output="year2020.nc")
cdo.selmon("1,2,3", input="input.nc", output="jan_mar.nc")
cdo.seltimestep("1/10", input="input.nc", output="first10.nc")

# Set time axis
cdo.settaxis("2020-01-15,12:00,1mon", input="input.nc", output="redate.nc")

# Merge time series
cdo.mergetime(input="jan.nc feb.nc mar.nc", output="q1.nc")

# Split by month / year
cdo.splitmon(input="input.nc", output="monthly_")
cdo.splityear(input="input.nc", output="yearly_")

# Detrend
cdo.detrend(input="input.nc", output="detrended.nc")
```

#### Field Statistics

```python
cdo = Cdo(options="-O")

# Field (spatial) statistics
cdo.fldmean(input="input.nc", output="fldmean.nc")
cdo.fldstd(input="input.nc", output="fldstd.nc")
cdo.fldmin(input="input.nc", output="fldmin.nc")
cdo.fldmax(input="input.nc", output="fldmax.nc")
cdo.fldsum(input="input.nc", output="fldsum.nc")
```

#### Arithmetic

```python
cdo = Cdo(options="-O")

# Scalar operations
cdo.mulc("2.0", input="input.nc", output="doubled.nc")
cdo.addc("273.15", input="celsius.nc", output="kelvin.nc")
cdo.divc("100", input="input.nc", output="divided.nc")
cdo.subc("273.15", input="kelvin.nc", output="celsius.nc")

# Unary operations
cdo.abs(input="input.nc", output="absolute.nc")
cdo.sqrt(input="positive.nc", output="sqrt.nc")

# Two-file operations
cdo.add(input="a.nc b.nc", output="sum.nc")
cdo.sub(input="a.nc b.nc", output="diff.nc")
cdo.mul(input="a.nc b.nc", output="product.nc")
cdo.div(input="a.nc b.nc", output="ratio.nc")

# Custom expression
cdo.expr("'new_var=temp*2+precip;'", input="input.nc", output="computed.nc")
```

#### Vertical Level Operations

```python
cdo = Cdo(options="-O")

# Select specific levels
cdo.sellevel("85000,50000,20000", input="input.nc", output="levels.nc")

# Interpolate to new levels
cdo.intlevel("90000,85000,70000,50000,30000", input="input.nc", output="interp.nc")

# Show levels
print(cdo.showlevel(input="input.nc"))
```

#### Grid Type Conversion & Spectral Transforms

```python
cdo = Cdo(options="-O")

# Convert grid type
cdo.setgridtype("regular", input="input.nc", output="regular.nc")

# Spectral ↔ Grid-point transforms
cdo.sp2gpl(input="spectral.nc", output="gaussian_linear.nc")
cdo.sp2gp(input="spectral.nc", output="gaussian.nc")
cdo.gp2sp(input="gaussian.nc", output="spectral.nc")

# Complex chain: spectral to regular grid in NetCDF4
cdo("-O -f nc4 -sp2gpl -setgridtype,regular spectral.nc regular.nc")
```

#### Format Conversion

```python
cdo = Cdo(options="-O")

# NetCDF formats
cdo.copy(input="input.nc", output="out.nc", options="-f nc4")     # NetCDF4
cdo.copy(input="input.nc", output="out.nc", options="-f nc4c")    # NetCDF4 Classic
cdo.copy(input="input.nc", output="out.nc", options="-f nc2")     # NetCDF 64-bit

# GRIB formats
cdo.copy(input="input.nc", output="out.grb", options="-f grb")    # GRIB1
cdo.copy(input="input.nc", output="out.grb2", options="-f grb2")  # GRIB2
```

#### Variable Metadata

```python
cdo = Cdo(options="-O")

# Rename variable
cdo.chname("old_name,new_name", input="input.nc", output="renamed.nc")

# Set attributes
cdo.setattribute("varname@units=kg/m2", input="input.nc", output="units.nc")

# Set missing value
cdo.setmissval("-999", input="input.nc", output="newmiss.nc")
```

#### Ensemble Operations

```python
cdo = Cdo(options="-O")

# Ensemble mean (multiple realizations)
cdo.ensmean(input="run_*.nc", output="ensemble_mean.nc")

# Ensemble standard deviation
cdo.ensstd(input="run_*.nc", output="ensemble_std.nc")

# Ensemble variance
cdo.ensvar(input="run_001.nc run_002.nc run_003.nc", output="ensemble_var.nc")

# Ensemble sum
cdo.enssum(input=["member_1.nc", "member_2.nc", "member_3.nc"], output="ens_sum.nc")

# Ensemble percentiles
cdo.enspctl("25,50,75", input="run_*.nc", output="ensemble_quartiles.nc")
```

#### Grid Information and Modification

```python
cdo = Cdo(options="-O")

# Calculate grid cell areas
cdo.gridarea(input="input.nc", output="grid_areas.nc")

# Calculate grid weights (for weighted averaging)
cdo.gridweights(input="input.nc", output="weights.nc")

# Set grid type
cdo.setgridtype("regular", input="curvilinear.nc", output="regular.nc")

# Show grid description
print(cdo.griddes(input="input.nc"))

# Invert latitude direction (flip N-S)
cdo.invertlat(input="input.nc", output="flipped.nc")
```

#### Advanced Statistical Operations

```python
cdo = Cdo(options="-O")

# Field percentiles
cdo.fldpctl("10,50,90", input="input.nc", output="field_percentiles.nc")

# Field range (max - min)
cdo.fldrange(input="input.nc", output="field_range.nc")

# Time series percentiles
cdo.timpctl("25,50,75", input="timeseries.nc", output="time_percentiles.nc")

# Running mean (e.g., 5-timestep window)
cdo.runmean("5", input="input.nc", output="smoothed.nc")

# Trend removal
cdo.detrend(input="input.nc", output="detrended.nc")
```

### 8. Error Handling

skyborn-cdo captures all CDO errors and raises them as `CdoError` exceptions with full diagnostic information.

```python
from skyborn_cdo import Cdo, CdoError

cdo = Cdo()

try:
    cdo.sellonlatbox("0,30,0,30", input="nonexistent.nc", output="out.nc")
except CdoError as e:
    print(f"Error: {e}")
    print(f"Return code: {e.returncode}")
    print(f"CDO stderr: {e.stderr}")
    print(f"Command: {e.cmd}")
```

Common errors:
- **File not found**: `Open failed on >file.nc< No such file or directory`
- **Output exists**: `Outputfile out.nc already exists!` → Add `-O` option
- **Invalid parameters**: `Float parameter >abc< contains invalid character`
- **Timeout**: `CDO command timed out after 60s` → Increase timeout

### 9. Timeout

```python
# Global timeout for all commands (seconds)
cdo = Cdo(timeout=300)

# Per-call timeout override
cdo.remapcon("r3600x1800", input="input.nc", output="hires.nc", timeout=600)
cdo("cdo -O remapcon,r3600x1800 input.nc hires.nc", timeout=600)
```

### 10. Debug Mode

```python
# Print executed commands and CDO stderr output
cdo = Cdo(debug=True)
cdo.sellonlatbox("0,30,0,30", input="input.nc", output="output.nc")
# [skyborn-cdo] Running: /path/to/cdo -sellonlatbox,0,30,0,30 input.nc output.nc
# [skyborn-cdo] stderr: cdo    sellonlatbox: ...
```

### 11. Getting Help

#### Python API — `cdo.help()`

```python
from skyborn_cdo import Cdo

cdo = Cdo()

# Help for a specific operator
print(cdo.help("sellonlatbox"))
print(cdo.help("mergetime"))
print(cdo.help("remapbil"))

# General usage summary
print(cdo.help())

# List all available operators
print(cdo.operators())

# Check if an operator exists
cdo.has_operator("sellonlatbox")  # True
```

#### CLI — `skyborn-cdo -h` / `--help`

```bash
# Show general help
skyborn-cdo --help

# Show help for a specific CDO operator (passed to CDO directly)
skyborn-cdo -h sellonlatbox
skyborn-cdo -h mergetime
skyborn-cdo -h remapbil

# List all operators
skyborn-cdo --operators

# Also works via python -m
python -m skyborn_cdo --help
python -m skyborn_cdo -h sellonlatbox
```

## CLI

The package provides a `skyborn-cdo` command-line tool (also available as `python -m skyborn_cdo`):

```bash
# Show installation info and CDO version
skyborn-cdo --info

# Show help
skyborn-cdo --help

# Operator-level help (forwarded to CDO)
skyborn-cdo -h sellonlatbox

# List all operators
skyborn-cdo --operators

# Pass-through to CDO (any CDO command works)
skyborn-cdo -O -f nc4 copy input.nc output.nc
skyborn-cdo mergetime in1.nc in2.nc out.nc
skyborn-cdo -O -f nc4 -sp2gpl -setgridtype,regular spectral.nc regular.nc
```

## CDO Operators

skyborn-cdo provides access to **938 CDO operators** covering all aspects of climate data processing. Check operator availability:

```python
cdo = Cdo()
print(len(cdo.operators()))  # 938
cdo.has_operator("sellonlatbox")  # True
print(cdo.help("sellonlatbox"))  # Get operator documentation
```

### Operator Categories

| Category | Count | Examples | Use Cases |
|----------|-------|----------|-----------|
| **Time Operations** | 208 | `timmean`, `timstd`, `mergetime`, `seldate`, `monmean`, `yearsum` | Time series analysis, temporal statistics, date/time selection |
| **Statistical** | 222 | `fldmean`, `fldstd`, `ensmean`, `timavg`, `zonmean`, `mermean` | Spatial/temporal statistics, ensemble analysis |
| **Spatial Selection** | 35 | `sellonlatbox`, `selindexbox`, `selgrid`, `sellevel`, `selcode` | Regional extraction, level selection |
| **Grid/Remapping** | 64 | `remapbil`, `remapcon`, `remapnn`, `setgrid`, `gridarea` | Grid conversion, interpolation, regridding |
| **Arithmetic** | 60 | `add`, `mul`, `expr`, `addc`, `sqrt`, `log` | Mathematical operations, calculations |
| **Vertical Levels** | 13 | `ml2pl`, `intlevel`, `pressure`, `sealevelpressure` | Model level conversion, vertical interpolation |
| **Format/IO** | 47 | `copy`, `merge`, `split`, `cat`, `import_*`, `output*` | File operations, format conversion |
| **Info/Query** | 50 | `info`, `showname`, `griddes`, `ntime`, `showdate` | File inspection, metadata extraction |
| **Spectral/Grid Transform** | 28 | `sp2gp`, `gp2sp`, `sp2gpl`, `fourier` | Spectral transforms, Fourier analysis |
| **Others** | 211 | Specialized operators for specific domains | Advanced climate computations |

**Total: 938 operators**

### Are All Operators Useful?

**Commonly used** (daily work): ~80-100 operators
- Space: `sellonlatbox`, `remapbil`, `remapcon`, `zonmean`, `mermean`
- Time: `mergetime`, `timmean`, `monmean`, `yearsum`, `seldate`, `selyear`
- Statistics: `fldmean`, `fldstd`, `fldmin`, `fldmax`, `ensmean`
- Arithmetic: `add`, `sub`, `mul`, `div`, `addc`, `mulc`, `expr`
- Info: `sinfo`, `showname`, `griddes`, `ntime`

**Specialized operators** (400+): Domain-specific
- Meteorology: `sealevelpressure`, `ml2pl`, `geopotheight`
- Oceanography: `detrend`, `dmean`, `seasmean`
- Climate indices: `eca_*` (extreme climate events), `ydrun*` (running means)
- Statistical analysis: `trend`, `regres`, `corr`, `eof`
- Spectral analysis: `sp2gp`, `gp2sp`, `dft`, `filter`

**Legacy/Niche operators** (400+): Less frequently used but available
- Format-specific imports (`import_cmsaf`, `import_grads`)
- Experimental features (`remapavgtest`, `remapcon2test`)
- Highly specialized computations

**How to decide if you need an operator:**
```python
# Search operators by keyword
cdo = Cdo()
time_ops = [op for op in cdo.operators() if 'time' in op]
print(f"Time-related: {len(time_ops)} operators")  # 208

# Get detailed help for any operator
print(cdo.help("operator_name"))
```

## CDO Version

This package bundles **CDO 2.6.0** with the following libraries:
- NetCDF-C 4.9.x
- HDF5 1.14.x
- ecCodes 2.40+
- FFTW3 3.3.10
- PROJ 9.5.x
- UDUNITS2 2.2.28

Exact library versions vary by platform (Linux/macOS build from source, Windows uses MSYS2 packages).

### What's New in CDO 2.6.0

#### New operators

| Operator | Description |
|----------|-------------|
| `varsskew` | Ensemble skewness across input files (VarsStat group) |
| `varskurt` | Ensemble kurtosis across input files (VarsStat group) |
| `varsmedian` | Ensemble median across input files (VarsStat group) |
| `varspctl` | Ensemble percentile across input files, e.g. `varspctl,90` |
| `symmetrize` | Mirrors data at the equator (creates symmetric fields) |
| `splitensemble` | Splits GRIB2 ensemble members into separate files |

#### New global options

| Option | Description |
|--------|-------------|
| `--query` | Pre-selects a subset of the data cube from a dataset |
| `--async_read true\|false` | Reads input data asynchronously; available for `diff`, `info`, `trend`, `detrend`, `Timstat` operators |

#### Performance improvements

- Significant performance improvement for reading HEALPix zarr datasets with NCZARR.

#### Bug fixes

- `fillmiss`: wrong result when fewer than 4 valid neighbours available [Bug #12341]
- `chparam`: failed since release 2.5.3 [Bug #12328]

#### Usage examples

```python
cdo = Cdo()

# VarsStat: compute skewness / kurtosis / median across an ensemble
cdo.varsskew(input=["mem1.nc", "mem2.nc", "mem3.nc"], output="skew.nc")
cdo.varskurt(input=["mem1.nc", "mem2.nc", "mem3.nc"], output="kurt.nc")
cdo.varsmedian(input=["mem1.nc", "mem2.nc", "mem3.nc"], output="median.nc")
cdo.varspctl("90", input=["mem1.nc", "mem2.nc", "mem3.nc"], output="pct90.nc")

# Symmetrize: mirror a field at the equator
cdo.symmetrize(input="topo.nc", output="topo_symmetric.nc")

# Global async read option (speeds up supported operators)
cdo("--async_read true -diff file1.nc file2.nc")
```

## Operator Reference

Complete listing of all operators available in CDO 2.6.0 (bundled with skyborn-cdo), grouped by category. Use `cdo.help("operator_name")` for detailed usage.

### Information

| Operator | Description |
|----------|-------------|
| `info` | Dataset information listed by identifier |
| `infon` | Dataset information listed by name |
| `cinfo` | Compact information listed by name |
| `map` | Dataset information and simple map |
| `sinfo` | Short information listed by identifier |
| `sinfon` | Short information listed by name |
| `xsinfo` | Extra short information listed by name |
| `xsinfop` | Extra short information listed by identifier |
| `diff` | Compare two datasets listed by identifier |
| `diffn` | Compare two datasets listed by name |
| `npar` | Number of parameters |
| `nlevel` | Number of levels |
| `nyear` | Number of years |
| `nmon` | Number of months |
| `ndate` | Number of dates |
| `ntime` | Number of timesteps |
| `ngridpoints` | Number of gridpoints |
| `ngrids` | Number of horizontal grids |
| `showformat` | Show file format |
| `showcode` | Show code numbers |
| `showname` | Show variable names |
| `showstdname` | Show standard names |
| `showlevel` | Show levels |
| `showltype` | Show GRIB level types |
| `showyear` | Show years |
| `showmon` | Show months |
| `showdate` | Show date information |
| `showtime` | Show time information |
| `showtimestamp` | Show timestamp |
| `showchunkspec` | Show chunk specification |
| `showfilter` | Show filter specification |
| `showattribute` | Show attributes |
| `partab` | Parameter table |
| `codetab` | Parameter code table |
| `griddes` | Grid description |
| `zaxisdes` | Z-axis description |
| `vct` | Vertical coordinate table |

### File Operations

| Operator | Description |
|----------|-------------|
| `copy` | Copy datasets |
| `clone` | Clone datasets |
| `cat` | Concatenate datasets |
| `tee` | Duplicate a data stream and write it to file |
| `pack` | Pack data |
| `unpack` | Unpack data |
| `setchunkspec` | Specify chunking |
| `setfilter` | Specify filter |
| `bitrounding` | Bit rounding |
| `replace` | Replace variables |
| `duplicate` | Duplicates a dataset |
| `mergegrid` | Merge grid |
| `merge` | Merge datasets with different fields |
| `mergetime` | Merge datasets sorted by date and time |
| `splitcode` | Split code numbers |
| `splitparam` | Split parameter identifiers |
| `splitname` | Split variable names |
| `splitlevel` | Split levels |
| `splitgrid` | Split grids |
| `splitzaxis` | Split z-axes |
| `splittabnum` | Split parameter table numbers |
| `splithour` | Split hours |
| `splitday` | Split days |
| `splitseas` | Split seasons |
| `splityear` | Split years |
| `splityearmon` | Split in years and months |
| `splitmon` | Split months |
| `splitsel` | Split selected timesteps |
| `splitdate` | Splits a file into dates |
| `distgrid` | Distribute horizontal grid |
| `collgrid` | Collect horizontal grid |

### Selection

| Operator | Description |
|----------|-------------|
| `select` | Select fields |
| `delete` | Delete fields |
| `selmulti` | Select multiple fields |
| `delmulti` | Delete multiple fields |
| `changemulti` | Change identification of multiple fields |
| `selparam` | Select parameters by identifier |
| `delparam` | Delete parameters by identifier |
| `selcode` | Select parameters by code number |
| `delcode` | Delete parameters by code number |
| `selname` | Select parameters by name |
| `delname` | Delete parameters by name |
| `selstdname` | Select parameters by standard name |
| `sellevel` | Select levels |
| `sellevidx` | Select levels by index |
| `selgrid` | Select grids |
| `selzaxis` | Select z-axes |
| `selzaxisname` | Select z-axes by name |
| `selltype` | Select GRIB level types |
| `seltabnum` | Select parameter table numbers |
| `seltimestep` | Select timesteps |
| `seltime` | Select times |
| `selhour` | Select hours |
| `selday` | Select days |
| `selmonth` | Select months (alias: `selmon`) |
| `selyear` | Select years |
| `selseason` | Select seasons |
| `seldate` | Select dates |
| `selsmon` | Select single month |
| `sellonlatbox` | Select a longitude/latitude box |
| `selindexbox` | Select an index box |
| `selregion` | Select cells inside regions |
| `selcircle` | Select cells inside a circle |
| `selgridcell` | Select grid cells |
| `delgridcell` | Delete grid cells |
| `samplegrid` | Resample grid cells |
| `selyearidx` | Select year by index |
| `seltimeidx` | Select timestep by index |
| `bottomvalue` | Extract bottom level |
| `topvalue` | Extract top level |
| `isosurface` | Extract isosurface |

### Conditional

| Operator | Description |
|----------|-------------|
| `ifthen` | If then |
| `ifnotthen` | If not then |
| `ifthenelse` | Conditional selection |
| `ifthenc` | If then constant |
| `ifnotthenc` | If not then constant |
| `reducegrid` | Reduce fields to user-defined mask |

### Comparison

| Operator | Description |
|----------|-------------|
| `eq` | Equal |
| `ne` | Not equal |
| `le` | Less equal |
| `lt` | Less than |
| `ge` | Greater equal |
| `gt` | Greater than |
| `eqc` | Equal constant |
| `nec` | Not equal constant |
| `lec` | Less equal constant |
| `ltc` | Less than constant |
| `gec` | Greater equal constant |
| `gtc` | Greater than constant |
| `ymoneq` | Compare time series with Equal (multi-year monthly) |
| `ymonne` | Compare time series with NotEqual (multi-year monthly) |
| `ymonle` | Compare time series with LessEqual (multi-year monthly) |
| `ymonlt` | Compare time series with LessThan (multi-year monthly) |
| `ymonge` | Compare time series with GreaterEqual (multi-year monthly) |
| `ymongt` | Compare time series with GreaterThan (multi-year monthly) |
| `yseaseq` | Compare time series with Equal (multi-year seasonal) |
| `yseasne` | Compare time series with NotEqual (multi-year seasonal) |
| `yseasle` | Compare time series with LessEqual (multi-year seasonal) |
| `yseaslt` | Compare time series with LessThan (multi-year seasonal) |
| `yseasge` | Compare time series with GreaterEqual (multi-year seasonal) |
| `yseasgt` | Compare time series with GreaterThan (multi-year seasonal) |

### Modification

| Operator | Description |
|----------|-------------|
| `setattribute` | Set attributes |
| `delattribute` | Delete attributes |
| `setpartabp` | Set parameter table (by path) |
| `setpartabn` | Set parameter table (by name) |
| `setcodetab` | Set parameter code table |
| `setcode` | Set code number |
| `setparam` | Set parameter identifier |
| `setname` | Set variable name |
| `setstdname` | Set standard name |
| `setunit` | Set variable unit |
| `setlevel` | Set level |
| `setltype` | Set GRIB level type |
| `setmaxsteps` | Set max timesteps |
| `setdate` | Set date |
| `settime` | Set time of the day |
| `setday` | Set day |
| `setmon` | Set month |
| `setyear` | Set year |
| `settunits` | Set time units |
| `settaxis` | Set time axis |
| `settbounds` | Set time bounds |
| `setreftime` | Set reference time |
| `setcalendar` | Set calendar |
| `shifttime` | Shift timesteps |
| `chcode` | Change code number |
| `chparam` | Change parameter identifier |
| `chname` | Change variable or coordinate name |
| `chunit` | Change variable unit |
| `chlevel` | Change level |
| `chlevelc` | Change level of one code |
| `chlevelv` | Change level of one variable |
| `setgrid` | Set grid |
| `setgridtype` | Set grid type |
| `setgridarea` | Set grid cell area |
| `setgridmask` | Set grid mask |
| `setprojparams` | Set proj params |
| `setzaxis` | Set z-axis |
| `genlevelbounds` | Generate level bounds |
| `invertlat` | Invert latitudes |
| `invertlev` | Invert levels |
| `shiftx` | Shift x |
| `shifty` | Shift y |
| `maskregion` | Mask regions |
| `masklonlatbox` | Mask a longitude/latitude box |
| `maskindexbox` | Mask an index box |
| `setclonlatbox` | Set a longitude/latitude box to constant |
| `setcindexbox` | Set an index box to constant |
| `enlarge` | Enlarge fields |
| `setmissval` | Set a new missing value |
| `setctomiss` | Set constant to missing value |
| `setmisstoc` | Set missing value to constant |
| `setrtomiss` | Set range to missing value |
| `setvrange` | Set valid range |
| `setmisstonn` | Set missing value to nearest neighbor |
| `setmisstodis` | Set missing value to distance-weighted average |
| `vertfillmiss` | Vertical filling of missing values |
| `timfillmiss` | Temporal filling of missing values |
| `setgridcell` | Set the value of a grid cell |

### Arithmetic

| Operator | Description |
|----------|-------------|
| `expr` | Evaluate expressions |
| `exprf` | Evaluate expressions script |
| `aexpr` | Evaluate expressions and append results |
| `aexprf` | Evaluate expression script and append results |
| `abs` | Absolute value |
| `int` | Integer value |
| `nint` | Nearest integer value |
| `pow` | Power |
| `sqr` | Square |
| `sqrt` | Square root |
| `exp` | Exponential |
| `ln` | Natural logarithm |
| `log10` | Base 10 logarithm |
| `sin` | Sine |
| `cos` | Cosine |
| `tan` | Tangent |
| `asin` | Arc sine |
| `acos` | Arc cosine |
| `atan` | Arc tangent |
| `reci` | Reciprocal value |
| `not` | Logical NOT |
| `addc` | Add a constant |
| `subc` | Subtract a constant |
| `mulc` | Multiply with a constant |
| `divc` | Divide by a constant |
| `minc` | Minimum of a field and a constant |
| `maxc` | Maximum of a field and a constant |
| `add` | Add two fields |
| `sub` | Subtract two fields |
| `mul` | Multiply two fields |
| `div` | Divide two fields |
| `min` | Minimum of two fields |
| `max` | Maximum of two fields |
| `atan2` | Arc tangent of two fields |
| `setmiss` | Set missing values |
| `dayadd` | Add daily time series |
| `daysub` | Subtract daily time series |
| `daymul` | Multiply daily time series |
| `daydiv` | Divide daily time series |
| `monadd` | Add monthly time series |
| `monsub` | Subtract monthly time series |
| `monmul` | Multiply monthly time series |
| `mondiv` | Divide monthly time series |
| `yearadd` | Add yearly time series |
| `yearsub` | Subtract yearly time series |
| `yearmul` | Multiply yearly time series |
| `yeardiv` | Divide yearly time series |
| `yhouradd` | Add multi-year hourly time series |
| `yhoursub` | Subtract multi-year hourly time series |
| `yhourmul` | Multiply multi-year hourly time series |
| `yhourdiv` | Divide multi-year hourly time series |
| `ydayadd` | Add multi-year daily time series |
| `ydaysub` | Subtract multi-year daily time series |
| `ydaymul` | Multiply multi-year daily time series |
| `ydaydiv` | Divide multi-year daily time series |
| `ymonadd` | Add multi-year monthly time series |
| `ymonsub` | Subtract multi-year monthly time series |
| `ymonmul` | Multiply multi-year monthly time series |
| `ymondiv` | Divide multi-year monthly time series |
| `yseasadd` | Add multi-year seasonal time series |
| `yseassub` | Subtract multi-year seasonal time series |
| `yseasmul` | Multiply multi-year seasonal time series |
| `yseasdiv` | Divide multi-year seasonal time series |
| `muldpm` | Multiply with days per month |
| `divdpm` | Divide by days per month |
| `muldpy` | Multiply with days per year |
| `divdpy` | Divide by days per year |
| `mulcoslat` | Multiply with the cosine of the latitude |
| `divcoslat` | Divide by cosine of the latitude |

### Statistics

#### Cumulative / Consecutive

| Operator | Description |
|----------|-------------|
| `timcumsum` | Cumulative sum over all timesteps |
| `consecsum` | Consecutive sum |
| `consects` | Consecutive timesteps |

#### Variables Statistics (across multiple input files)

| Operator | Description |
|----------|-------------|
| `varsmin` | Variables minimum |
| `varsmax` | Variables maximum |
| `varsrange` | Variables range |
| `varssum` | Variables sum |
| `varsmean` | Variables mean |
| `varsavg` | Variables average |
| `varsstd` | Variables standard deviation |
| `varsstd1` | Variables standard deviation (n-1) |
| `varsvar` | Variables variance |
| `varsvar1` | Variables variance (n-1) |
| `varsskew` | Variables skewness |
| `varskurt` | Variables kurtosis |
| `varsmedian` | Variables median |
| `varspctl` | Variables percentile |

#### Ensemble Statistics

| Operator | Description |
|----------|-------------|
| `ensmin` | Ensemble minimum |
| `ensmax` | Ensemble maximum |
| `ensrange` | Ensemble range |
| `enssum` | Ensemble sum |
| `ensmean` | Ensemble mean |
| `ensavg` | Ensemble average |
| `ensstd` | Ensemble standard deviation |
| `ensstd1` | Ensemble standard deviation (n-1) |
| `ensvar` | Ensemble variance |
| `ensvar1` | Ensemble variance (n-1) |
| `ensskew` | Ensemble skewness |
| `enskurt` | Ensemble kurtosis |
| `ensmedian` | Ensemble median |
| `enspctl` | Ensemble percentile |
| `ensrkhistspace` | Ranked Histogram averaged over space |
| `ensrkhisttime` | Ranked Histogram averaged over time |
| `ensroc` | Ensemble Receiver Operating characteristics |
| `enscrps` | Ensemble CRPS and decomposition |
| `ensbrs` | Ensemble Brier score |

#### Field Statistics (spatial reduction to scalar per timestep)

| Operator | Description |
|----------|-------------|
| `fldmin` | Field minimum |
| `fldmax` | Field maximum |
| `fldrange` | Field range |
| `fldsum` | Field sum |
| `fldint` | Field integral |
| `fldmean` | Field mean |
| `fldavg` | Field average |
| `fldstd` | Field standard deviation |
| `fldstd1` | Field standard deviation (n-1) |
| `fldvar` | Field variance |
| `fldvar1` | Field variance (n-1) |
| `fldskew` | Field skewness |
| `fldkurt` | Field kurtosis |
| `fldmedian` | Field median |
| `fldcount` | Field count |
| `fldpctl` | Field percentile |

#### Zonal Statistics (average along longitude axis)

| Operator | Description |
|----------|-------------|
| `zonmin` | Zonal minimum |
| `zonmax` | Zonal maximum |
| `zonrange` | Zonal range |
| `zonsum` | Zonal sum |
| `zonmean` | Zonal mean |
| `zonavg` | Zonal average |
| `zonstd` | Zonal standard deviation |
| `zonstd1` | Zonal standard deviation (n-1) |
| `zonvar` | Zonal variance |
| `zonvar1` | Zonal variance (n-1) |
| `zonskew` | Zonal skewness |
| `zonkurt` | Zonal kurtosis |
| `zonmedian` | Zonal median |
| `zonpctl` | Zonal percentile |

#### Meridional Statistics (average along latitude axis)

| Operator | Description |
|----------|-------------|
| `mermin` | Meridional minimum |
| `mermax` | Meridional maximum |
| `merrange` | Meridional range |
| `mersum` | Meridional sum |
| `mermean` | Meridional mean |
| `meravg` | Meridional average |
| `merstd` | Meridional standard deviation |
| `merstd1` | Meridional standard deviation (n-1) |
| `mervar` | Meridional variance |
| `mervar1` | Meridional variance (n-1) |
| `merskew` | Meridional skewness |
| `merkurt` | Meridional kurtosis |
| `mermedian` | Meridional median |
| `merpctl` | Meridional percentile |

#### Gridbox Statistics (aggregate rectangular blocks of grid cells)

| Operator | Description |
|----------|-------------|
| `gridboxmin` | Gridbox minimum |
| `gridboxmax` | Gridbox maximum |
| `gridboxrange` | Gridbox range |
| `gridboxsum` | Gridbox sum |
| `gridboxmean` | Gridbox mean |
| `gridboxavg` | Gridbox average |
| `gridboxstd` | Gridbox standard deviation |
| `gridboxstd1` | Gridbox standard deviation (n-1) |
| `gridboxvar` | Gridbox variance |
| `gridboxvar1` | Gridbox variance (n-1) |
| `gridboxskew` | Gridbox skewness |
| `gridboxkurt` | Gridbox kurtosis |
| `gridboxmedian` | Gridbox median |

#### Remap Statistics

| Operator | Description |
|----------|-------------|
| `remapmin` | Remap minimum |
| `remapmax` | Remap maximum |
| `remaprange` | Remap range |
| `remapsum` | Remap sum |
| `remapmean` | Remap mean |
| `remapavg` | Remap average |
| `remapstd` | Remap standard deviation |
| `remapstd1` | Remap standard deviation (n-1) |
| `remapvar` | Remap variance |
| `remapvar1` | Remap variance (n-1) |
| `remapskew` | Remap skewness |
| `remapkurt` | Remap kurtosis |
| `remapmedian` | Remap median |

#### Vertical Statistics

| Operator | Description |
|----------|-------------|
| `vertmin` | Vertical minimum |
| `vertmax` | Vertical maximum |
| `vertrange` | Vertical range |
| `vertsum` | Vertical sum |
| `vertmean` | Vertical mean |
| `vertavg` | Vertical average |
| `vertstd` | Vertical standard deviation |
| `vertstd1` | Vertical standard deviation (n-1) |
| `vertvar` | Vertical variance |
| `vertvar1` | Vertical variance (n-1) |

#### Time Selection Statistics

| Operator | Description |
|----------|-------------|
| `timselmin` | Time selection minimum |
| `timselmax` | Time selection maximum |
| `timselrange` | Time selection range |
| `timselsum` | Time selection sum |
| `timselmean` | Time selection mean |
| `timselavg` | Time selection average |
| `timselstd` | Time selection standard deviation |
| `timselstd1` | Time selection standard deviation (n-1) |
| `timselvar` | Time selection variance |
| `timselvar1` | Time selection variance (n-1) |
| `timselpctl` | Time range percentile |

#### Running Statistics

| Operator | Description |
|----------|-------------|
| `runmin` | Running minimum |
| `runmax` | Running maximum |
| `runrange` | Running range |
| `runsum` | Running sum |
| `runmean` | Running mean |
| `runavg` | Running average |
| `runstd` | Running standard deviation |
| `runstd1` | Running standard deviation (n-1) |
| `runvar` | Running variance |
| `runvar1` | Running variance (n-1) |
| `runpctl` | Running percentile |

#### Time Statistics (reduce all timesteps to one)

| Operator | Description |
|----------|-------------|
| `timmin` | Time minimum |
| `timmax` | Time maximum |
| `timminidx` | Index of time minimum |
| `timmaxidx` | Index of time maximum |
| `timrange` | Time range |
| `timsum` | Time sum |
| `timmean` | Time mean |
| `timavg` | Time average |
| `timstd` | Time standard deviation |
| `timstd1` | Time standard deviation (n-1) |
| `timvar` | Time variance |
| `timvar1` | Time variance (n-1) |
| `timpctl` | Temporal percentile |

#### Hourly Statistics

| Operator | Description |
|----------|-------------|
| `hourmin` | Hourly minimum |
| `hourmax` | Hourly maximum |
| `hourrange` | Hourly range |
| `hoursum` | Hourly sum |
| `hourmean` | Hourly mean |
| `houravg` | Hourly average |
| `hourstd` | Hourly standard deviation |
| `hourstd1` | Hourly standard deviation (n-1) |
| `hourvar` | Hourly variance |
| `hourvar1` | Hourly variance (n-1) |
| `hourpctl` | Hourly percentile |

#### Daily Statistics

| Operator | Description |
|----------|-------------|
| `daymin` | Daily minimum |
| `daymax` | Daily maximum |
| `dayrange` | Daily range |
| `daysum` | Daily sum |
| `daymean` | Daily mean |
| `dayavg` | Daily average |
| `daystd` | Daily standard deviation |
| `daystd1` | Daily standard deviation (n-1) |
| `dayvar` | Daily variance |
| `dayvar1` | Daily variance (n-1) |
| `daypctl` | Daily percentile |

#### Monthly Statistics

| Operator | Description |
|----------|-------------|
| `monmin` | Monthly minimum |
| `monmax` | Monthly maximum |
| `monrange` | Monthly range |
| `monsum` | Monthly sum |
| `monmean` | Monthly mean |
| `monavg` | Monthly average |
| `monstd` | Monthly standard deviation |
| `monstd1` | Monthly standard deviation (n-1) |
| `monvar` | Monthly variance |
| `monvar1` | Monthly variance (n-1) |
| `monpctl` | Monthly percentile |

#### Yearly Statistics

| Operator | Description |
|----------|-------------|
| `yearmonmean` | Yearly mean from monthly data |
| `yearmin` | Yearly minimum |
| `yearmax` | Yearly maximum |
| `yearminidx` | Index of yearly minimum |
| `yearmaxidx` | Index of yearly maximum |
| `yearrange` | Yearly range |
| `yearsum` | Yearly sum |
| `yearmean` | Yearly mean |
| `yearavg` | Yearly average |
| `yearstd` | Yearly standard deviation |
| `yearstd1` | Yearly standard deviation (n-1) |
| `yearvar` | Yearly variance |
| `yearvar1` | Yearly variance (n-1) |
| `yearpctl` | Yearly percentile |

#### Seasonal Statistics

| Operator | Description |
|----------|-------------|
| `seasmin` | Seasonal minimum |
| `seasmax` | Seasonal maximum |
| `seasrange` | Seasonal range |
| `seassum` | Seasonal sum |
| `seasmean` | Seasonal mean |
| `seasavg` | Seasonal average |
| `seasstd` | Seasonal standard deviation |
| `seasstd1` | Seasonal standard deviation (n-1) |
| `seasvar` | Seasonal variance |
| `seasvar1` | Seasonal variance (n-1) |
| `seaspctl` | Seasonal percentile |

#### Multi-year Hourly Statistics

| Operator | Description |
|----------|-------------|
| `yhourmin` | Multi-year hourly minimum |
| `yhourmax` | Multi-year hourly maximum |
| `yhourrange` | Multi-year hourly range |
| `yhoursum` | Multi-year hourly sum |
| `yhourmean` | Multi-year hourly mean |
| `yhouravg` | Multi-year hourly average |
| `yhourstd` | Multi-year hourly standard deviation |
| `yhourstd1` | Multi-year hourly standard deviation (n-1) |
| `yhourvar` | Multi-year hourly variance |
| `yhourvar1` | Multi-year hourly variance (n-1) |

#### Multi-day Hourly Statistics

| Operator | Description |
|----------|-------------|
| `dhourmin` | Multi-day hourly minimum |
| `dhourmax` | Multi-day hourly maximum |
| `dhourrange` | Multi-day hourly range |
| `dhoursum` | Multi-day hourly sum |
| `dhourmean` | Multi-day hourly mean |
| `dhouravg` | Multi-day hourly average |
| `dhourstd` | Multi-day hourly standard deviation |
| `dhourstd1` | Multi-day hourly standard deviation (n-1) |
| `dhourvar` | Multi-day hourly variance |
| `dhourvar1` | Multi-day hourly variance (n-1) |
| `dminutemin` | Multi-day by the minute minimum |
| `dminutemax` | Multi-day by the minute maximum |
| `dminuterange` | Multi-day by the minute range |
| `dminutesum` | Multi-day by the minute sum |
| `dminutemean` | Multi-day by the minute mean |
| `dminuteavg` | Multi-day by the minute average |
| `dminutestd` | Multi-day by the minute standard deviation |
| `dminutestd1` | Multi-day by the minute standard deviation (n-1) |
| `dminutevar` | Multi-day by the minute variance |
| `dminutevar1` | Multi-day by the minute variance (n-1) |

#### Multi-year Daily Statistics

| Operator | Description |
|----------|-------------|
| `ydaymin` | Multi-year daily minimum |
| `ydaymax` | Multi-year daily maximum |
| `ydayrange` | Multi-year daily range |
| `ydaysum` | Multi-year daily sum |
| `ydaymean` | Multi-year daily mean |
| `ydayavg` | Multi-year daily average |
| `ydaystd` | Multi-year daily standard deviation |
| `ydaystd1` | Multi-year daily standard deviation (n-1) |
| `ydayvar` | Multi-year daily variance |
| `ydayvar1` | Multi-year daily variance (n-1) |
| `ydaypctl` | Multi-year daily percentile |

#### Multi-year Monthly Statistics

| Operator | Description |
|----------|-------------|
| `ymonmin` | Multi-year monthly minimum |
| `ymonmax` | Multi-year monthly maximum |
| `ymonrange` | Multi-year monthly range |
| `ymonsum` | Multi-year monthly sum |
| `ymonmean` | Multi-year monthly mean |
| `ymonavg` | Multi-year monthly average |
| `ymonstd` | Multi-year monthly standard deviation |
| `ymonstd1` | Multi-year monthly standard deviation (n-1) |
| `ymonvar` | Multi-year monthly variance |
| `ymonvar1` | Multi-year monthly variance (n-1) |
| `ymonpctl` | Multi-year monthly percentile |

#### Multi-year Seasonal Statistics

| Operator | Description |
|----------|-------------|
| `yseasmin` | Multi-year seasonal minimum |
| `yseasmax` | Multi-year seasonal maximum |
| `yseasrange` | Multi-year seasonal range |
| `yseassum` | Multi-year seasonal sum |
| `yseasmean` | Multi-year seasonal mean |
| `yseasavg` | Multi-year seasonal average |
| `yseasstd` | Multi-year seasonal standard deviation |
| `yseasstd1` | Multi-year seasonal standard deviation (n-1) |
| `yseasvar` | Multi-year seasonal variance |
| `yseasvar1` | Multi-year seasonal variance (n-1) |
| `yseaspctl` | Multi-year seasonal percentile |

#### Multi-year Daily Running Statistics

| Operator | Description |
|----------|-------------|
| `ydrunmin` | Multi-year daily running minimum |
| `ydrunmax` | Multi-year daily running maximum |
| `ydrunsum` | Multi-year daily running sum |
| `ydrunmean` | Multi-year daily running mean |
| `ydrunavg` | Multi-year daily running average |
| `ydrunstd` | Multi-year daily running standard deviation |
| `ydrunstd1` | Multi-year daily running standard deviation (n-1) |
| `ydrunvar` | Multi-year daily running variance |
| `ydrunvar1` | Multi-year daily running variance (n-1) |
| `ydrunpctl` | Multi-year daily running percentile |

### Correlation and Regression

| Operator | Description |
|----------|-------------|
| `fldcor` | Correlation in grid space |
| `timcor` | Correlation over time |
| `fldcovar` | Covariance in grid space |
| `timcovar` | Covariance over time |
| `regres` | Regression |
| `detrend` | Detrend time series |
| `trend` | Trend of time series |
| `addtrend` | Add trend |
| `subtrend` | Subtract trend |

### EOFs

| Operator | Description |
|----------|-------------|
| `eof` | Calculate EOFs in spatial or time space |
| `eoftime` | Calculate EOFs in time space |
| `eofspatial` | Calculate EOFs in spatial space |
| `eof3d` | Calculate 3-Dimensional EOFs in time space |
| `eofcoeff` | Principal coefficients of EOFs |

### Interpolation

| Operator | Description |
|----------|-------------|
| `remapbil` | Bilinear interpolation |
| `genbil` | Generate bilinear interpolation weights |
| `remapbic` | Bicubic interpolation |
| `genbic` | Generate bicubic interpolation weights |
| `remapknn` | k-nearest neighbor remapping |
| `remapnn` | Nearest neighbor remapping |
| `remapdis` | Distance weighted average remapping |
| `genknn` | Generate k-nearest neighbor remap weights |
| `gennn` | Generate nearest neighbor remap weights |
| `gendis` | Generate distance weighted average remap weights |
| `remapcon` | First order conservative remapping |
| `gencon` | Generate 1st order conservative remap weights |
| `remaplaf` | Largest area fraction remapping |
| `genlaf` | Generate largest area fraction remap weights |
| `remap` | Grid remapping (with pre-computed weights) |
| `remapeta` | Remap vertical hybrid levels |
| `ml2pl` | Model to pressure level interpolation |
| `ap2pl` | Vertical pressure interpolation |
| `gh2hl` | Geometric height interpolation |
| `intlevel` | Linear level interpolation |
| `intlevel3d` | Linear level interpolation from/to 3D vertical coordinates |
| `inttime` | Interpolation between timesteps |
| `intntime` | Interpolation between n timesteps |
| `intyear` | Interpolation between two years |

### Transformation

| Operator | Description |
|----------|-------------|
| `sp2gp` | Spectral to gridpoint |
| `gp2sp` | Gridpoint to spectral |
| `sp2sp` | Spectral to spectral (truncation change) |
| `sp2gpl` | Spectral to Gaussian gridpoint (linear) |
| `dv2ps` | Divergence and vorticity to velocity potential and stream function |
| `dv2uv` | Divergence and vorticity to U and V wind |
| `uv2dv` | U and V wind to divergence and vorticity |
| `fourier` | Fourier transformation |

### Import / Export

| Operator | Description |
|----------|-------------|
| `import_binary` | Import binary data sets |
| `import_cmsaf` | Import CM-SAF HDF5 files |
| `input` | ASCII input |
| `inputsrv` | SERVICE ASCII input |
| `inputext` | EXTRA ASCII input |
| `output` | ASCII output |
| `outputf` | Formatted output |
| `outputint` | Integer output |
| `outputsrv` | SERVICE ASCII output |
| `outputext` | EXTRA ASCII output |
| `outputtab` | Table output |
| `gmtxyz` | GMT xyz format |
| `gmtcells` | GMT multiple segment format |

### Miscellaneous

| Operator | Description |
|----------|-------------|
| `gradsdes` | GrADS data descriptor file |
| `after` | ECHAM standard post processor |
| `bandpass` | Bandpass filtering |
| `lowpass` | Lowpass filtering |
| `highpass` | Highpass filtering |
| `gridarea` | Grid cell area |
| `gridweights` | Grid cell weights |
| `smooth` | Smooth grid points |
| `smooth9` | 9 point smoothing |
| `deltat` | Difference between timesteps |
| `setvals` | Set list of old values to new values |
| `setrtoc` | Set range to constant |
| `setrtoc2` | Set range to constant, others to constant2 |
| `gridcellindex` | Get grid cell index |
| `const` | Create a constant field |
| `random` | Create a field with random numbers |
| `topo` | Create a field with topography |
| `seq` | Create a time series |
| `stdatm` | Create values for hydrostatic atmosphere |
| `timsort` | Temporal sorting |
| `uvDestag` | Destaggering of u/v wind components |
| `rotuvNorth` | Rotate u/v wind to North pole |
| `projuvLatLon` | Cylindrical Equidistant projection |
| `rotuvb` | Backward wind rotation |
| `mrotuvb` | Backward rotation of MPIOM data |
| `mastrfu` | Mass stream function |
| `pressure_half` | Pressure on half-levels |
| `pressure` | Pressure on full-levels |
| `delta_pressure` | Pressure difference of half-levels |
| `sealevelpressure` | Sea level pressure |
| `gheight` | Geopotential height on full-levels |
| `gheight_half` | Geopotential height on half-levels |
| `air_density` | Air density |
| `adisit` | Potential temperature to in-situ temperature |
| `adipot` | In-situ temperature to potential temperature |
| `rhopot` | Calculates potential density |
| `histcount` | Histogram count |
| `histsum` | Histogram sum |
| `histmean` | Histogram mean |
| `histfreq` | Histogram frequency |
| `sethalo` | Set the bounds of a field |
| `wct` | Windchill temperature |
| `fdns` | Frost days where no snow index per time period |
| `strwin` | Strong wind days index per time period |
| `strbre` | Strong breeze days index per time period |
| `strgal` | Strong gale days index per time period |
| `hurr` | Hurricane days index per time period |
| `cmorlite` | CMOR lite |
| `verifygrid` | Verify grid coordinates |
| `hpdegrade` | Degrade HEALPix grid |
| `hpupgrade` | Upgrade HEALPix grid |
| `symmetrize` | Mirrors data at the equator |
| `uv2vr_cfd` | U and V wind to relative vorticity (curvilinear finite differences) |
| `uv2dv_cfd` | U and V wind to divergence (curvilinear finite differences) |
| `cmor` | Climate Model Output Rewriting (CMIP-compliant) |

### Climate Indices (ETCCDI / ECA&D)

| Operator | Description |
|----------|-------------|
| `eca_cdd` / `etccdi_cdd` | Consecutive dry days index per time period |
| `eca_cfd` | Consecutive frost days index per time period |
| `eca_csu` | Consecutive summer days index per time period |
| `eca_cwd` / `etccdi_cwd` | Consecutive wet days index per time period |
| `eca_cwdi` | Cold wave duration index wrt mean of reference period |
| `eca_cwfi` / `etccdi_csdi` | Cold-spell days index wrt 10th percentile of reference period |
| `eca_etr` | Intra-period extreme temperature range |
| `eca_fd` / `etccdi_fd` | Frost days index per time period |
| `eca_gsl` | Thermal growing season length index |
| `eca_hd` | Heating degree days per time period |
| `eca_hwdi` | Heat wave duration index wrt mean of reference period |
| `eca_hwfi` / `etccdi_wsdi` | Warm spell days index wrt 90th percentile of reference period |
| `eca_id` / `etccdi_id` | Ice days index per time period |
| `eca_r75p` | Moderate wet days wrt 75th percentile of reference period |
| `eca_r75ptot` | Precipitation percent due to R75p days |
| `eca_r90p` | Wet days wrt 90th percentile of reference period |
| `eca_r90ptot` | Precipitation percent due to R90p days |
| `eca_r95p` | Very wet days wrt 95th percentile of reference period |
| `eca_r95ptot` | Precipitation percent due to R95p days |
| `eca_r99p` | Extremely wet days wrt 99th percentile of reference period |
| `eca_r99ptot` | Precipitation percent due to R99p days |
| `eca_pd` / `etccdi_r1mm` | Precipitation days index per time period |
| `eca_r10mm` | Heavy precipitation days index per time period |
| `eca_r20mm` | Very heavy precipitation days index per time period |
| `eca_rr1` | Wet days index per time period |
| `eca_rx1day` / `etccdi_rx1day` | Highest one day precipitation amount per time period |
| `eca_rx5day` / `etccdi_rx5day` | Highest five-day precipitation amount per time period |
| `eca_sdii` | Simple daily intensity index per time period |
| `eca_su` / `etccdi_su` | Summer days index per time period |
| `eca_tg10p` | Cold days percent wrt 10th percentile of reference period |
| `eca_tg90p` | Warm days percent wrt 90th percentile of reference period |
| `eca_tn10p` | Cold nights percent wrt 10th percentile of reference period |
| `eca_tn90p` | Warm nights percent wrt 90th percentile of reference period |
| `eca_tr` / `etccdi_tr` | Tropical nights index per time period |
| `eca_tx10p` | Very cold days percent wrt 10th percentile of reference period |
| `eca_tx90p` | Very warm days percent wrt 90th percentile of reference period |
| `etccdi_tx90p` | Percentage of days when daily max temperature > 90th percentile |
| `etccdi_tx10p` | Percentage of days when daily max temperature < 10th percentile |
| `etccdi_tn90p` | Percentage of days when daily min temperature > 90th percentile |
| `etccdi_tn10p` | Percentage of days when daily min temperature < 10th percentile |

## API Reference

### `Cdo` class

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cdo_path` | `str` | `None` | Path to CDO binary (auto-discovered if None) |
| `options` | `str` | `""` | Default CDO options (e.g. `"-O -s -f nc4"`) |
| `env` | `dict` | `None` | Custom environment variables |
| `debug` | `bool` | `False` | Print commands and stderr |
| `timeout` | `int` | `None` | Default timeout in seconds |

### Operator method parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| First positional arg | `str` | Operator parameters (e.g. `"0,30,0,30"`) |
| `input` | `str` or `list` | Input file(s) |
| `output` | `str` | Output file path |
| `options` | `str` | Additional CDO options for this call |
| `timeout` | `int` | Timeout override |
| `returnXArray` | `bool` | Return as xarray.Dataset |
| `returnCdf` | `bool` | Return as netCDF4.Dataset |
| `returnArray` | `bool` | Return as numpy.ndarray |
| `returnMaArray` | `bool` | Return as masked numpy.ndarray |

### Utility methods

| Method | Returns | Description |
|--------|---------|-------------|
| `cdo.help()` | `str` | General usage summary |
| `cdo.help("operator")` | `str` | CDO help text for a specific operator |
| `cdo.version()` | `str` | CDO version string |
| `cdo.operators()` | `set` | All available CDO operator names |
| `cdo.has_operator(name)` | `bool` | Check if an operator exists |
| `cdo.cleanup()` | — | Remove temporary files |

### `CdoError` exception

Import: `from skyborn_cdo import CdoError`

| Attribute | Type | Description |
|-----------|------|-------------|
| `returncode` | `int` | CDO exit code |
| `stderr` | `str` | Full CDO error output |
| `cmd` | `str` | Command that failed |

## Development

```bash
git clone --recurse-submodules https://github.com/QianyeSu/skyborn-cdo.git
cd skyborn-cdo
pip install -e ".[test]"
pytest tests/
```

## Author

**Qianye Su**  
Email: suqianye2000@gmail.com  
GitHub: [@QianyeSu](https://github.com/QianyeSu)

## License

This Python wrapper is licensed under **BSD-3-Clause**.

CDO itself is licensed under **BSD-3-Clause** by MPI für Meteorologie.
See [LICENSE](LICENSE) for details.
