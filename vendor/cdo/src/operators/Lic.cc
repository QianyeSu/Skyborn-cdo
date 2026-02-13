/*
  This file is part of CDO. CDO is a collection of Operators to manipulate and analyse Climate model Data.

  Author: Uwe Schulzweida

*/

#include <cstdlib>
#include <cdi.h>

#include "c_wrapper.h"
#include "cdo_options.h"
#include "process_int.h"
#include "color.h"

CPT cpt;

// http://www.zhanpingliu.org/research/flowvis/LIC/LIC_Source.htm

///////////////////////////////////////////////////////////////////////////
//               Line Integral Convolution for Flow Visualization
//                                                   Initial  Version
//
//                                                     May 15, 1999
//
//                                                                  by
//
//                             Zhanping Liu
//
//                      (zhanping@erc.msstate.edu)
//                                                          while with
//
//                                             Graphics  Laboratory
//
//                                              Peking  University
//
//                                                     P. R.  China
//
//----------------------------------------------------------------------///
//                                                   Later  Condensed
//
//                             May 4,  2002
//
//          VAIL (Visualization Analysis & Imaging Laboratory)
//                  ERC  (Engineering Research Center)
//                     Mississippi State University
///////////////////////////////////////////////////////////////////////////

#define SQUARE_FLOW_FIELD_SZ 400
#define DISCRETE_FILTER_SIZE 2048 * 2
#define LOWPASS_FILTR_LENGTH 10.00000f
#define LINE_SQUARE_CLIP_MAX 100000.0f
#define VECTOR_COMPONENT_MIN 0.050000f

///		normalize the vector field     ///
static void
NormalizVectrs(int n_xres, int n_yres, float *pVectr)
{
  for (int j = 0; j < n_yres; ++j)
    for (int i = 0; i < n_xres; ++i)
    {
      int index = (j * n_xres + i) << 1;
      float vcMag = (float) (std::sqrt((double) (pVectr[index] * pVectr[index] + pVectr[index + 1] * pVectr[index + 1])));

      float scale = (vcMag == 0.0f) ? 0.0f : 1.0f / vcMag;
      pVectr[index] *= scale;
      pVectr[index + 1] *= scale;
    }
}

///		make white noise as the LIC input texture     ///
static void
MakeWhiteNoise(int n_xres, int n_yres, std::vector<unsigned char> &pNoise)
{
  for (int j = 0; j < n_yres; ++j)
    for (int i = 0; i < n_xres; ++i)
    {
      int r = std::rand();
      r = ((r & 0xff) + ((r & 0xff00) >> 8)) & 0xff;
      pNoise[j * n_xres + i] = (unsigned char) r;
    }
}

///		generate box filter LUTs     ///
void
GenBoxFiltrLUT(int LUTsiz, std::vector<float> &p_LUT0, std::vector<float> &p_LUT1)
{
  for (int i = 0; i < LUTsiz; ++i) p_LUT0[i] = p_LUT1[i] = i;
}

// write the LIC image to a PPM file
/*
static void
WriteImage2PPM(int n_xres, int n_yres, std::vector<unsigned char> &pImage, const char *f_name)
{
  std::FILE *o_file;
  if ((o_file = std::fopen(f_name, "w")) == nullptr)
    {
      printf("Can't open output file\n");
      return;
    }

  std::fprintf(o_file, "P6\n%d %d\n255\n", n_xres, n_yres);

  for (int j = 0; j < n_yres; ++j)
    for (int i = 0; i < n_xres; ++i)
      {
        unsigned char unchar = pImage[j * n_xres + i];
        std::fprintf(o_file, "%c%c%c", unchar, unchar, unchar);
      }

  std::fclose(o_file);
}
*/

// write the LIC image to a PNG file

#ifdef HAVE_LIBPNG
#include <png.h>

static void
setRGBmag(png_byte *ptr, float val, float mag)
{
  int r = 0, g = 0, b = 0, n;

  for (n = 0; n < cpt.ncolors; ++n)
    if (mag > cpt.lut[n].z_low && mag <= cpt.lut[n].z_high) break;

  if (n == cpt.ncolors)
  {
    r = cpt.bfn[0].rgb[0];
    g = cpt.bfn[0].rgb[1];
    b = cpt.bfn[0].rgb[2];
  }
  else
  {
    r = cpt.lut[n].rgb_high[0];
    g = cpt.lut[n].rgb_high[1];
    b = cpt.lut[n].rgb_high[2];
  }

  float foregroundAlpha = .8;
  r = (val * foregroundAlpha) + (r * (1.0 - foregroundAlpha));
  g = (val * foregroundAlpha) + (g * (1.0 - foregroundAlpha));
  b = (val * foregroundAlpha) + (b * (1.0 - foregroundAlpha));
  ptr[0] = r;
  ptr[1] = g;
  ptr[2] = b;
}

static int
WriteImage2PNG(int width, int height, std::vector<unsigned char> const &pImage, float *mag, const char *filename)
{
  int code = 0;
  png_structp png_ptr = nullptr;
  png_infop info_ptr = nullptr;
  png_bytep row = nullptr;

  // Open file for writing (binary mode)
  auto fobj = c_fopen(filename, "wb");
  if (fobj.get() == nullptr)
  {
    std::fprintf(stderr, "Could not open file %s for writing\n", filename);
    code = 1;
    goto finalise;
  }

  // Initialize write structure
  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr);
  if (png_ptr == nullptr)
  {
    std::fprintf(stderr, "Could not allocate write struct\n");
    code = 1;
    goto finalise;
  }

  // Initialize info structure
  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == nullptr)
  {
    std::fprintf(stderr, "Could not allocate info struct\n");
    code = 1;
    goto finalise;
  }

  // Setup Exception handling
  if (setjmp(png_jmpbuf(png_ptr)))
  {
    std::fprintf(stderr, "Error during png creation\n");
    code = 1;
    goto finalise;
  }

  png_init_io(png_ptr, fobj.get());

  // Write header (8 bit colour depth)
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE,
               PNG_FILTER_TYPE_BASE);

  // Set title
  /*
  if (title != nullptr) {
          png_text title_text;
          title_text.compression = PNG_TEXT_COMPRESSION_NONE;
          title_text.key = "Title";
          title_text.text = title;
          png_set_text(png_ptr, info_ptr, &title_text, 1);
  }
  */
  png_write_info(png_ptr, info_ptr);

  // Allocate memory for one row (3 bytes per pixel - RGB)
  std::vector<png_bytep> row(3 * width);

  // Write image data
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++) { setRGBmag(&(row[x * 3]), pImage[y * width + x], mag[y * width + x]); }
    png_write_row(png_ptr, row.data());
  }

  // End write
  png_write_end(png_ptr, nullptr);

finalise:
  if (info_ptr != nullptr) png_free_data(png_ptr, info_ptr, PNG_FREE_ALL, -1);
  if (png_ptr != nullptr) png_destroy_write_struct(&png_ptr, (png_infopp) nullptr);

  return code;
}
#endif

///		flow imaging (visualization) through Line Integral Convolution     ///
static void
FlowImagingLIC(int n_xres, int n_yres, float *pVectr, std::vector<unsigned char> const &pNoise, std::vector<unsigned char> &pImage,
               std::vector<float> &p_LUT0, std::vector<float> &p_LUT1, float krnlen)
{

  /// for each pixel in the 2D output LIC image///

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (int j = 0; j < n_yres; ++j)
  {
    int vec_id;                       /// ID in the VECtor buffer (for the input flow field)
    int advDir;                       /// ADVection DIRection (0: positive;  1: negative)
    int advcts;                       /// number of ADVeCTion stepS per direction (a step counter)
    int ADVCTS = (int) (krnlen * 3);  /// MAXIMUM number of advection steps per direction to break dead loops

    float vctr_x;     /// x-component  of the VeCToR at the forefront point
    float vctr_y;     /// y-component  of the VeCToR at the forefront point
    float clp0_x;     /// x-coordinate of CLiP point 0 (current)
    float clp0_y;     /// y-coordinate of CLiP point 0	(current)
    float clp1_x;     /// x-coordinate of CLiP point 1 (next   )
    float clp1_y;     /// y-coordinate of CLiP point 1 (next   )
    float samp_x;     /// x-coordinate of the SAMPle in the current pixel
    float samp_y;     /// y-coordinate of the SAMPle in the current pixel
    float tmpLen;     /// TeMPorary LENgth of a trial clipped-segment
    float segLen;     /// SEGment   LENgth
    float curLen;     /// CURrent   LENgth of the streamline
    float prvLen;     /// PReVious  LENgth of the streamline
    float W_ACUM;     /// ACcuMulated Weight from the seed to the current streamline forefront
    float texVal;     /// TEXture VALue
    float smpWgt;     /// WeiGhT of the current SaMPle
    float t_acum[2];  /// two ACcUMulated composite Textures for the two directions, perspectively
    float w_acum[2];  /// two ACcUMulated Weighting values   for the two directions, perspectively
    float len2ID = (DISCRETE_FILTER_SIZE - 1) / krnlen;  /// map a curve LENgth TO an ID in the LUT
    for (int i = 0; i < n_xres; ++i)
    {
      /// init the composite texture accumulators and the weight accumulators///
      t_acum[0] = t_acum[1] = w_acum[0] = w_acum[1] = 0.0f;

      /// for either advection direction///
      for (advDir = 0; advDir < 2; advDir++)
      {
        /// init the step counter, curve-length measurer, and streamline seed///
        advcts = 0;
        curLen = 0.0f;
        clp0_x = i + 0.5f;
        clp0_y = j + 0.5f;

        /// access the target filter LUT///
        /// WeiGhT Look Up Table pointing to the target filter LUT
        auto const &wgtLUT = (advDir == 0) ? p_LUT0 : p_LUT1;

        /// until the streamline is advected long enough or a tightly  spiralling center / focus is encountered///
        while (curLen < krnlen && advcts < ADVCTS)
        {
          /// access the vector at the sample///
          vec_id = ((int) (clp0_y) *n_xres + (int) (clp0_x)) << 1;
          vctr_x = pVectr[vec_id];
          vctr_y = pVectr[vec_id + 1];

          /// in case of a critical point///
          if (vctr_x == 0.0f && vctr_y == 0.0f)
          {
            t_acum[advDir] = (advcts == 0) ? 0.0f : t_acum[advDir];  /// this line is indeed unnecessary
            w_acum[advDir] = (advcts == 0) ? 1.0f : w_acum[advDir];
            break;
          }

          /// negate the vector for the backward-advection case///
          vctr_x = (advDir == 0) ? vctr_x : -vctr_x;
          vctr_y = (advDir == 0) ? vctr_y : -vctr_y;

          /// clip the segment against the pixel boundaries --- find the shorter from the two clipped segments///
          /// replace  all  if-statements  whenever  possible  as  they  might  affect the computational speed///
          segLen = LINE_SQUARE_CLIP_MAX;
          segLen = (vctr_x < -VECTOR_COMPONENT_MIN) ? ((int) (clp0_x) -clp0_x) / vctr_x : segLen;
          segLen = (vctr_x > VECTOR_COMPONENT_MIN) ? ((int) ((int) (clp0_x) + 1.5f) - clp0_x) / vctr_x : segLen;
          segLen = (vctr_y < -VECTOR_COMPONENT_MIN) ? (((tmpLen = ((int) (clp0_y) -clp0_y) / vctr_y) < segLen) ? tmpLen : segLen)
                                                    : segLen;
          segLen = (vctr_y > VECTOR_COMPONENT_MIN)
                       ? (((tmpLen = ((int) ((int) (clp0_y) + 1.5f) - clp0_y) / vctr_y) < segLen) ? tmpLen : segLen)
                       : segLen;

          /// update the curve-length measurers///
          prvLen = curLen;
          curLen += segLen;
          segLen += 0.0004f;

          /// check if the filter has reached either end///
          segLen = (curLen > krnlen) ? ((curLen = krnlen) - prvLen) : segLen;

          /// obtain the next clip point///
          clp1_x = clp0_x + vctr_x * segLen;
          clp1_y = clp0_y + vctr_y * segLen;

          /// obtain the middle point of the segment as the texture-contributing sample///
          samp_x = (clp0_x + clp1_x) * 0.5f;
          samp_y = (clp0_y + clp1_y) * 0.5f;

          /// obtain the texture value of the sample///
          if (samp_x >= n_xres) samp_x = n_xres - 1;  // bug fix
          if (samp_y >= n_yres) samp_y = n_yres - 1;  // bug fix
          texVal = pNoise[(int) (samp_y) *n_xres + (int) (samp_x)];

          /// update the accumulated weight and the accumulated composite texture (texture x weight)///
          W_ACUM = wgtLUT[(int) (curLen * len2ID)];
          smpWgt = W_ACUM - w_acum[advDir];
          w_acum[advDir] = W_ACUM;
          t_acum[advDir] += texVal * smpWgt;

          /// update the step counter and the "current" clip point///
          advcts++;
          clp0_x = clp1_x;
          clp0_y = clp1_y;

          /// check if the streamline has gone beyond the flow field///
          if (clp0_x < 0.0f || clp0_x >= n_xres || clp0_y < 0.0f || clp0_y >= n_yres) break;
        }
      }

      /// normalize the accumulated composite texture///
      texVal = (t_acum[0] + t_acum[1]) / (w_acum[0] + w_acum[1]);

      /// clamp the texture value against the displayable intensity range [0, 255]
      texVal = (texVal < 0.0f) ? 0.0f : texVal;
      texVal = (texVal > 255.0f) ? 255.0f : texVal;
      pImage[j * n_xres + i] = (unsigned char) texVal;
    }
  }
}

static void
lic1(const char *obasename, int num, int n_xres, int n_yres, float *pVectr)
{
  std::vector<float> p_LUT0(DISCRETE_FILTER_SIZE);
  std::vector<float> p_LUT1(DISCRETE_FILTER_SIZE);
  std::vector<unsigned char> pNoise(n_xres * n_yres);
  std::vector<unsigned char> pImage(n_xres * n_yres);

  std::vector<float> mag(n_xres * n_yres);
  int n = n_xres * n_yres;

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (int i = 0; i < n; ++i)
    mag[i] = (float) std::sqrt((double) pVectr[i * 2] * pVectr[i * 2] + (double) pVectr[i * 2 + 1] * pVectr[i * 2 + 1]);

  if (Options::cdoVerbose)
  {
    float vmin = 1.e33;
    float vmax = -1.e33;
    for (int i = 0; i < n_xres * n_yres; ++i)
    {
      if (mag[i] < vmin) vmin = mag[i];
      if (mag[i] > vmax) vmax = mag[i];
    }
    cdo_print("ts=%d minval=%g  maxval=%g", num + 1, vmin, vmax);
  }

  NormalizVectrs(n_xres, n_yres, pVectr);
  MakeWhiteNoise(n_xres, n_yres, pNoise);
  GenBoxFiltrLUT(DISCRETE_FILTER_SIZE, p_LUT0, p_LUT1);
  FlowImagingLIC(n_xres, n_yres, pVectr, pNoise, pImage, p_LUT0, p_LUT1, LOWPASS_FILTR_LENGTH);
  // const char *fname = "LIC.ppm";
  // WriteImage2PPM(n_xres, n_yres, pImage, fname);
  char filename[8192];
  std::snprintf(filename, sizeof(filename), "%s%04d.png", obasename, num);

#ifdef HAVE_LIBPNG
  WriteImage2PNG(n_xres, n_yres, pImage.data(), mag.data(), filename);
#else
  cdo_warning("PNG support not compiled in!");
#endif
}

class Lic : public Process
{
public:
  using Process::Process;
  inline static CdoModule module = {
    .name = "Lic",
    .operators = { { "lic" } },
    .aliases = {},
    .mode = EXPOSED,     // Module mode: 0:intern 1:extern
    .number = CDI_REAL,  // Allowed number type
    .constraints = { 1, 1, NoRestriction },
  };
  inline static RegisterEntry<Lic> registration = RegisterEntry<Lic>();

  int numLevels = 0;
  size_t numMissVals;

  CdoStreamID streamID1;
  CdoStreamID streamID2;

  int varID1 = 0;
  int varID2 = 1;

  bool invertlat;

  size_t nx;
  size_t ny;

  size_t gridsize;

  int gridID;
  int nstart;

  const char *obasename;

  Varray<double> xvals;
  Varray<double> yvals;

public:
  void
  init() override
  {
    // int LIC  = cdo_operator_add("lic",  0, 0, nullptr);

    // int operatorID = cdo_operator_id();

    if (cdo_operator_argc() < 1 || cdo_operator_argc() > 2) operator_check_argc(1);
    auto const &cptFile = cdo_operator_argv(0);
    nstart = (cdo_operator_argc() == 2) ? std::stoi(cdo_operator_argv(1)) : 0;

    {
      auto fobj = c_fopen(cptFile, "r");
      if (fobj.get() == nullptr) cdo_abort("Open failed on color palette table %s", cptFile);

      auto status = cpt_read(fobj.get(), &cpt);
      if (status != 0) cdo_abort("Error during read of color palette table %s", cptFile);
    }

    obasename = cdo_get_stream_name(1);
    streamID1 = cdo_open_read(0);

    auto vlistID1 = cdo_stream_inq_vlist(streamID1);
    VarList varList1(vlistID1);

    // int taxisID1 = vlistInqTaxis(vlistID1);

    /* find variables */
    auto numVars = varList1.numVars();
    if (numVars != 2) cdo_abort("Need 2 input variable, found %d\n", numVars);

    auto numGrids = vlistNumGrids(vlistID1);
    if (numGrids != 1) cdo_abort("Need 1 grid, found %d\n", numGrids);

    gridID = varList1.vars[varID1].gridID;
    // int zaxisID = varList1.vars[varID1].zaxisID;

    gridsize = varList1.gridsizeMax();

    Varray<double> ivar1, ivar2, ovar1;
    if (varID1 != -1 && varID2 != -1)
    {
      numLevels = varList1.vars[varID1].nlevels;
      gridsize = varList1.vars[varID1].gridsize;
      ivar1.resize(numLevels * gridsize);
      ivar2.resize(numLevels * gridsize);
      ovar1.resize(numLevels * gridsize);
    }

    auto gridtype = gridInqType(gridID);
    if (gridtype != GRID_LONLAT && gridtype != GRID_GAUSSIAN) cdo_abort("Unsupported grid!\n");

    nx = gridInqXsize(gridID);
    ny = gridInqYsize(gridID);

    xvals.resize(nx);
    yvals.resize(ny);
    gridInqYvals(gridID, yvals.data());
    invertlat = (yvals[0] < yvals[ny - 1]);
  }

  void
  run() override
  {
    Varray<double> array1(gridsize);
    std::vector<float> varray(2 * gridsize);

    int tsID = 0;
    while (true)
    {
      auto numFields = cdo_stream_inq_timestep(streamID1, tsID);
      if (numFields == 0) break;
      /*
      cdo_taxis_copy_timestep(taxisID2, taxisID1);
      cdo_def_timestep(streamID2, tsID);
      */
      for (int fieldID = 0; fieldID < numFields; ++fieldID)
      {
        auto [varID, levelID] = cdo_inq_field(streamID1);

        if (varID == varID1 || varID == varID2)
        {
          cdo_read_field(streamID1, array1.data(), &numMissVals);
          if (numMissVals) cdo_abort("Missing values unsupported!");

          gridsize = gridInqSize(gridID);

          if (invertlat)
          {
            for (size_t j = 0; j < ny / 2; ++j)
            {
              for (size_t i = 0; i < nx; ++i) xvals[i] = array1[j * nx + i];
              for (size_t i = 0; i < nx; ++i) array1[j * nx + i] = array1[(ny - j - 1) * nx + i];
              for (size_t i = 0; i < nx; ++i) array1[(ny - j - 1) * nx + i] = xvals[i];
            }
          }

          if (varID == varID1)
            for (size_t i = 0; i < gridsize; ++i) varray[i * 2] = (float) array1[i];
          else if (varID == varID2)
            for (size_t i = 0; i < gridsize; ++i) varray[i * 2 + 1] = (float) array1[i];
        }
      }

      lic1(obasename, nstart + tsID, nx, ny, varray.data());

      tsID++;
    }
  }

  void
  close() override
  {
    cdo_stream_close(streamID1);
  }
};

// mencoder
// png + alpha canel libpng
// ls -1v | grep JPG > files.txt
// mencoder -nosound -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=21600000 -o windowsill_flowers_7.avi -mf type=jpeg:fps=24
// mf://@files.txt -vf scale=1920:1080 convert imageA.png imageB.png -alpha off -compose CopyOpacity -composite out.png
// makecpt -Cjet -T-10/30/.2  > cjet.cpt
// makecpt -Chot -T0/30/.2 -I > chot.cpt
// ffmpeg -r 1/5 -start_number 0 -i "image%4d.png" -c:v libx264 -video_size 4k -vf "fps=25,format=yuv420p" movie.mp4
