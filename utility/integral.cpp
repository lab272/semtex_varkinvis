/*****************************************************************************
 * integral: utility which returns the domain integral of all 
 * fields in dump file.
 *
 * Usage
 * -----
 * integral [-h] [-v] [-c] session [file]
 *
 * Synopsis 
 * --------
 * Read in file, first print up area of domain.  If 3D perform Fourier
 * transform to get mean value into plane zero for each field.  Then
 * return integral (and centroidal x,y locations) for each scalar
 * field.  For 3D, values are multiplied by domain length, to produce
 * volume integrals of each scalar.
 *
 * If the coordinate system is cylindrical, then the integrals
 * (including the area) are weighted by the radius (hence domain
 * volume = "area" * TWOPI/BETA: "area" is the true area * centroidal
 * radius). Use -c switch to turn this off.
 *
 * Integrals are approximated by GLL quadrature at order supplied in
 * session file.
 *
 * @file utility/integral.cpp
 * @ingroup group_utility
 *****************************************************************************/
// Copyright (c) 1999+, Hugh M Blackburn

#include <sem.h>

static char  prog[]  = "integral";
static int_t verbose = 0;
static void  getargs  (int, char**, char*&, char*&, bool&);
static bool  getDump  (istream&, vector<AuxField*>&, vector<Element*>&,
		       const int_t, const int_t, const int_t);
static bool  doSwap   (const char*);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char               *session = 0, *dump = 0;
  istream            *fldfile;
  int_t              NP, NZ,  NEL;
  int_t              np, nel, ntot, i;
  real_t             Lz, Area = 0.0, integral;
  Vector             centroid;
  const real_t       *z;
  FEML*              F;
  Mesh*              M;
  Geometry::CoordSys space;
  bool               cylind = true;
  vector<Element*>   Esys;
  vector<AuxField*>  u;

  // -- Initialize.

  Femlib::init ();
  
  getargs (argc, argv, session, dump, cylind);
  
  cout.precision (8);

  if (dump) {
    fldfile = new ifstream (dump);
    if (fldfile -> fail()) Veclib::alert (prog, "no field file", ERROR);
  } else fldfile = &cin;

  // -- Set up 2D mesh information.
  
  F   = new FEML (session);
  M   = new Mesh (F);

  NEL = M -> nEl();  
  NP  = Femlib::ivalue ("N_P");
  NZ  = Femlib::ivalue ("N_Z");
  Lz  = (NZ > 1) ? Femlib::value ("TWOPI / BETA") : 1.;
  space = (Femlib::ivalue ("CYLINDRICAL") && cylind) ? 
    Geometry::Cylindrical : Geometry::Cartesian;

  Geometry::set (NP, NZ, NEL, space);
  Esys.resize   (NEL);

  // -- Note: for cylindrical coords, this "Area" is radius-weighted.
  for (i = 0; i < NEL; i++) {
    Esys[i] = new Element (i, NP, M);
    Area   += Esys[i] -> area();
  }
  cout << Area << endl;
  
  // -- Load field file, Gauss--Lobatto integrate all variables within it.

  while (getDump (*fldfile, u, Esys, NP, NZ, NEL)) {
    for (i = 0; i < u.size(); i++) {
      u[i] -> transform (FORWARD); // -- Go back to Fourier space.
      centroid = u[i] -> centroid (0);
      integral = u[i] -> integral (0);
      cout << u[i] -> name() << ": " << Lz * integral 
	   << " , centroid: " << centroid.x << " , " << centroid.y << endl;
    }
  }

  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session,
		     char*& dump   ,
		     bool&  cylind )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "Usage: integral [options] session [dump]\n"
    "options:\n"
    "-h ... print this message\n"
    "-v ... verbose output\n"
    "-c ... switch cylindrical coordinates off, if defined in session\n";
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cout << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'c':
      cylind = false;
      break;
    default:
      cerr << usage;
      exit (EXIT_FAILURE);
      break;
    }

  if      (argc == 1)   session = argv[0];
  else if (argc == 2) { session = argv[0]; dump = argv[1]; }
  else                  Veclib::alert (prog, usage, ERROR);
}


static bool getDump (istream&           file,
		     vector<AuxField*>& u   ,
		     vector<Element*>&  Esys,
		     const int_t        np  ,
		     const int_t        nz  ,
		     const int_t        nel )
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  char    buf[StrMax], fields[StrMax];
  int_t   i, swab, nf, npnew, nznew, nelnew;
  real_t* alloc;

  if (file.getline(buf, StrMax).eof()) return 0;
  
  if (!strstr (buf, "Session"))
    Veclib::alert (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> npnew >> nznew >> nznew >> nelnew;
  file.getline (buf, StrMax);
  
  if (np != npnew || nz != nznew || nel != nelnew)
    Veclib::alert (prog, "size of dump mismatch with session file", ERROR);

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  nf = strlen  (fields);
  file.getline (buf, StrMax);

  // -- Arrange for byte-swapping if required.

  file.getline  (buf, StrMax);
  swab = doSwap (buf);

  // -- Create AuxFields on first pass.

  if (u.size() == 0) {
    u.resize (nf);
    for (i = 0; i < nf; i++) {
      alloc = new real_t [Geometry::nTotProc()];
      u[i]  = new AuxField (alloc, nz, Esys, fields[i]);
    }
  } else if (u.size() != nf) 
    Veclib::alert
      (prog, "number of fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[i];
    if (swab) u[i] -> reverse();
  }

  return file.good();
}


static bool doSwap (const char* ffmt)
// ---------------------------------------------------------------------------
// Figure out if byte-swapping is required to make sense of binary input.
// ---------------------------------------------------------------------------
{
  char mfmt[StrMax];

  Veclib::describeFormat (mfmt);   

  if (!strstr (ffmt, "binary"))
    Veclib::alert (prog, "input field file not in binary format", ERROR);
  else if (!strstr (ffmt, "endian"))
    Veclib::alert (prog, "input field file in unknown binary format", WARNING);

  return (strstr (ffmt, "big") && strstr (mfmt, "little")) || 
         (strstr (mfmt, "big") && strstr (ffmt, "little"));
}
