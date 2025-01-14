/*****************************************************************************
 * stressdiv: from a field file containing Reynolds stresses, utility
 * computes their NEGATIVE divergence to give mean-flow forcing terms.
 *
 * Usage
 * -----
 * stressdiv [options] session session.rey
 *   options:
 *   -h        ... print this message
 *
 * Output to standard output.
 *
 * Synopsis
 * --------
 * NB: the input field file is assumed to contain only velocity and
 * appropriate Reynolds stress data (and so is an AVERAGE=2 .avg file
 * that has already been passed through rstress).  The naming
 * conventions employed for example in addfield.cpp are broken (names
 * are here used for different variables).
 *
 * Input field names
 * -----------------
 *
 * u -- x velocity component (cylindrical: axial)
 * v -- y velocity component (cylindrical: radial)
 * w -- z velocity component (cylindrical: azimuthal)
 * p -- pressure/density
 *
 * Names for components of the symmetric Reynolds stresses:
 *
 *                      / uu uv uw \     /  A  B  D \
 *                      | .  vv vw |  =  |  .  C  E |
 *                      \ .  .  ww /     \  .  .  F /
 *
 * Output field names
 * ------------------
 *
 * u -- x momentum equation forcing
 * v -- y momentum equation forcing
 * w -- z momentum equation forcing
 *
 * Cartesian coordinates/3D
 * ------------------------
 *
 * u = -[d(A)/dx + d(B)/dy + d(D)/dz]
 * v = -[d(B)/dx + d(C)/dy + d(E)/dz]
 * w = -[d(D)/dx + d(E)/dy + d(F)/dz]
 *
 * Cartesian coordinates/2D
 * ------------------------
 *
 * u = -[d(A)/dx + d(B)/dy]
 * v = -[d(B)/dx + d(C)/dy]
 *
 * Cylindrical coordinates/3D
 * --------------------------
 *
 * u = -[d(A)/dx + 1/y*d(y*B)/dy + 1/y*d(D)/dz]
 *   = -[d(A)/dx + d(B)/dy + B/y + 1/y*d(D)/dz]
 * v = -[d(B)/dx + 1/y*d(y*C)/dy + 1/y*d(E)/dz - F/y]
 *   = -[d(B)/dx + d(C)/dy + C/y + 1/y*d(E)/dz - F/y]
 * w = -[d(D)/dx + d(E)/dy       + 1/y*d(F)/dz + 2*E/y]
 *
 * Cylindrical coordinates/2D
 * --------------------------
 *
 * u = -[d(A)/dx + 1/y*d(y*B)/dy = d(A)/dx + d(B)/dy + B/y]
 * v = -[d(B)/dx + 1/y*d(y*C)/dy = d(B)/dx + d(C)/dy + C/y]
 *
 * The dimensionality of the problem is determined from
 *
 * 1. The number of z planes in the input data file (not the
 * session file). If this is 1 then only A, B, C are considered.
 *
 * 2. If nz > 1 then also the input must contain the terms D, E, F. In
 * this case, the value of BETA in the session file is used in
 * computing z-derivatives.
 *
 * NB. The output fields are smoothed (made C_0 across element boundaries).
 *
 * @file utility/stressdiv.cpp
 * @ingroup group_utility
 *****************************************************************************/
// Copyright (c) 2010+, Hugh M Blackburn

#include <sem.h>

static char  prog[] = "stressdiv";
static void  getargs (int,char**,const char*&,const char*&);
static int_t preScan (ifstream&, int_t&);
static void  getMesh (const char*,vector<Element*>&,
		      vector<int_t>&, vector<real_t>&, const int_t);
static void  makeBuf (map<char,AuxField*>&, AuxField*&, vector<Element*>&, 
		      const int_t);
static bool  getDump (ifstream&,map<char, AuxField*>&,vector<Element*>&);
static bool  doSwap  (const char*);
static void  stress  (map<char,AuxField*>&, map<char,AuxField*>&, AuxField*, 
		      const vector<int_t>&, const vector<real_t>&, const int);
static const char* fieldNames(map<char, AuxField*>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver -- adapted from eneq.cpp.
// ---------------------------------------------------------------------------
{
  const char           *session, *dump;
  ifstream             avgfile;
  vector<Element*>     elmt;
  map<char, AuxField*> input, output;
  vector<AuxField*>    outbuf;
  AuxField*            work;
  vector<int_t>        bmap;
  vector<real_t>       massInv;
  int_t                NDIM;
  int_t                NCOM;

  Femlib::init ();
  
  getargs (argc, argv, session, dump);

  avgfile.open (dump, ios::in);
  if (!avgfile) Veclib::alert (prog, "no field file", ERROR);
  NDIM = preScan (avgfile, NCOM);
  
  getMesh (session,      elmt, bmap, massInv, NDIM);
  makeBuf (output, work, elmt, NCOM);

  // -- Need to link outbuf to output so we can use writeField.
  //    Maybe in the longer term we should overload writeField.

  outbuf.resize (output.size()); int_t i = 0;
  for (map<char,AuxField*>::iterator k = output.begin();
       k != output.end(); k++, i++) outbuf[i] = k -> second;
  
  while (getDump (avgfile, input, elmt)) {
    stress     (input, output, work, bmap, massInv, NCOM);
    writeField (cout, session, 0, 0.0, outbuf);
  }
  
  return EXIT_SUCCESS;
}


static void getargs (int          argc   ,
		     char**       argv   ,
		     const char*& session,
		     const char*& dump   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] session dump.avg\n"
    "options:\n"
    "  -h ... print this message \n";
              
  char buf[StrMax];
 
  while (--argc && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
      break;
    }

  if (argc != 2) {
    sprintf (buf, usage, prog); cerr << buf; exit (EXIT_FAILURE);
  } else {
    --argc; session = *argv;
    --argc; dump    = *++argv;
  }
}


static int_t preScan (ifstream& file, int_t& ncom)
// ---------------------------------------------------------------------------
// Read ahead in the avg file and find the dimensionality of the
// stress (2 if nz = 1 or 3 if nz > 1).  Also check that the file
// contains the relevant Reynolds stress information. If OK, rewind
// and return, else die.
// ---------------------------------------------------------------------------
{
  char  buf[StrMax], fields[StrMax];
  int_t nz, ndim = 0;

  if (file.getline(buf, StrMax).eof()) 
    Veclib::alert (prog, "stress file is empty", ERROR);
  
  if (!strstr (buf, "Session")) 
    Veclib::alert (prog, "stress file not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> nz >> nz >> nz;
  file.getline (buf, StrMax);

  ndim = (nz == 1) ? 2 : 3;

  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);
  file.getline (buf, StrMax);

  // -- Input field names, assumed to be written without intervening spaces.

  file >> fields;
  file.getline (buf, StrMax);

  if (strstr(fields, "ABCDEF")) ncom = 3;
  else if (strstr(fields, "ABC")) ncom = 2;
  
  if (ndim == 3 and ncom != 3)
    Veclib::alert (prog,  "stress file must contain fields ABCDEF", ERROR);
    
  // -- All clear: rewind & return.

  file.seekg (0);

  return ndim;
}


static void getMesh (const char*       session,
		     vector<Element*>& elmt   ,
		     vector<int_t>&    bmap   ,
		     vector<real_t>&   imass  ,
		     const int_t       ndim   )
// ---------------------------------------------------------------------------
// Set up 2D mesh information. Note that parser tokens and Geometry
// are set here, too.  NDIM value from stress file over-rides session
// file, and from here forward, Geometry::nZ() can serve as a proxy
// for NDIM.
// ---------------------------------------------------------------------------
{
  FEML* F = new FEML (session);
  Mesh* M = new Mesh (F);

  const int_t nel = M -> nEl();  
  const int_t np  = Femlib::ivalue ("N_P");
  const int_t nz  = (ndim == 2) ? 1 : Femlib::ivalue ("N_Z");
  int_t       k, *gid;
  
  Geometry::CoordSys space = (Femlib::ivalue ("CYLINDRICAL")) ?
    Geometry::Cylindrical : Geometry::Cartesian;
  
  Geometry::set (np, nz, nel, space);
  elmt.resize   (nel);
  
  for (k = 0; k < nel; k++) elmt[k] = new Element (k, np, M);

  const int_t nbndry = Geometry::nBnode();
  
  bmap.resize (nbndry);
  M -> buildAssemblyMap (Geometry::nP(), &bmap[0]);
  const int_t nglobal = bmap[Veclib::imax (nbndry, &bmap[0], 1)] + 1;
  imass.resize (nglobal);

  const int_t npnp = Geometry::nTotElmt();
  const int_t next = Geometry::nExtElmt();
  
  vector <real_t> unity (npnp, 1.0);
  
  Veclib::zero (nglobal, &imass[0], 1);

  for (gid = &bmap[0], k = 0; k < nel; k++, gid += next)
    elmt[k] -> bndryDsSum (gid, &unity[0], &imass[0]);
  Veclib::vrecp (nglobal, &imass[0], 1, &imass[0], 1);
}


static void makeBuf (map<char, AuxField*>& output,
		     AuxField*&            work  ,
		     vector<Element*>&     elmt  ,
		     const int             ncom  )
// ---------------------------------------------------------------------------
// Note that we only set up the output and work buffers here. The
// input buffers get created in getDump.
// ---------------------------------------------------------------------------
{
  const int_t nz   = Geometry::nZ();
  const int_t ntot = Geometry::nTotal();

  if (ncom == 2) {			// -- Set by stress file info.
    output['u'] = new AuxField (new real_t[ntot], nz, elmt, 'u');
    output['v'] = new AuxField (new real_t[ntot], nz, elmt, 'v');
  } else {
    output['u'] = new AuxField (new real_t[ntot], nz, elmt, 'u');
    output['v'] = new AuxField (new real_t[ntot], nz, elmt, 'v');
    output['w'] = new AuxField (new real_t[ntot], nz, elmt, 'w');
  }

  work = new AuxField (new real_t[ntot], nz, elmt, '\0');
}


static const char* fieldNames (map<char, AuxField*>& u)
// ---------------------------------------------------------------------------
// Return string containing single-character names of fields.
// ---------------------------------------------------------------------------
{
  int_t i;
  static char buf[StrMax];
  map<char, AuxField*>::iterator k;

  for (i = 0, k = u.begin(); k != u.end(); k++, i++) buf[i] = k -> first;
  buf[i] = '\0';

  return buf;
}


static bool getDump (ifstream&             file,
		     map<char, AuxField*>& u   ,
		     vector<Element*>&     elmt)
// ---------------------------------------------------------------------------
// Load data from field dump, with byte-swapping if required.
// If there is more than one dump in file, it is required that the
// structure of each dump is the same as the first.
// ---------------------------------------------------------------------------
{
  const int_t ntot = Geometry::nTotal();
  char        buf[StrMax], fields[StrMax];
  int_t       i, nf, np, nz, nel;
  bool        swab;

  if (file.getline(buf, StrMax).eof()) return false;
  
  if (!strstr (buf, "Session")) Veclib::alert (prog, "not a field file", ERROR);
  file.getline (buf, StrMax);

  // -- Input numerical description of field sizes.

  file >> np >> nz >> nz >> nel;
  file.getline (buf, StrMax);
  
  if (np != Geometry::nP() || nz != Geometry::nZ() || nel != Geometry::nElmt())
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
    for (i = 0; i < nf; i++)
      u[fields[i]] = new AuxField (new real_t[ntot], nz, elmt, fields[i]);
  } else if (strcmp (fieldNames (u), fields) != 0)
    Veclib::alert (prog, "fields mismatch with first dump in file", ERROR);

  // -- Read binary field data.

  for (i = 0; i < nf; i++) {
    file >> *u[fields[i]];
    if (swab) u[fields[i]] -> reverse();
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


static void stress  (map<char, AuxField*>& in   ,
		     map<char, AuxField*>& out  ,
		     AuxField*             work ,
		     const vector<int_t>&  bmap ,
		     const vector<real_t>& imass,
		     const int             ncom )
// ---------------------------------------------------------------------------
// This does the actual work of building stress divergence terms.
//
// Note again that the input file is assumed to already have correctly
// reduced Reynolds stress terms in it (e.g. as computed by rstress).
//
// Cartesian coordinates/3D:
// ---------------------------
//
// u = -[d(A)/dx + d(B)/dy + d(D)/dz]
// v = -[d(B)/dx + d(C)/dy + d(E)/dz]
// w = -[d(D)/dx + d(E)/dy + d(F)/dz]
//
// Cartesian coordinates/2D:
// ---------------------------
//
// u = -[d(A)/dx + d(B)/dy]
// v = -[d(B)/dx + d(C)/dy]
//
// Cylindrical coordinates/3D:
// ---------------------------
//
// u = -[d(A)/dx + 1/y*d(y*B)/dy + 1/y*d(D)/dz]
//   = -[d(A)/dx + d(B)/dy + B/y + 1/y*d(D)/dz]
// v = -[d(B)/dx + 1/y*d(y*C)/dy + 1/y*d(E)/dz - F/y]
//   = -[d(B)/dx + d(C)/dy + C/y + 1/y*d(E)/dz - F/y]
// w = -[d(D)/dx + d(E)/dy       + 1/y*d(F)/dz + 2*E/y]
//
// TA: double checked second lines
//
// Cylindrical coordinates/2D:
// ---------------------------
//
// u = -[d(A)/dx + 1/y*d(y*B)/dy] = -[d(A)/dx + d(B)/dy + B/y]
// v = -[d(B)/dx + 1/y*d(y*C)/dy] = -[d(B)/dx + d(C)/dy + C/y]
// if 2D3C:
// w = -[d(D)/dx + d(E)/dy + 2*E/y]

// NB: Cylindrical/3D is implemented but yet to be validated!
// ---------------------------------------------------------------------------
{
  if (Geometry::cylindrical()) {
    if (Geometry::nZ() == 1) {
                                // -- 2D2C
      (*out['u'] = *in['A']) . gradient(0);
      *out['u'] += (*work = *in['B']) . gradient(1);
      (*out['v'] = *in['B']) . gradient(0);
      *out['v'] += (*work = *in['C']) . gradient(1);

      *out['u'] += (*work = *in['B']) . divY();
      *out['v'] += (*work = *in['C']) . divY();
      *out['v'] -= (*work = *in['F']) . divY();

      *out['u'] *= -1.0;
      *out['v'] *= -1.0;

      if (ncom == 3) {          // -- 2D3C
        (*out['w'] = *in['D']) . gradient(0);
        *out['w'] += (*work = *in['E']) . gradient(1);
        *work = (*in['E']) . divY();
        *out['w'] += (*work *= 2.);
        *out['w'] *= -1.0;
      }
    } else {			// -- 3D

      // -- TA: checked that stressdiv and project -w -z1 commute.

      // u  = -[d(A)/dx + d(B)/dy + B/y + 1/y*d(D)/dz]
      (*out['u'] = *in['D']).transform(FORWARD).gradient(2).transform(INVERSE);
      *out['u'] += *in['B'];
      out['u'] -> divY();
      *out['u'] += (*work = *in['B']) . gradient(1);
      *out['u'] += (*work = *in['A']) . gradient(0);
      *out['u'] *= -1.0;

      // v  = -[d(B)/dx + d(C)/dy + C/y + 1/y*d(E)/dz - F/y]
      (*out['v'] = *in['E']).transform(FORWARD).gradient(2).transform(INVERSE);
      *out['v'] -= *in['F'];
      *out['v'] += *in['C'];
      out['v'] -> divY();
      *out['v'] += (*work = *in['C']) . gradient(1);
      *out['v'] += (*work = *in['B']) . gradient(0);
      *out['v'] *= -1.0;

      // w = -[d(D)/dx + d(E)/dy       + 1/y*d(F)/dz + 2*E/y]
      (*out['w'] = *in['F']).transform(FORWARD).gradient(2).transform(INVERSE);
      *out['w'] += *in['E'];
      *out['w'] += *in['E'];
      out['w'] -> divY();
      *out['w'] += (*work = *in['E']) . gradient(1);
      *out['w'] += (*work = *in['D']) . gradient(0);
      *out['w'] *= -1.0;
    }
  } else {			// -- Cartesian.
    if (Geometry::nZ() == 1) {
                                // -- 2D2C
      (*out['u'] = *in['A']) . gradient(0);
      *out['u'] += (*work = *in['B']) . gradient(1);

      (*out['v'] = *in['B']) . gradient(0);
      *out['v'] += (*work = *in['C']) . gradient(1);

      *out['u'] *= -1.0;
      *out['v'] *= -1.0;
      if (ncom == 3) {          // -- 2D3C
        (*out['w'] = *in['D']) . gradient(0);
        *out['w'] += (*work = *in['E']) . gradient(1);
        *out['w'] *= -1.0;
      }
    } else {			// -- 3D
      (*out['u'] = *in['A']) . gradient(0);
      *out['u'] += (*work = *in['B']) . gradient(1);
      *out['u'] += (*work = *in['D']) .
	transform(FORWARD).gradient(2).transform(INVERSE);

      (*out['v'] = *in['B']) . gradient(0);
      *out['v'] += (*work = *in['C']) . gradient(1);
      *out['v'] += (*work = *in['E'])
	. transform(FORWARD).gradient(2).transform(INVERSE);

      (*out['w'] = *in['D']) . gradient(0);
      *out['w'] += (*work = *in['E']) . gradient(1);
      *out['w'] += (*work = *in['F'])
	. transform(FORWARD).gradient(2).transform(INVERSE);

      *out['u'] *= -1.0;
      *out['v'] *= -1.0;
      *out['w'] *= -1.0;
    }
  }

  out['u'] -> smooth (imass.size(), &bmap[0], &imass[0]);
  out['v'] -> smooth (imass.size(), &bmap[0], &imass[0]);
  if (ncom == 3) out['w'] -> smooth (imass.size(), &bmap[0], &imass[0]);
}

