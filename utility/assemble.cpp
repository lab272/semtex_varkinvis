/*****************************************************************************
 * assemble: diagnostic utility to generate global mesh numbering from
 * session file.
 *
 * NB: now that this numbering is computed at runtime by top-level
 * codes such as dns, this utility is, as stated, diagnostic, and its
 * output is not required by those codes.  Its output format matches
 * that of the old enumerate utility, to which it can be compared for
 * checking purposes.
 *
 * Usage
 * -----
 * assemble [options] file
 *   options:
 *   -h       ... display this message
 *   -v       ... set verbose output
 *   -n N     ... override element order to be N
 *   -O [0-3] ... set level of bandwidth optimization (default is 3)
 *
 * Synopsis
 * --------
 * Determine, from BCs section of FEML file, list of fields for which
 * numbering schemes are to be constructed.
 *
 * Generate a BC mask and initial numbering scheme for first field,
 * using Mesh class routines.  Optimize numbering scheme according to
 * selected level.
 *
 * For each succeeding field, first generate a BC mask and, if it
 * matches a mask previously generated, add the field's name to the
 * previous field's name vector but take no further action.
 * Otherwise, generate and optimize a new numbering system.
 *
 * Print up the masks and numbering schemes on cout.
 *
 * @file utility/assemble.cpp
 * @ingroup group_utility
 *****************************************************************************/

#include <cstdlib>
#include <cstring>
#include <climits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

using namespace std;

#include <cfemdef.h>
#include <femlib.h>
#include <utility.h>
#include <veclib.h>

#include <geometry.h>
#include <mesh.h>
#include <assemblymap.h>
#include <numbersys.h>

static char        prog[] = "assemble";
static int_t       verb   = 0;
static const int_t FldMax = 16;

static void getargs   (int, char**, char*&, int_t&, int_t&, int_t&);
static char axial     (FEML*);
static void getfields (FEML*, char*);
static void checkVBCs (FEML*, const char*);
static void checkABCs (FEML*, const char);
static void printup   (const char*, vector<AssemblyMap*>&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// See synopsis in file header.
// ---------------------------------------------------------------------------
{
  char  *session = 0, field[StrMax], axistag, name;
  FEML* file;
  Mesh* mesh;
  int_t mode, strat, i, j, k = 0, nz, nel, np = 0, opt = 0;
  bool  found, cylind;

  Geometry::CoordSys   space;
  AssemblyMap*         N;
  vector<AssemblyMap*> allMappings; // Complete set of AssemblyMaps for domain.
  
  Femlib::init ();
  
  getargs (argc, argv, session, verb, np, opt);
  
  if (verb) Femlib::ivalue ("VERBOSE", verb);

  file = new FEML (session);
  mesh = new Mesh (file);

  if (np) Femlib::ivalue ("N_P", np); // -- Over-ride on command-line request.
  
  nz     = Femlib::ivalue ("N_Z");
  np     = Femlib::ivalue ("N_P");
  cylind = static_cast<bool>(Femlib::ivalue ("CYLINDRICAL"));
  space  = (cylind) ? Geometry::Cylindrical : Geometry::Cartesian;
  nel    = mesh -> nEl();  

  Geometry::set (np, nz, nel, space);

  axistag = axial (file) && cylind && Geometry::nDim() > 2;

  getfields (file, field);
  
  if (axistag)                        checkABCs (file, axistag);
  if (cylind && Geometry::nDim() > 2) checkVBCs (file, field);

  vector<int_t> bmap (Geometry::nBnode());
  vector<int_t> mask (Geometry::nBnode());

  mesh -> buildAssemblyMap (np, &bmap[0]); // -- Naive version.
  strat = Femlib::ivalue ("ENUMERATION");

  // -- What follows inside the test is a cut+paste copy of code in
  //    Domain::makeAssemblyMaps().
  
  if (axistag) {

    // -- (Cylindrical coords.) Number systems are different for modes
    //    0, 1, 2+ owing to presence of axial BCs (Blackburn & Sherwin
    //    JCP 197 2004).

    for (i = 0; i < strlen(field); i++) {
      name = field[i];
      
      for (mode = 0; mode < 3; mode++) {

	mesh -> buildLiftMask (Geometry::nP(), name, mode, &mask[0]);
	
	for (found = false, j = 0; !found && j < allMappings.size(); j++)
	  if (found = allMappings[j] -> willMatch (mask)) {
	    allMappings[j] -> addTag (name, mode);
	    break;
	  }
	if (!found) {
	  N = new AssemblyMap (Geometry::nP(), Geometry::nElmt(),
			       strat, bmap, mask, name, mode);
	  allMappings.push_back (N);
	}
      }
    }
    
  } else {
    
    // -- Number systems for all Fourier modes are identical.
    
    for (i = 0; i < strlen (field); i++) {
      name = field[i];

      mesh -> buildLiftMask (Geometry::nP(), name, 0, &mask[0]);

      for (found = false, j = 0; !found && j < allMappings.size(); j++)
	if (found = allMappings[j] -> willMatch (mask)) {
	  allMappings[j] -> addTag (name, 0);
	  allMappings[j] -> addTag (name, 1);
	  allMappings[j] -> addTag (name, 2);
	  break;
	}
      if (!found) {
	N = new AssemblyMap (Geometry::nP(), Geometry::nElmt(),
			     strat, bmap, mask, name, 0);
	N -> addTag (name, 1);
	N -> addTag (name, 2);

	allMappings.push_back (N);
      }
    }
  }

  printup (field, allMappings);

  return EXIT_SUCCESS;
}


static void getargs (int    argc   , 
		     char** argv   ,
		     char*& session,
		     int_t& verb   ,
		     int_t& np     ,
		     int_t& opt    )
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] = "usage: assemble [options] session\n"
                 "options:\n"
                 "  -h       ... display this message\n"
                 "  -v       ... set verbose output\n"
		 "  -n N     ... override number of element knots to be N\n"
		 "  -O [0-3] ... bandwidth optimization level [Default: 3]\n";
  char err[StrMax];
  char c;

  while (--argc && **++argv == '-')
    switch (c = *++argv[0]) {
    case 'h':
      cerr << usage;
      exit (EXIT_SUCCESS);
      break;
    case 'v':
      for (verb = 1; *++argv[0] == 'v'; verb++);
      break;
    case 'n':
      if (*++argv[0])
	np = atoi (*argv);
      else {
	--argc;
	np = atoi (*++argv);
      }
      break;
    case 'O':
      if (*++argv[0])
	opt = atoi (*argv);
      else {
	--argc;
	opt = atoi (*++argv);
      }
      break;
    default:
      sprintf (err, "getargs: illegal option: %c\n", c);
      Veclib::alert (prog, err, ERROR);
      break;
    }

  if   (argc == 1) session = *argv;
  else             Veclib::alert (prog, "must provide session file", ERROR);
}


static char axial (FEML* file)
// ---------------------------------------------------------------------------
// Return the character tag corresponding to "axis" group, if that exists.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("GROUPS")) return '\0';

  int_t       i;
  char        nextc, buf[StrMax];
  const int_t N (file->attribute ("GROUPS", "NUMBER"));
  
  for (i = 0; i < N; i++) {
    while ((nextc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');
    file->stream() >> buf >> nextc >> buf;
    if (strstr (buf, "axis")) return nextc;
  }

  return '\0';
}


static void getfields (FEML* file ,
		       char* field)
// ---------------------------------------------------------------------------
// Set up the list of fields according to the names found in the
// 'FIELDS' section of the FEML file, e.g.
//
// <FIELDS>
// # optional comment lines...
//   u v w c p
// </FIELDS>
// ---------------------------------------------------------------------------
{
  int_t i = 0;
  char  t[StrMax];

  // -- Set up string for the fields listed.

  if (file->seek ("FIELDS")) {

    file->stream().ignore (StrMax, '\n');
    
    while (file->stream().peek() == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    do {
      file->stream() >> field[i++];
    } while (field[i - 1] != '<' && i < StrMax);

    if (field[--i] == '<') {
      field[i] = '\0';
      file->stream() >> t;
      if (!(strstr (t,   "/FIELDS")))
	   Veclib::alert (prog, "FIELDS section not closed", ERROR);
    } else Veclib::alert (prog, "FIELDS section not closed", ERROR);
  } else   Veclib::alert (prog, "FIELDS section not found",  ERROR);
}


static void checkVBCs (FEML*       file ,
		       const char* field)
// ---------------------------------------------------------------------------
// For cylindrical 3D fluids problems, the declared boundary condition types
// for velocity fields v & w must be the same for all groups, to allow
// for coupling of these fields (which uncouples the viscous substep).
//
// Check it out by running through each group's BCs and checking for 
// tag agreement on v & w BCs.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("BCS")) return;
  if (!strchr (field, 'u')) return;
  if (!strchr (field, 'v') || !strchr (field, 'w')) return;

  int_t       i, j, id, nbcs;
  char        vtag, wtag, groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;
    vtag = wtag = '\0';

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	Veclib::alert (prog, err, ERROR);
      }

      file->stream() >> fieldc;
      if      (fieldc == 'v') vtag = tagc;
      else if (fieldc == 'w') wtag = tagc;
      file->stream().ignore (StrMax, '\n');
    }
    
    if (!(vtag && wtag)) {
      sprintf (err, "group %c: BCs for fields 'v' & 'w' needed", groupc);
      Veclib::alert (prog, err, ERROR);
    }
    if (vtag != wtag) {
      sprintf (err, "group %c, fields 'v' & 'w': BC type mismatch", groupc);
      Veclib::alert (prog, err, ERROR);
    }
  }
}


static void checkABCs (FEML*      file ,
		       const char atag)
// ---------------------------------------------------------------------------
// Run through and ensure that for "axis" group, all BCs are of type <A>.
// ---------------------------------------------------------------------------
{
  if (!file->seek ("BCS")) return;

  int_t       i, j, id, nbcs;
  char        groupc, fieldc, tagc, tag[StrMax], err[StrMax];
  const int_t N (file->attribute ("BCS", "NUMBER"));

  for (i = 0; i < N; i++) {

    while ((groupc = file->stream().peek()) == '#') // -- Skip comments.
      file->stream().ignore (StrMax, '\n');

    file->stream() >> id >> groupc >> nbcs;

    for (j = 0; j < nbcs; j++) {

      file->stream() >> tag;
      if (strchr (tag, '<') && strchr (tag, '>') && (strlen (tag) == 3))
	tagc = tag[1];
      else {
	sprintf (err, "unrecognized BC tag format: %s", tag);
	Veclib::alert (prog, err, ERROR);
      }

      file->stream() >> fieldc;

      if (groupc == atag && tagc != 'A') {
	sprintf (err, "group '%c': field '%c' needs axis BC", groupc, fieldc);
	Veclib::alert (prog, err, ERROR);
      }
      file->stream().ignore (StrMax, '\n');
    }
  }
}


void printup (const char*           F,
	      vector<AssemblyMap*>& S)
// ---------------------------------------------------------------------------
// Print up summary info followed by map & mask for each system.
// ---------------------------------------------------------------------------
{
  int_t       i, j, k, side, soff;
  const int_t nedge = Geometry::nP() - 1;
  const int_t nSys  = S.size();
  char        mField[StrMax];
  
  cout << "# FIELDS         :  " << F << endl;

  cout << "# ----------------";
  for (j = 0; j < nSys; j++) cout << "  -------------";
  cout << endl;

  cout << "# " << nSys << " NUMBER SETS  :";
  for (j = 0; j < nSys; j++) {
    S[j] -> printTags (mField);
    cout << setw (15);
    cout << mField;
  }
  cout << endl;

  cout << "# NEL            :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nElmt();
  }
  cout << endl;
  
  cout << "# NP_MAX         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nP();
  }
  cout << endl;
  
  cout << "# NEXT_MAX       :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nExtElmt();
  }
  cout << endl;
  
  cout << "# NINT_MAX       :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nIntElmt();
  }
  cout << endl;
  
  cout << "# NTOTAL         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nPlane();
  }
  cout << endl;
  
  cout << "# NBOUNDARY      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Geometry::nBnode();
  }
  cout << endl;

  cout << "# NGLOBAL        :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nGlobal();
  }
  cout << endl;

  cout << "# NSOLVE         :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nSolve();
  }
  cout << endl;

  cout << "# OPTIMIZATION   :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << Femlib::ivalue ("ENUMERATION");
  }
  cout << endl;

  cout << "# BANDWIDTH      :";
  for (j = 0; j < nSys; j++) {
    cout << setw (15);
    cout << S[j] -> nBand();
  }
  cout << endl;

  cout << "# ----------------";
  for (j = 0; j < nSys; j++) cout << "  -------------";
  cout << endl;

  cout << "# elmt  side offst";
  for (j = 0; j < nSys; j++) cout << "     bmap  mask";
  cout << endl;

  for (i = 0, k = 1; k <= Geometry::nElmt(); k++)
    for (side = 1; side <= 4; side++)
      for (soff = 0; soff < nedge; soff++, i++) {
	cout << setw (6) << k << setw (6) << side << setw (6) << soff;
	for (j = 0; j < nSys; j++)
	  cout 
	    << setw (9) << S[j] -> btog ()[i]
	    << setw (6) << S[j] -> bmask()[i];
	cout << endl;
      }
}
