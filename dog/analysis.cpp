//////////////////////////////////////////////////////////////////////////////
// analysis.C: implement Analyser class for NS-type solvers.
//
// This deals with output of runtime information such as step numbers,
// CFL estimation, modal energies, etc.  If set, also output history
// point and particle track information.
//
// It is assumed that the first 2 or 3 (for 3D) entries in the Domain
// u vector are velocity fields.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <unistd.h>
#include <sem.h>


Analyser::Analyser (Domain* D   ,
		    FEML*   file) :
// ---------------------------------------------------------------------------
// Set up.
//
// Fluid particle files (session.par) are dealt with here.
// Each line is of form
// #     tag  time  x     y      z
//       1    0.0   1.0   10.0   0.5.
// Output is of the same form, called session.trk.
//
// NB: Particle tracking is broken for multiprocessor application.
//
// History points are also set up here.  They are nominated in the
// optional HISTORY section of the session file.  Output is to
// session.his.
// ---------------------------------------------------------------------------
  _src (D)
{
  const char routine[] = "Analyser::Analyser";
  char       str[StrMax];
  time_t     tp (time (0));

  cout << setprecision (6);

  // -- Set up for history points: open files, create points.

  if (file -> seek ("HISTORY")) {
    int_t          i, id, num = 0;
    const int_t    NH = file -> attribute ("HISTORY", "NUMBER");
    const Element* E;
    HistoryPoint*  H;
    real_t         r, s, x, y, z;
    
    for (i = 0; i < NH; i++) {
      file -> stream() >> id >> x >> y >> z;
      if ((E = HistoryPoint::locate (x, y, D -> elmt, r, s))) {
	H = new HistoryPoint (id, E, r, s, x, y, z);
	_history.insert (_history.end(), H);
	num++;
      } else {
	sprintf (str, "History point at (%f, %f, %f) not in mesh", x, y, z);
	Veclib::alert (routine, str, WARNING);
      }
    }
    
    _his_strm.open (strcat (strcpy (str, _src -> name), ".his"));
    _his_strm.setf (ios::scientific, ios::floatfield);
    _his_strm.precision (6);
    if (!_his_strm) Veclib::alert (routine, "can't open history file", ERROR);
  }

  // -- Set up for output of modal energies every IO_CFL steps.

  _mdl_strm.open (strcat (strcpy (str, _src -> name), ".mdl"), ios::out); 
  _mdl_strm << "#     Time Mode         Energy" << endl
	    << "# ----------------------------" << endl;

  // -- Dump run information to file.

  ofstream runfile (strcat (strcpy (str, _src -> name), ".run"), ios::out);
  gethostname (str, StrMax);
  runfile << "-- Host                    : " << str << endl;
  runfile << "   PID                     : " << getpid() << endl;

  strftime (str, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  runfile << "   Date                    : " << str << endl;

  D -> report (runfile);
  
  runfile.close();
}


void Analyser::analyse (AuxField** work,
			AuxField** not_used)
// ---------------------------------------------------------------------------
// Step-by-step processing.  If SPAWN was set, add more particles at
// original absolute positions.
// ---------------------------------------------------------------------------
{
  const int_t cflstep = Femlib::ivalue ("IO_CFL");

  // -- Run information update.

  cout << "Step: " << _src -> step << "  Time: " << _src -> time << endl;

  // -- CFL, energy, divergence information.

  if (cflstep && !(_src -> step % cflstep)) {
    modalEnergy ();
    //    estimateCFL (work[0]);  // -- No real point in doing this.
    divergence  (work);
  }

  // -- Periodic dumps and global information.
  
  const bool periodic = !(_src -> step %  Femlib::ivalue("IO_HIS")) ||
                        !(_src -> step %  Femlib::ivalue("IO_FLD")) ;
  const bool final    =   _src -> step == Femlib::ivalue("N_STEP");
  const bool state    = periodic || final;

  if (state) {

    // -- Output history point data.
      
     int_t    i, j;
    const int_t       NH = _history.size();
    const int_t       NF = _src -> u.size();
    HistoryPoint*     H;
    vector<real_t>    tmp (NF);
    vector<AuxField*> u   (NF);

    for (i = 0; i < NF; i++)
      u[i] = _src -> u[i];

    for (i = 0; i < NH; i++) {
      H = _history[i];

      H -> extract (u, &tmp[0]);

      _his_strm << setw(4) << H->ID() << " " << setw(14) << _src->time << " ";
      for (j = 0; j < NF; j++) _his_strm << setw(15) << tmp[j];
      _his_strm << endl;
    }
  }

  _src -> dump();
}


void Analyser::modalEnergy ()
// ---------------------------------------------------------------------------
// Print out modal energies per unit area, output by root processor.
// ---------------------------------------------------------------------------
{
  const int_t NC = Geometry::nPert();
  real_t      ek = 0.0;

  for (int_t i = 0; i < NC; i++) ek += _src -> u[i] -> mode_L2 (0);

  _mdl_strm << setw(10) << _src -> time 
	    << setw( 5) << 1
	    << setw(15) << ek
	    << endl;
}


void Analyser::divergence (AuxField** Us) const
// ---------------------------------------------------------------------------
// Print out the velocity field's divergence energy per unit area.  Us
// is used as work area.
// ---------------------------------------------------------------------------
{
  const char routine[] = "Analyser::divergence";
  const int_t NC = Geometry::nPert();
  int_t       i;
  real_t      L2;

  if (Geometry::system() == Geometry::Cartesian) {
    for (i = 0; i < NC; i++) {
      *Us[i] = *_src -> u[i];
      Us[i] -> gradient (i);
    }
  } else {
    for (i = 0; i < NC; i++) *Us[i] = *_src -> u[i];
    Us[1] -> mulY();
    for (i = 0; i < NC; i++)  Us[i] -> gradient (i);
    Us[1] -> divY();
    if (NC == 3) Us[2] -> divY();
  }

  if (Geometry::problem() == Geometry::O2_3D_SYMM) *Us[2] *= -1.0;

  for (i = 1; i < NC; i++) *Us[0] += *Us[i];

  L2 = Us[0] -> mode_L2 (0);

  cout << "-- Divergence Energy: " << L2 << endl;

  // -- Crash stop.

  if (L2 != L2) Veclib::alert (routine, "forcing termination on NaN.", ERROR);
}


void Analyser::estimateCFL (AuxField* work) const
// ---------------------------------------------------------------------------
// Estimate and print the peak CFL number, based on zero-mode velocities.
// ---------------------------------------------------------------------------
{
  const int_t    pid     = Geometry::procID();
  const int_t    nProc   = Geometry::nProc();
  const real_t   dt      = Femlib::value ("D_T");
  real_t         CFL_dt, dt_max;
  int_t          i, percent;
  char           vcmpt;
  real_t         CFL_i[3];
  real_t         cmpt_i;
  int_t          elmt_i, elmt_j, elmt_k;

  CFL_i[0] = (*work = *_src -> u[0]) . CFL(0, elmt_i);
  CFL_i[1] = (*work = *_src -> u[1]) . CFL(1, elmt_j);

  CFL_dt = max(CFL_i[0], CFL_i[1]);
  cmpt_i = (CFL_i[0] > CFL_i[1]) ? 0.0 : 1.0;
  elmt_i = (CFL_i[0] > CFL_i[1]) ? elmt_i : elmt_j;

  if (_src -> nField() > 3) {
    CFL_i[2] = (*work = *_src -> u[2]) . CFL(2, elmt_k);

    if (CFL_i[2] > CFL_dt) {
      CFL_dt = CFL_i[2];
      cmpt_i = 2.0;
      elmt_i = elmt_k;
    }
  }
  cout << setprecision (3)
       << "# CFL: "     << CFL_dt * dt
       << ", dt (max): " << dt_max
       << ", dt (set): " << dt
       << " ("           << percent
       << "%), field: "  << vcmpt 
       << ", elmt: "     << elmt_i + 1 << endl;
  // -- 1-based indexing as in session file.
}
