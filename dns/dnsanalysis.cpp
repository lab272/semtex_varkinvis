///////////////////////////////////////////////////////////////////////////////
// This version of analysis.cpp is specialized so that it computes and
// prints out forces (and if appropriate, integrated scalar flux)
// exerted on "wall" boundary group.
//
// Copyright (c) 1994+, Hugh M Blackburn
///////////////////////////////////////////////////////////////////////////////

#include <dns.h>

static bool WALLED;


DNSAnalyser::DNSAnalyser (Domain* D   ,
			  BCmgr*  B   ,
			  FEML*   feml) :
// ---------------------------------------------------------------------------
// Extensions to Analyser class.
// ---------------------------------------------------------------------------
  Analyser (D, feml),
  _wss (Femlib::ivalue ("IO_WSS") && B -> nWall())
{
  const char routine[] = "DNSAnalyser::DNSAnalyser";
  char       str[StrMax];

  WALLED = B -> nWall() > 0; // -- Only true if we have a "wall" group.

  if (WALLED) {
    ROOTONLY {
      // -- Open state-variable file.

      _flx_strm.open (strcat (strcpy (str, _src -> name), ".flx"));
      if (!_flx_strm) Veclib::alert (routine, "can't open flux file", ERROR);

      _flx_strm << "# DNS state information file"      << endl;
      if (_src -> hasScalar()) {
	_flx_strm << "# Step Time Flux [Fpre Fvis Ftot]-axis" << endl;
	_flx_strm << "# ------------------------------------" << endl;

      } else {
	_flx_strm << "# Step Time [Fpre Fvis Ftot]-axis" << endl;
	_flx_strm << "# -------------------------------" << endl;
      }
    }

    if (_wss) {
      // -- Set up to compute wall shear stresses.    
    
      const int_t npr = Geometry::nProc();
      const int_t np  = Geometry::nP();
      const int_t nz  = Geometry::nZProc();

      // -- Allocate storage area: 3 = 1 normal component + 2 tangential.
    
      _nwall = B -> nWall();
      _nline = np * _nwall;
      _npad  = 3  * _nline;

      // -- Round up length for Fourier transform/exchange.

      if   (npr > 1) _npad += 2 * npr - _npad % (2 * npr);
      else           _npad += _npad % 2;

      _work.resize (_npad * nz);

      // -- Open file.

      ROOTONLY {
	_wss_strm.open (strcat (strcpy (str, _src -> name), ".wss"));
	if (!_wss_strm) Veclib::alert (routine, "can't open WSS file", ERROR);
      }
    }
  }
}


void DNSAnalyser::analyse (AuxField** work0,
			   AuxField** work1)
// ---------------------------------------------------------------------------
// Step-by-step processing.  The default action is to carry out a standard
// analysis.  Then, if we have a SURFACE that is in a "wall" GROUP, output wall
// momentum/scalar flux information, either in integral or pointwise form.
//
// Thie routine is designed to be called at the end of a time step, when
// the domain D (or _src) has the velocity (scalar) and pressure data.
// The AuxField storage arrays (typically, the top level of Us and Uf)
// are available as scratch.  
// ---------------------------------------------------------------------------
{

  Analyser::analyse (work0, work1);

  if (WALLED) {
    
    const char  routine[] = "DNSAnalyser::analyse";
    const int_t NVEL = _src -> nVelCmpt();
    const int_t NADV = _src -> nAdvect();
    bool        periodic = !(_src->step %  Femlib::ivalue ("IO_HIS")) ||
                           !(_src->step %  Femlib::ivalue ("IO_FLD"));
    bool        final    =   _src->step == Femlib::ivalue ("N_STEP");
    bool        state    = periodic || final;

    if (state) ROOTONLY {
	real_t flux;
	Vector pfor, vfor, tfor;
	char   s[StrMax];

	if (NVEL == 3) {
	  pfor   = Field::normTraction (_src -> u[NADV]);
	  vfor   = Field::tangTraction (_src -> u[0], _src -> u[1], _src->u[2]);
	  tfor.x = pfor.x + vfor.x;
	  tfor.y = pfor.y + vfor.y;
	  tfor.z = pfor.z + vfor.z;
	} else {
	  pfor   = Field::normTraction (_src -> u[NADV]);
	  vfor   = Field::tangTraction (_src -> u[0], _src -> u[1]);
	  tfor.x = pfor.x + vfor.x;
	  tfor.y = pfor.y + vfor.y;
	  tfor.z = pfor.z = vfor.z = 0.0;
	}

	if (_src -> hasScalar()) {
	  flux = Field::scalarFlux (_src -> u[NVEL]);
	  sprintf (s,
		   "%6d %#10.6g %10.6g "
		   "%#10.6g %#10.6g %#10.6g "
		   "%#10.6g %#10.6g %#10.6g "
		   "%#10.6g %#10.6g %#10.6g",
		   _src -> step, _src -> time, flux,
		   pfor.x,   vfor.x,   tfor.x,
		   pfor.y,   vfor.y,   tfor.y,
		   pfor.z,   vfor.z,   tfor.z);
	} else {
	  sprintf (s,
		   "%6d %#10.6g "
		   "%#10.6g %#10.6g %#10.6g "
		   "%#10.6g %#10.6g %#10.6g "
		   "%#10.6g %#10.6g %#10.6g",
		   _src -> step, _src -> time,
		   pfor.x,   vfor.x,   tfor.x,
		   pfor.y,   vfor.y,   tfor.y,
		   pfor.z,   vfor.z,   tfor.z);
	}

	_flx_strm << s << endl;
      }

    if (_wss) {
      periodic = !(_src->step % Femlib::ivalue ("IO_WSS")) ||
	!(_src->step % Femlib::ivalue ("IO_FLD")) ;
      state    = periodic || final;

      if (state) {
	const int_t    nP  = Geometry::nP();
	const int_t    nZ  = Geometry::nZ();
	const int_t    nZP = Geometry::nZProc();
	const int_t    nPR = Geometry::nProc();
	const int_t    nPP = _npad / nPR;
	int_t          i, j, k;
	real_t*        plane;
	vector<real_t> buffer (_nline);

	// -- Load the local storage area.

	Veclib::zero (_work.size(), &_work[0], 1);

	if (NVEL == 3)
	  Field::traction (&_work[0], &_work[_nline], &_work[2*_nline],
			   _nwall, _npad,
			   _src->u[NADV], _src->u[0], _src->u[1], _src->u[2]);
	else
	  Field::traction (&_work[0], &_work[_nline], &_work[2*_nline],
			   _nwall, _npad,
			   _src->u[NADV], _src->u[0], _src->u[1]);

	// -- Inverse Fourier transform (like Field::bTransform).

	if (nPR == 1) {
	  if (nZ > 1)
	    if (nZ == 2)
	      Veclib::copy (_npad, &_work[0], 1, &_work[_npad], 1);
	    else
	      Femlib::DFTr (&_work[0], nZ, _npad, INVERSE);
	} else {
	  Message::exchange (&_work[0], nZP, _npad, FORWARD);
	  Femlib::DFTr      (&_work[0], nZ,    nPP, INVERSE);
	  Message::exchange (&_work[0], nZP, _npad, INVERSE);
	}

	// -- Write to file.

	// -- Header: this will be a lot like a standard header.
	//    Output normal and tangential tractions, 'n', 't', 's'.

	ROOTONLY {
	  const char *Hdr_Fmt[] = { 
	    "%-25s "                "Session\n",
	    "%-25s "                "Created\n",
	    "%-5d1    %-5d %-10d"   "Nr, Ns, Nz, Elements\n",
	    "%-25d "                "Step\n",
	    "%-25.6g "              "Time\n",
	    "%-25.6g "              "Time step\n",
	    "%-25.6g "              "Kinvis\n",
	    "%-25.6g "              "Beta\n",
	    "%-25s "                "Fields written\n",
	    "%-25s "                "Format\n"
	  };

	  char   s1[StrMax], s2[StrMax];
	  time_t tp (time (0));

	  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
	
	  sprintf (s1, Hdr_Fmt[0], _src->name);               _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[1], s2);                       _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[2], nP, nZ, _nwall);           _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[3], _src->step);               _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[4], _src->time);               _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[5], Femlib::value ("D_T"));    _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[6], Femlib::value ("KINVIS")); _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[7], Femlib::value ("BETA"));   _wss_strm << s1;
	  sprintf (s1, Hdr_Fmt[8], "nts");                     _wss_strm << s1;
	  sprintf (s2, "binary "); Veclib::describeFormat  (s2 + strlen (s2));
	  sprintf (s1, Hdr_Fmt[9], s2);                       _wss_strm << s1;

	  if (!_wss_strm) Veclib::alert (routine,
					 "failed writing WSS header", ERROR);
	  _wss_strm << flush;
	}

	// -- Data.

	if (nPR > 1) {		// -- Parallel.
	  for (j = 0; j < 3; j++)	// -- Reminder: there are 3 components.
	    ROOTONLY {
	      for (i = 0; i < nZP; i++) {
		plane = &_work[i*_npad + j*_nline];
		_wss_strm.write(reinterpret_cast<char*>(plane),
				static_cast<int_t>(_nline * sizeof (real_t))); 
		if (_wss_strm.bad())
		  Veclib::alert (routine,
				 "unable to write binary output", ERROR);
	      }
	      for (k = 1; k < nPR; k++)
		for (i = 0; i < nZP; i++) {
		  Message::recv (&buffer[0], _nline, k);
		  _wss_strm.write(reinterpret_cast<char*>(&buffer[0]),
				  static_cast<int_t>(_nline * sizeof(real_t))); 
		  if (_wss_strm.bad())
		    Veclib::alert (routine,
				   "unable to write binary output", ERROR);
		}
	      _wss_strm << flush;
	    } else			// -- Not on root process.
		for (i = 0; i < nZP; i++) {
		  plane = &_work[i*_npad + j*_nline];
		  Message::send (plane, _nline, 0);
		}
	} else {			// -- Serial.
	  for (j = 0; j < 3; j++)
	    for (i = 0; i < nZ; i++) {
	      plane = &_work[i*_npad + j*_nline];
	      _wss_strm.write (reinterpret_cast<char*>(plane),
			       static_cast<int_t>(_nline * sizeof (real_t))); 
	      if (_wss_strm.bad())
		Veclib::alert (routine,
			       "unable to write binary output", ERROR);
	    }
	  _wss_strm << flush;
	}
      }
    }
  }
}
