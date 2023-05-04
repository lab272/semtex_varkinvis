/*****************************************************************************
 * calc.cpp: a simple calculator based on the femlib function parser.
 *
 * Usage
 * -----
 * calc [-h] [file]
 *
 * @file utility/calc.cpp
 * @ingroup group_utility
 *****************************************************************************/
// Copyright (c) 1994+, Hugh M Blackburn

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include <cfemdef.h>
#include <utility.h>
#include <veclib.h>
#include <femlib.h>

static char prog[] = "calc";
static void getargs (int, char**, istream*&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  char     buf[StrMax];
  istream* input;

  Femlib::init ();
  
  getargs (argc, argv, input);

  while (input -> getline (buf, FILENAME_MAX))
    if (strstr (buf, "="))
      Femlib::value (buf);
    else
      cout << setprecision(17) << Femlib::value (buf) << endl;
  
  return EXIT_SUCCESS;
}


static void getargs (int       argc ,
		     char**    argv ,
		     istream*& input)
// ---------------------------------------------------------------------------
// Parse command-line arguments.
// ---------------------------------------------------------------------------
{
  char buf[StrMax];
  char usage[] = "Usage: %s [-h] [file]\n";
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      cerr << "-- Preset internal variables:"  << endl;
      yy_show ();
      cerr << endl;
      cerr << "-- Calculator operators, functions and procedures:" << endl;
      yy_help ();
      exit (EXIT_SUCCESS);
      break;
    default:
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_FAILURE);
      break;
    }

  if (argc == 1) {
    input = new ifstream (*argv);
    if (input -> fail()) {
      cerr << usage;
      sprintf (buf, "unable to open file: %s", *argv);
      Veclib::messg (prog, buf, ERROR);
    }
  } else input = &cin;
}

