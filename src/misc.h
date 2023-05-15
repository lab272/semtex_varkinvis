#ifndef MISC_H
#define MISC_H

// char*    upperCase   (char *);

// -- Routines from misc.C:

ostream& printVector (ostream&, const char*, const int_t, ... );
void     readField   (istream&, vector<AuxField*>&);
void     writeField  (ostream&, const char*, const int_t, const real_t,
		      vector<AuxField*>&);

#endif
