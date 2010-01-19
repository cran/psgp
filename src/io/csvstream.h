#ifndef CSVSTREAM_H_
#define CSVSTREAM_H_

// Uncomment following line to enable debugging messages
// #define DEBUG
#include <iostream>
#include <string>
#include <vector>
#include <itpp/itbase.h>
#include <cassert>

#define debug_msg(msg)

using namespace std;
using namespace itpp;

/**
 * This class provides support for reading from and
 * writing to CSV (comma separated values) files. At the moment,
 * it only supports reading to and writing from matrices, and does
 * not (yet) work like a proper C++ stream, although support for 
 * streamed input/output will be added later.
 */
class csvstream
{
public:
	csvstream();
	virtual ~csvstream();
	
	int read(mat &matrix, const string filename);
	int write(const mat matrix, const string filename, int decimals = 5);

};

#endif /*CSVSTREAM_H_*/
