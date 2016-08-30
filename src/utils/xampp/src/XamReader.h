#ifndef _XAMREADER_H
#define _XAMREADER_H

#include "Xam.h"
#include "XamHeader.h"
#include "sam.h" // from htslib
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <queue>
#include <functional>

using namespace std;

struct ChromOrder : public std::binary_function<Xam*, Xam*, bool>
{
    bool operator()(const Xam* lhs, const Xam* rhs) const
    {
        if (lhs->ChromId() == rhs->ChromId())
        {
            return !(lhs->Start() < rhs->Start());
        }
        else 
        {
            return !(lhs->ChromId() < rhs->ChromId());
        }
    }
};

class XamReader {

  public:

    bool Open(const std::string& filename);
    bool Open(const vector<string> & filenames);
    ~XamReader(void);

    string Chrom(uint32_t filenum, uint32_t chrom_id);
    string GetHeaderString(int filenum);
    XamHeader GetHeader(int filenum);
    RefVector GetRefData();
    bool IsCoordinateSorted();
    string GetFileName(int filenum);

    bool Next(Xam &rec);  // get the next record from the file(s)
    Xam *_rec;            // pointer to the record loaded by next(). user must delete


  private:
  	void _prime_buffer();
    vector<samFile*>   _file_ptrs;       ///< vector of pointers to the internal file structures of the sam/bam/cram files
    vector<bam_hdr_t*> _hdr_ptrs;        ///< vector of pointers to the internal header structures of the sam/bam/cram files
    vector<XamHeader*> _xamheader_ptrs;  ///< vector of pointers to XamHeader objects created for the sam/bam/cram files
    vector<string> _filenames;
    bool _has_been_read;
    int _num_files;
    priority_queue<Xam*, vector<Xam*>, ChromOrder> _buffer;
};

#endif