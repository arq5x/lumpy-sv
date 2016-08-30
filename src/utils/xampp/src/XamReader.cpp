#include "XamReader.h"

XamReader::~XamReader()
{}

bool XamReader::Open(const string &filename)
{
	_filenames.push_back(filename);
	_file_ptrs.push_back(sam_open(filename.empty() ? "-" : filename.c_str(), "r"));
	if ( _file_ptrs[0] == nullptr ) {
		fprintf(stderr, "Couldn't open \"%s\"\n", filename.c_str());
		return false;
	}
	else {
		bam_hdr_t *_hdr_ptr = sam_hdr_read(_file_ptrs[0]);
    	if (_hdr_ptr == NULL) {
	        fprintf(stderr, "Couldn't read header for \"%s\"\n", filename.c_str());
	        return false;
    	}
    	else 
    	{
    		_hdr_ptrs.push_back(_hdr_ptr);
    		_xamheader_ptrs.push_back(new XamHeader(_hdr_ptr));
    	}
	}
	_num_files = 1;
	_prime_buffer();
	return true;
}

bool XamReader::Open(const vector<string> &filenames)
{
	_filenames = filenames;

	vector<string>::const_iterator fIt  = filenames.begin();
    vector<string>::const_iterator fEnd = filenames.end();
    for (; fIt != fEnd; ++fIt) {
    	samFile *_file_ptr = sam_open(fIt->empty() ? "-" : (*fIt).c_str(), "r");
    	if ( _file_ptr == nullptr ) {
			fprintf(stderr, "Couldn't open \"%s\"\n", (*fIt).c_str());
			return false;
		}
		else {
			_file_ptrs.push_back(_file_ptr);
			bam_hdr_t *_hdr_ptr = sam_hdr_read(_file_ptr);
    		if (_hdr_ptr == NULL) {
	        	fprintf(stderr, "Couldn't read header for \"%s\"\n", (*fIt).c_str());
	        	return false;
    		}
    		else {
    			_hdr_ptrs.push_back(_hdr_ptr);
    			_xamheader_ptrs.push_back(new XamHeader(_hdr_ptr));
    			_num_files++;
    		}
		}
    }
    _prime_buffer();
    return true;
}

void XamReader::_prime_buffer()
{
	for (int i = 0; i < _num_files; ++i)
	{
		bam1_t *rec_ptr = bam_init1();
		if (sam_read1(_file_ptrs[i], _hdr_ptrs[i], rec_ptr) >= 0)
		{
			Xam *r = new Xam(i, this, rec_ptr);
			_buffer.push(r);
		}
		else 
		{
			fprintf(stderr, "Couldn't read record for file %s\"\n", _filenames[i].c_str());
		}
	}
}

bool XamReader::Next(Xam &rec)
{
	if (_buffer.empty())
	{
	 	return false;
	}
	_rec = _buffer.top();
	rec = *_rec;
	_buffer.pop();

	//// read next record from file that was popped
	bam1_t *rec_ptr = bam_init1();
	if (sam_read1(_file_ptrs[_rec->filenum], _hdr_ptrs[_rec->filenum], rec_ptr) >= 0)
	{
		Xam *r = new Xam(_rec->filenum, this, rec_ptr);
		_buffer.push(r);
	}
	return true;
}

string XamReader::Chrom(uint32_t filenum, uint32_t chrom_id)
{
	return _xamheader_ptrs[filenum]->sequence_name(chrom_id);
}

string XamReader::GetHeaderString(int filenum)
{
	return _xamheader_ptrs[filenum]->HeaderText();
}

XamHeader XamReader::GetHeader(int filenum)
{
	return *_xamheader_ptrs[filenum];
}

RefVector XamReader::GetRefData()
{
	return _xamheader_ptrs[0]->GetRefData();	
}

string XamReader::GetFileName(int filenum)
{
	return _filenames[filenum];
}

bool XamReader::IsCoordinateSorted()
{
	for (int i = 0; i < _num_files; ++i)
	{
		if (!_xamheader_ptrs[i]->IsCoordinateSorted())
		{
			return false;
		}
	}
	return true;
}

