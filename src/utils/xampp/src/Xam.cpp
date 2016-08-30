#include "Xam.h"
#include "XamReader.h"
//#include "cigar.h"
#include <iostream>
#include <string>

using namespace std;

const string MATE_CIGAR_TAG = "MC";

Xam::Xam()
:
	filenum(-1),
	m_reader(NULL),
	m_body(NULL)
{}

Xam::Xam(int filenum, XamReader *reader, bam1_t *body)
:
	filenum(filenum),
	m_reader(reader),
	m_body(body)
{}

Xam::~Xam()
{
    // use htslib's free
    bam_destroy1(m_body);
}

Xam::Xam(const Xam& other) :
  m_reader ( other.m_reader ),
  // deep copy here instead?
  m_body ( other.m_body )
{}

Xam& Xam::operator=(const Xam& other) {
  if ( &other == this )  
    return *this;
  filenum = other.filenum;
  m_reader = other.m_reader;      ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_body = other.m_body;     ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  return *this;
}

XamTag<int32_t> Xam::GetIntegerTag(const std::string& tag_name) const {
	const uint8_t *aux_ptr = bam_aux_get(m_body, tag_name.c_str());
	return aux_ptr == nullptr ? XamTag<int32_t>(tag_name, 0, true) :
	                          XamTag<int32_t>(tag_name, bam_aux2i(aux_ptr));

  // TODO: htslib bam_aux2i() returns 0 if the tag is not of integer type,
  //       but 0 is a valid value for a correctly-typed tag. This should be patched.
}

XamTag<double> Xam::double_tag(const std::string& tag_name) const {
	const uint8_t *aux_ptr = bam_aux_get(m_body, tag_name.c_str());
	return aux_ptr == nullptr ? XamTag<double>(tag_name, 0.0, true) :
	                          XamTag<double>(tag_name, bam_aux2f(aux_ptr));

  // TODO: htslib bam_aux2f() returns 0.0 if the tag is not of double/float type,
  //       but 0.0 is a valid value for a correctly-typed tag. This should be patched.
}

XamTag<char> Xam::char_tag(const std::string& tag_name) const {
	const uint8_t *aux_ptr = bam_aux_get(m_body, tag_name.c_str());
	if ( aux_ptr == nullptr )  // tag doesn't exist
		return XamTag<char>(tag_name, '\0', true);

	const char char_val = bam_aux2A(aux_ptr);
	if ( char_val == '\0' )    // tag not of char type
		return XamTag<char>(tag_name, '\0', true);

	return XamTag<char>(tag_name, char_val);
}

XamTag<std::string> Xam::GetStringTag(const std::string& tag_name) const 
{
	const uint8_t *aux_ptr = bam_aux_get(m_body, tag_name.c_str());
	if ( aux_ptr == nullptr )  // tag doesn't exist
		return XamTag<string>(tag_name, "", true);

	const char *str_ptr = bam_aux2Z(aux_ptr);
	if ( str_ptr == nullptr )  // tag not of string type
		return XamTag<string>(tag_name, "", true);

	return XamTag<string>(tag_name, string(str_ptr));
}

bool Xam::HasTag(const std::string& tag_name) const
{
	const uint8_t *aux_ptr = bam_aux_get(m_body, tag_name.c_str());
	if ( aux_ptr == nullptr )  // tag doesn't exist
		return false;
	return true;
}

void Xam::AddStringTag(string label, char type, string value) const
{
	bam_aux_append(m_body, label.c_str(), type, value.size() + 1, (uint8_t*)value.c_str());
}

string Xam::Chrom() const
{
	return m_reader->Chrom(filenum, m_body->core.tid);
}

string Xam::MateChrom() const
{ 
	return m_reader->Chrom(filenum, m_body->core.mtid);
}

string Xam::Filename() const 
{ 
	return m_reader->GetFileName(filenum); 
}

vector<CigarOp> Xam::CigarData() const
{
	uint32_t *cigar = bam_get_cigar(m_body);  
	uint32_t num_cigar_elements = m_body->core.n_cigar;
	vector<CigarOp> cd;

	for ( uint32_t i = 0; i < num_cigar_elements; ++i )
	{
        char type = cigar_ops_as_chars[static_cast<int>(bam_cigar_op(cigar[i]))];
        uint32_t len = bam_cigar_oplen(cigar[i]);
		cd.push_back(CigarOp(type, len));
	}
	return cd;
}

string Xam::bases() const 
{ 
	uint8_t* bases = bam_get_seq(m_body);
    uint32_t num_bases = m_body->core.l_qseq;
	stringstream ss;
    for ( uint32_t i = 0; i < num_bases; i++ ) 
    {
    	switch(bam_seqi(bases, i))
    	{
	    	case 1  :
	       		ss << "A"; break;
	       	case 2  :
	       		ss << "C"; break;
	       	case 4  :
	       		ss << "G"; break;
	       	case 8  :
	       		ss << "T"; break;
	       	case 15  :
	       		ss << "N"; break;
	    }
    }
    return ss.str(); 
}


