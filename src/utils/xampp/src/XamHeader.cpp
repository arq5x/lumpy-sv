#include "XamHeader.h"

//#include "../utils/hts_memory.h"

using namespace std;


/**
 * @brief creates a SamHeader object that points to htslib memory already allocated
 *
 * @note the resulting SamHeader object shares ownership of the pre-allocated memory via
 *       shared_ptr reference counting
 */
XamHeader::XamHeader() :
  m_header(NULL)
{}

XamHeader::XamHeader(bam_hdr_t *header) :
  m_header(header)
{}

XamHeader::XamHeader(const XamHeader& other) :
    //deep copy?
    m_header (other.m_header)
{}

XamHeader& XamHeader::operator=(const XamHeader& other) 
{
    if ( &other == this )  
      return *this;
    //deep copy?
    m_header = other.m_header; ///< shared_ptr assignment will take care of deallocating old sam record if necessary
    return *this;
  }
/**
 * @brief Returns the length of the given sequence as stored in the \@SQ tag in the BAM header, or 0 if the sequence
 * name is not found.
 */
// uint32_t XamHeader::sequence_length(const std::string& sequence_name) const {
//   char *c = sequence_name.c_str();
//   for (int i = 0; i < m_header->n_targets; i++) {
// 	  if (strcmp(c,m_header->target_name[i]) == 0) {
// 		  return m_header->target_len[i];
// 	  }
//   }
//   return 0;
// }

/**
 * @brief extracts read group objects from a SAM header
 */
vector<ReadGroup> XamHeader::read_groups() const {
    const static auto RG_TAG = "@RG";
    const static auto NOT_FOUND = string::npos;
    auto result = vector<ReadGroup>();
    auto text = HeaderText();

    for (auto rg_start=text.find(RG_TAG), rg_end=rg_start; rg_start!=NOT_FOUND; rg_start=text.find(RG_TAG,rg_end+1) ) 
    {
        rg_end = text.find('\n', rg_start+1);
        auto rg_record = text.substr(rg_start, rg_end-rg_start);
        result.push_back(ReadGroup(rg_record));
    }
    return result;
}

bool XamHeader::IsCoordinateSorted() const {
    if (m_header->l_text > 3) 
    {
        if (strstr(m_header->text, "SO:coordinate") != NULL) 
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

RefVector XamHeader::GetRefData() const {
    RefVector rv;
    for (int i = 0; i < n_sequences(); ++i)
    {
        rv.push_back(RefData(sequence_name(i), sequence_length(i)));
    }
    return rv;
}

