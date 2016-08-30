#ifndef _XAMHEADER_H
#define _XAMHEADER_H

#include "sam.h"
#include "ReadGroup.h"

#include <memory>
#include <string>
#include <vector>

struct RefData {
   
    std::string RefName;    //!< name of reference sequence
    int32_t     RefLength;  //!< length of reference sequence
    
    //! constructor
    RefData(const std::string& name = "",
            const int32_t& length = 0)
        : RefName(name)
        , RefLength(length)
    { }
};

typedef std::vector<RefData> RefVector;

class XamHeader {
  public:
    XamHeader();                                            ///< @brief initializes a null SamHeader @warning if you need to create a SamHeader from scratch, use the builder instead
    XamHeader(bam_hdr_t *header);                                     ///< @brief creates a SamHeader given htslib object. @note used by all iterators
    XamHeader(const XamHeader& other);                                ///< @brief makes a deep copy of a SamHeader. Shared pointers maintain state to all other associated objects correctly.
    XamHeader& operator=(const XamHeader& other);                     ///< @brief deep copy assignment of a SamHeader. Shared pointers maintain state to all other associated objects correctly.
    XamHeader(XamHeader&& other) = default;                           ///< @brief moves SamHeader accordingly. Shared pointers maintain state to all other associated objects correctly.
    XamHeader& operator=(XamHeader&& other) = default;                ///< @brief move assignment of a SamHeader. Shared pointers maintain state to all other associated objects correctly.
    uint32_t n_sequences() const {return m_header->n_targets;}        ///< @brief Returns the number of reference sequences in the header
    uint32_t sequence_length(const std::string& sequence_name) const; ///< @brief Returns the length of the given reference sequence as stored in the \@SQ tag in the BAM header.
    uint32_t sequence_length(const uint32_t sequence_index) const { return m_header->target_len[sequence_index]; }                ///< @brief Returns the length of the given reference sequence as stored in the \@SQ tag in the BAM header.
    std::string sequence_name(const uint32_t sequence_index) const { return std::string(m_header->target_name[sequence_index]); } ///< @brief Returns the sequence name for the sequence with the given zero-based index
    std::vector<ReadGroup> read_groups() const;
    std::string HeaderText() const { return std::string(m_header->text, m_header->l_text); } ///< @brief Returns the text of the SAM header for parsing as an std::string. @note htslib only parses @SQ records, so gamgee is not simply a wrapper for C code.
    bool IsCoordinateSorted() const;
    RefVector GetRefData() const;
    bam_hdr_t * GetHtsLibHeader() const {return m_header; }
    private:
    bam_hdr_t *m_header;

  //friend class SamWriter;
  //friend class SamBuilder;
};

#endif