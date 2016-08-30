#ifndef _CIGAR_H
#define _CIGAR_H

#include "sam.h"

#include <memory>
#include <string>
#include <vector>
#include <sstream>

/**
 * @brief comprehensive list of valid cigar operators
 * @note order of the operators in this enum must match the order in BAM_CIGAR_STR from htslib/sam.h
 */
enum class CigarOperator { M, I, D, N, S, H, P, EQ, X, B };

using CigarElement = uint32_t;

/**
 * @brief Utility class to manage the memory of the cigar structure
 */
class Cigar {
 public:
  Cigar(const bam1_t *sam_record);
  Cigar(const Cigar& other);
  Cigar(Cigar&& other) = default;
  Cigar& operator=(const Cigar& other);
  Cigar& operator=(Cigar&& other) = default;
  ~Cigar() = default; ///< default destruction is sufficient, since our shared_ptr will handle deallocation

  CigarElement operator[](const uint32_t index) const;       ///< use freely as you would an array.
  CigarElement& operator[](const uint32_t index);            ///< use freely as you would an array
  uint32_t size() const { return m_num_cigar_elements; } ///< number of base qualities in the container
  bool operator==(const Cigar& other) const;  ///< check for equality with another Cigar
  bool operator!=(const Cigar& other) const;  ///< check for inequality with another Cigar
  std::string to_string() const;  ///< produce a string representation of this Cigar

  /**
   * @brief gets the operator of an individual cigar element
   */
  inline static CigarOperator cigar_op(const CigarElement cigar_element) {
    return static_cast<CigarOperator>(bam_cigar_op(cigar_element));
  }

  /**
   * @brief gets the length of an individual cigar element
   */
  inline static uint32_t cigar_oplen(const CigarElement cigar_element) {
    return bam_cigar_oplen(cigar_element);
  }

  /**
   * @brief creates an encoded htslib cigar element suitable for direct insertion into a Cigar
   *        out of a length and a CigarOperator
   */
  inline static CigarElement make_cigar_element(const uint32_t oplen, const CigarOperator op) {
    return (oplen << BAM_CIGAR_SHIFT) | static_cast<uint32_t>(op);
  }
  
  /**
   * @brief Table used to parse chars representing cigar operations into their htslib encodings
   *
   * @note: This table should eventually be moved to htslib. Currently htslib dynamically allocates
   *        and fills the equivalent of this table on-demand, which makes little sense.
   */
  static const std::vector<int8_t> cigar_op_parse_table;

  /**
   * @brief utility function to parse cigar strings (in a stringstream) one element at a time
   *
   * This function is useful if you are parsing a cigar string into elements and operators. A simple way to iterate over all the elements is :
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * auto cigar_stream = stringstream{cigar_string};  // requires gcc 4.9's stdlibc++ if you can't use that, don't use auto on this line
   * while (cigar_stream.peek() != std::char_traits<char>::eof()) {
   *   const auto element = parse_next_cigar_element(cigar_stream);
   *   do_something_with(Cigar::cigar_op(element), Cigar::cigar_oplen(element));
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * @return a CigarElement with length and operator that can be extracted using the Cigar API (Cigar::cigar_op and Cigar::cigar_oplen) to actual elements.
   */
  static CigarElement parse_next_cigar_element (std::stringstream& cigar_stream);

  inline static bool consumes_read_bases(const CigarOperator op) {return bam_cigar_type(static_cast<int8_t>(op))&1;}      ///< returns true if operator is one of the following: Match (M), Insertion (I), Soft-Clip (S), Equal (=) or Different (X)
  inline static bool consumes_reference_bases(const CigarOperator op) {return bam_cigar_type(static_cast<int8_t>(op))&2;} ///< returns true if operator is one of the following: Match (M), Deletion (D), Reference-Skip (N), Equal (=) or Different (X)



 private:
  bam1_t *m_xam_record;                   ///< sam record containing our cigar, potentially co-owned by multiple other objects
  uint32_t *m_cigar;                      ///< pointer to the start of the cigar in m_sam_record, cached for efficiency
  uint32_t m_num_cigar_elements;          ///< number of elements in our cigar

  static const char cigar_ops_as_chars[]; ///< static lookup table to convert CigarOperator enum values to chars.

};

#endif
