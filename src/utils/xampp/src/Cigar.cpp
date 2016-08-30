#include "Cigar.h"
#include <string>
#include <sstream>
#include <stdexcept>

using namespace std;

/**
  * @warning the order of these operators are matching the defines in htslib, do not change!
  */
const char Cigar::cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };
const std::vector<int8_t> Cigar::cigar_op_parse_table = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,7,-1,-1,-1,-1,9,-1,2,-1,-1,-1,5,1,-1,-1,-1,0,3,-1,6,-1,-1,4,-1,-1,-1,-1,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

/**
  * @brief creates a Cigar object that points to htslib memory already allocated
  *
  * @note the resulting Cigar object shares ownership of the pre-allocated memory via
  *       shared_ptr reference counting
  */
Cigar::Cigar(const std::shared_ptr<bam1_t>& sam_record) :
  m_sam_record { sam_record },
  m_cigar { bam_get_cigar(sam_record.get()) },
  m_num_cigar_elements { (sam_record.get())->core.n_cigar }
{}

/**
  * @brief creates a deep copy of a Cigar object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
Cigar::Cigar(const Cigar& other) :
  m_sam_record { utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())) },
  m_cigar { bam_get_cigar(m_sam_record.get()) },
  m_num_cigar_elements { other.m_num_cigar_elements }
{}

/**
  * @brief creates a deep copy of a Cigar object
  *
  * @note the copy will have exclusive ownership over the newly-allocated htslib memory
  */
Cigar& Cigar::operator=(const Cigar& other) {
  if (&other == this) ///< check for self assignment
    return *this; 
  m_sam_record = utils::make_shared_sam(utils::sam_deep_copy(other.m_sam_record.get())); ///< shared_ptr assignment will take care of deallocating old sam record if necessary
  m_cigar = bam_get_cigar(m_sam_record.get());
  m_num_cigar_elements = other.m_num_cigar_elements;
  return *this;
}

/**
  * @brief access an individual cigar element by index
  *
  * @return cigar element at the specified index as an encoded uint32_t. use cigar_op() and
  *         cigar_oplen() to unpack the cigar operator and length
  */
CigarElement Cigar::operator[](const uint32_t index) const {
  utils::check_max_boundary(index, m_num_cigar_elements);
  return m_cigar[index];
}

/**
  * @brief access and/or modify an individual cigar element by index
  *
  * @return cigar element at the specified index as an encoded uint32_t. use cigar_op() and
  *         cigar_oplen() to unpack the cigar operator and length
  */
CigarElement& Cigar::operator[](const uint32_t index) {
  utils::check_max_boundary(index, m_num_cigar_elements);
  return m_cigar[index];
}

bool Cigar::operator==(const Cigar& other) const {
  if ( m_num_cigar_elements != other.m_num_cigar_elements )
    return false;

  for ( auto i = 0u; i < m_num_cigar_elements; ++i ) {
    if ( m_cigar[i] != other.m_cigar[i] )
      return false;
  }

  return true;
}

bool Cigar::operator!=(const Cigar& other) const {
  return !(*this == other);
}


/**
  * @brief returns a string representation of this cigar
  */
string Cigar::to_string() const {
  stringstream stream;

  for ( uint32_t i = 0; i < m_num_cigar_elements; ++i ) 
    stream << cigar_oplen(m_cigar[i]) << cigar_ops_as_chars[static_cast<int>(cigar_op(m_cigar[i]))];
  return stream.str();
}

CigarElement Cigar::parse_next_cigar_element (stringstream& cigar_stream) {
  auto element_length = uint32_t{0};
  auto element_op = uint8_t{0};
  cigar_stream >> element_length;
  if ( cigar_stream.fail() )
    throw invalid_argument(string("Error parsing cigar string: ") + cigar_stream.str());
  cigar_stream >> element_op;
  if ( cigar_stream.fail() || int(element_op) >= 128 )
    throw invalid_argument(string("Error parsing cigar string: ") + cigar_stream.str());
  const auto encoded_op = cigar_op_parse_table[element_op];
  if ( encoded_op < 0 )
    throw invalid_argument(string("Unrecognized operator ") + char(element_op) + " in cigar string: " + cigar_stream.str());
  return make_cigar_element(element_length, static_cast<CigarOperator>(encoded_op));
} 