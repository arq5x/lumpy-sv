#ifndef _XAMWRITER_H
#define _XAMWRITER_H

#include <string>
#include <memory>

#include "Xam.h"
#include "XamHeader.h"

#include "sam.h"


/**
 * @brief utility class to write out a SAM/BAM/CRAM file to any stream
 * @todo add serialization option
 */
class XamWriter {

 public: 

  /**
   * @brief Creates a new SamWriter using the specified output file name
   * @param output_fname file to write to. The default is stdout (as defined by htslib)
   * @param binary whether the output should be in BAM (true) or SAM format (false) 
   * @note the header is copied and managed internally
   */
  explicit XamWriter(const std::string& output_fname = "-", const bool binary = true);

  /**
   * @brief Creates a new SamWriter with the header extracted from a Sam record and using the specified output file name
   * @param header       SamHeader object to make a copy from
   * @param output_fname file to write to. The default is stdout  (as defined by htslib)
   * @param binary whether the output should be in BAM (true) or SAM format (false) 
   * @note the header is copied and managed internally
   */
  explicit XamWriter(const XamHeader& header, const std::string& output_fname = "-", const bool binary = true);

  /**
   * @brief a SamWriter cannot be copied safely, as it is iterating over a stream.
   */

  XamWriter(const XamWriter& other) = delete;
  XamWriter& operator=(const XamWriter& other) = delete;

  /**
   * @brief a SamWriter can be moved
   */

  XamWriter(XamWriter&& other) = default;
  XamWriter& operator=(XamWriter&& other) = default;

  /**
   * @brief Adds a record to the file stream
   * @param body the record
   */
  void AddRecord(const Xam& xam);

  /**
   * @brief Adds a header to the file stream.
   * @param header the header
   * @note the header is a requirement to add records
   */
  void AddHeader(const XamHeader& header);

  static htsFile* Open(const std::string& output_fname, const std::string& binary);
  bool Open(const XamHeader& header, const std::string& output_fname, const bool binary);
  bool Close();

 private:
  htsFile *m_out_file;  ///< the file or stream to write out to ("-" means stdout)
  XamHeader m_header;  ///< holds a copy of the header throughout the production of the output (necessary for every record that gets added)

  void WriteHeader() const;

};

#endif