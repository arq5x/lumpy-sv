#ifndef _READGROUP_H
#define _READGROUP_H

#include <string>
#include "split.h"
/**
 * @brief Helper struct to hold one read group record from a sam file header
 *
 * A read group, which is contained in a line starting with "@RG" in a sam file header,
 * contains information about the seuencing process: type of machine, time, center etc.
 * Most important for variant calling and quality control is the sample ID -- each sample
 * may be split into several read groups.
 */
class ReadGroup {

public:
  std::string id;
  std::string center;
  std::string description;
  std::string date_time;
  std::string flow_order;
  std::string key_sequence;
  std::string library;
  std::string programs;
  std::string median_insert_size;
  std::string platform;
  std::string platform_unit;
  std::string sample;

  ReadGroup() = default;
  ReadGroup(const std::string& header_line);

private:
  //static const char ID_TAG [];
  //static const char CENTER_TAG [];
  //static const char DESCRIPTION_TAG [];
  //static const char DATE_TIME_TAG [];
  //static const char FLOW_ORDER_TAG [];
  //static const char KEY_SEQUENCE_TAG [];
  //static const char LIBRARY_TAG [];
  //static const char PROGRAMS_TAG [];
  //static const char MEDIAN_INSERT_SIZE_TAG [];
  //static const char PLATFORM_TAG [];
  //static const char PLATFORM_UNIT_TAG [];
  //static const char SAMPLE_TAG [];

  static bool starts_with(const std::string& token, const char tag[2]) {
      return token[0] == tag[0] && token[1] == tag[1];
  }
};

#endif
