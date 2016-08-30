#include "XamWriter.h"


XamWriter::XamWriter(const std::string& output_fname, const bool binary) :
  m_out_file (Open(output_fname, binary ? "wb" : "w")),
  m_header (nullptr)
{}

XamWriter::XamWriter(const XamHeader& header, const std::string& output_fname, const bool binary) :
  m_out_file (Open(output_fname, binary ? "wb" : "w")),
  m_header (header)
{
  WriteHeader();
}

void XamWriter::WriteHeader() const {
  sam_hdr_write(m_out_file, m_header.GetHtsLibHeader());
}


void XamWriter::AddHeader(const XamHeader& header) { 
  m_header = header;
  WriteHeader();
}

void XamWriter::AddRecord(const Xam& xam) { 
  sam_write1(m_out_file, m_header.GetHtsLibHeader(), xam.GetHtsLibBody());
}

htsFile* XamWriter::Open(const std::string& output_fname, const std::string& mode) {
  return sam_open(output_fname.empty() ? "-" : output_fname.c_str(), mode.c_str());
}

bool XamWriter::Open(const XamHeader& header, const std::string& output_fname, const bool binary) 
{
    m_out_file = Open(output_fname, binary ? "wb" : "w");
    m_header = header;
    WriteHeader();
    return true;
}

bool XamWriter::Close() {
  return sam_close(m_out_file);
}
