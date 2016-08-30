#ifndef _XAM_H
#define _XAM_H

//#include "ReadBases.h"
//#include "BaseQuals.h"
//#include "Cigar.h"
#include "XamTag.h"

#include "sam.h"
#include <string>
#include <memory>
#include <vector>

struct CigarOp {
  
    char     Type;   //!< CIGAR operation type (MIDNSHPX=)
    uint32_t Length; //!< CIGAR operation length (number of bases)
    
    //! constructor
    CigarOp(const char type = '\0', 
            const uint32_t& length = 0)
        : Type(type)
        , Length(length) 
    { }
};
const char cigar_ops_as_chars[] = { 'M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B' };

class XamReader;
class Xam {

 public:
    Xam(); 
    Xam(int filenum, XamReader *reader, bam1_t *body); 
    ~Xam(); 

    /**
    * @brief creates a deep copy of a sam record
    *
    * @note the copy will have exclusive ownership over the newly-allocated htslib memory
    *       until a data field (cigar, bases, etc.) is accessed, after which it will be
    *       shared via reference counting with the Cigar, etc. objects
    * @note does not perform a deep copy of the sam header; to copy the header,
    *       first get it via the header() function and then copy it via the usual C++
    *       semantics
    */
    Xam(const Xam& other);

    /**
    * @copydoc Sam::Xam(const Xam&)  
    */ 
    Xam& operator=(const Xam& other);

    /**
    * @brief moves Xam and its header accordingly. Shared pointers maintain state to all other associated objects correctly.
    */
    Xam(Xam&& other) = default;

    /**
    * @brief move assignment of a Xam and it's header. Shared pointers maintain state to all other associated objects correctly.
    */
    Xam& operator=(Xam&& other) = default;


    /**
    * @brief Filename from which the read came.
    *
    * @return the filename string.
    */
    std::string Filename() const;

    /**
    * @brief chromosome index of the read.
    *
    * Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). Similar to Picards getReferenceIndex()
    *
    * @return the integer representation of the chromosome.
    */
    uint32_t ChromId() const { return uint32_t(m_body->core.tid);     }       


    /**
    * @brief chromosome name of the read.
    *   *
    * @return the name of the chromosome.
    */
    std::string Chrom() const; 
    
    /**
    * @brief chromosome name of the mate.
    *   *
    * @return the name of the mate's chromosome.
    */
    std::string MateChrom() const;       

    /**
    * @brief the reference position of the first base in the read
    * @note the internal encoding is 0-based to mimic that of the BAM files.
    * @return a (1-based and inclusive) alignment start position (as you would see in a Xam file).
    */
    uint32_t Start() const { return uint32_t(m_body->core.pos);   }       

    /**
    * @brief returns a (1-based and inclusive) alignment stop position. 
    * @note the internal encoding is 0-based to mimic that of the BAM files. 
    * @note htslib's bam_endpos returns the coordinate of the first base AFTER the alignment, 0-based, so that translates into the last base IN the 1-based alignment.
    */
    uint32_t End() const { return uint32_t(bam_endpos(m_body)); }   

    /**
    * @brief calculates the theoretical alignment start of a read that has soft/hard-clips preceding the alignment
    *
    * For example if the read has an alignment start of 100 but the first 4 bases
    * were clipped (hard or soft clipped) then this method will return 96.
    *
    * @return the alignment start (1-based, inclusive) adjusted for clipped bases.  
    *
    * Invalid to call on an unmapped read.
    * Invalid to call with cigar = null
    */
    uint32_t UnclippedStart() const ;

    /**
    * @brief calculates the theoretical alignment stop of a read that has soft/hard-clips preceding the alignment
    *
    * For example if the read has an alignment stop of 100 but the last 4 bases
    * were clipped (hard or soft clipped) then this method will return 104.
    *
    * @return the alignment stop (1-based, inclusive) adjusted for clipped bases.     
    * @warning Invalid to call on an unmapped read.
    * @warning Invalid to call with cigar = null
    */
    uint32_t UnclippedStop() const ;

    /**
    * @brief returns the integer representation of the mate's chromosome. 
    *
    * Notice that chromosomes are listed in index order with regards to the header (so a 0-based number). 
    */
    uint32_t MateChromId() const { return uint32_t(m_body->core.mtid);    }       

    /**
    * @brief returns a (1-based and inclusive) mate's alignment start position (as you would see in a Xam file). 
    * @note the internal encoding is 0-based to mimic that of the BAM files.
    */
    uint32_t MateStart() const { return uint32_t(m_body->core.mpos);  }       

    void SetMateChromId(uint32_t chromid) {m_body->core.mtid = chromid;}
    void SetMateStart(uint32_t start) {m_body->core.mpos = start;}

    /**
    * @brief returns a (1-based and inclusive) mate's alignment stop position.
    * @note the internal encoding is 0-based to mimic that of the BAM files. 
    * @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
    */
    uint32_t MateStop() const ;                                                

    /**
    * @brief returns a (1-based and inclusive) mate's alignment stop position. 
    *
    * This overload is for usage when the user checks for the existence of the
    * tag themselves and passes it in to avoid exception throwing. This is provided
    * for performance conscious use of this function. This way you will only create 
    * one XamTag object for the mate cigar tag, instead of potentially two when checking
    * for its availability and then calling this function. For example: 
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_alignment_stop(tag) << endl;  // this will reuse the tag you have already obtained
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * This is better than the alternative using the other overload where you
    * have to either get the Tag twice or check for the exception thrown:
    *
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_alignment_stop() << endl;  // this will obtain a new tag internally (unnecessary)
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Xam. 
    * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
    * @note the internal encoding is 0-based to mimic that of the BAM files. 
    */
    //uint32_t mate_alignment_stop(const XamTag<std::string>& mate_cigar_tag) const ;       

    /**
    * @brief returns a (1-based and inclusive) mate's unclipped alignment start position. 
    * @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
    */
    uint32_t MateUnclippedStart() const;                                                 

    /**
    * @brief returns a (1-based and inclusive) mate's unclipped alignment start position.
    *
    * This overload is for usage when the user checks for the existence of the
    * tag themselves and passes it in to avoid exception throwing. This is provided
    * for performance conscious use of this function. This way you will only create 
    * one XamTag object for the mate cigar tag. Instead of potentially two when checking
    * for it's availability and then calling this function. For example: 
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_unclipped_start(tag) << endl;  // this will reuse the tag you have already obtained
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * This is better than the alternative using the other overload where you
    * have to either get the Tag twice or check for the exception thrown:
    *
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_unclipped_start() << endl;  // this will obtain a new tag internally (unnecessary)
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Xam. 
    * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
    * @note the internal encoding is 0-based to mimic that of the BAM files.
    */
    //uint32_t mate_unclipped_start(const XamTag<std::string>& mate_cigar_tag) const;        

    /**
    * @brief returns a (1-based and inclusive) mate's unclipped alignment stop position. @throw std::invalid_argument if called on a record that doesn't contain the mate cigar ("MC") tag.
    */
    uint32_t MateUnclippedStop() const ;                                                

    /**
    * @brief returns a (1-based and inclusive) mate's unclipped alignment stop position. 
    *
    * This overload is for usage when the user checks for the existence of the
    * tag themselves and passes it in to avoid exception throwing. This is provided
    * for performance conscious use of this function. This way you will only create 
    * one XamTag object for the mate cigar tag. Instead of potentially two when checking
    * for it's availability and then calling this function. For example: 
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_unclipped_stop(tag) << endl;  // this will reuse the tag you have already obtained
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * This is better than the alternative using the other overload where you
    * have to either get the Tag twice or check for the exception thrown:
    *
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    * const auto tag = record.string_tag("MC");  // obtains the tag from the record (expensive operation)
    * if (!missing(tag))
    *   cout << record.mate_unclipped_stop() << endl;  // this will obtain a new tag internally (unnecessary)
    * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    *
    * @param mate_cigar_tag the MC tag as obtained via the string_tag("MC") API in Xam. 
    * @warning This overload DOES NOT throw an exception if the mate cigar tag is missing. Instead it returns mate_alignment_start(). Treat it as undefined behavior. 
    * @note the internal encoding is 0-based to mimic that of the BAM files.
    */
    //uint32_t mate_unclipped_stop(const XamTag<std::string>& mate_cigar_tag) const;        

    /**
    * @brief returns the mapping quality of this alignment
    */
    uint8_t MappingQuality() const { return uint8_t(m_body->core.qual); }

    /**
    * @brief inferred insert size as reported by the aligner
    *
    * This is the signed observed insert size. If all segments are mapped to the Xame reference, 
    * the unsigned observed template length equals the number of bases from the leftmost mapped 
    * base to the rightmost mapped base. The leftmost segment has a plus sign and the rightmost 
    * has a minus sign. The sign of segments in the middle is undefined. 
    *
    * It is set as 0 for single read or when the information is unavailable.
    * 
    * @return a signed insert size or zero if it can't be inferred.
    */
    int32_t InsertSize() const { return m_body->core.isize; }

    // modify non-variable length fields (things outside of the data member)
    void set_chromosome(const uint32_t chr)              { m_body->core.tid  = int32_t(chr);        } ///< @brief simple setter for the chromosome index. Index is 0-based.
    void set_alignment_start(const uint32_t start)       { m_body->core.pos  = int32_t(start-1);    } ///< @brief simple setter for the alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
    void set_mate_chromosome(const uint32_t mchr)        { m_body->core.mtid = int32_t(mchr);       } ///< @brief simple setter for the mate's chromosome index. Index is 0-based.
    void set_mate_alignment_start(const uint32_t mstart) { m_body->core.mpos = int32_t(mstart - 1); } ///< @brief simple setter for the mate's alignment start. @warning You should use (1-based and inclusive) alignment but internally this is stored 0-based to simplify BAM conversion.
    void set_mapping_qual(const uint8_t mapq)            { m_body->core.qual = mapq;                } ///< @brief simple setter for the alignment quality
    void set_insert_size(const int32_t isize)            { m_body->core.isize = isize;              } ///< @brief simple setter for the insert size

    // getters for fields inside the data field
    std::string Name() const { return bam_get_qname(m_body); }                      ///< @brief returns the read name
    std::vector<CigarOp> CigarData() const;                                ///< @brief returns the cigar. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.
    std::string bases() const;                                                      ///< @brief returns the read bases. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.
    //BaseQuals base_quals() const { return BaseQuals{m_body}; }                    ///< @brief returns the base qualities. @warning the objects returned by this member function will share underlying htslib memory with this object. @warning creates an object but doesn't copy the underlying values.

    // getters for tagged values within the aux part of the data field
    XamTag<int32_t> GetIntegerTag(const std::string& tag_name) const;      ///< @brief retrieve an integer-valued tag by name. @warning creates an object but doesn't copy the underlying values.
    XamTag<double> double_tag(const std::string& tag_name) const;        ///< @brief retrieve an double/float-valued tag by name. @warning creates an object but doesn't copy the underlying values.
    XamTag<char> char_tag(const std::string& tag_name) const;            ///< @brief retrieve a char-valued tag by name. @warning creates an object but doesn't copy the underlying values.
    XamTag<std::string> GetStringTag(const std::string& tag_name) const;   ///< @brief retrieve a string-valued tag by name. @warning creates an object but doesn't copy the underlying values.
    bool HasTag(const std::string& tag_name) const;

    void AddStringTag(std::string label, char type, std::string value) const;

    // convenience getters for tag data
    std::string ReadGroup() const { return GetStringTag("RG").value(); };
    std::string FileName() const { return GetStringTag("LS").value(); };

    uint32_t edit_dist() const { return GetIntegerTag("NM").value(); };
    uint32_t alignment_score() const { return GetIntegerTag("AS").value(); };

    // getters for flags 
    bool Paired() const { return m_body->core.flag & BAM_FPAIRED;        }          ///< @brief whether or not this read is paired
    bool ProperlyPaired() const { return m_body->core.flag & BAM_FPROPER_PAIR;   } ///< @brief whether or not this read is properly paired (see definition in BAM spec)
    bool Unmapped() const { return m_body->core.flag & BAM_FUNMAP;         }        ///< @brief whether or not this read is unmapped
    bool MateUnmapped() const { return m_body->core.flag & BAM_FMUNMAP;        }   ///< @brief whether or not the mate read is unmapped
    bool Reverse() const { return m_body->core.flag & BAM_FREVERSE;       }         ///< @brief whether or not this read is from the reverse strand
    bool MateReverse() const { return m_body->core.flag & BAM_FMREVERSE;      }    ///< @brief whether or not the mate read is from the reverse strand
    bool First() const { return m_body->core.flag & BAM_FREAD1;         }           ///< @brief whether or not this read is the first read in a pair (or multiple pairs)
    bool Last() const { return m_body->core.flag & BAM_FREAD2;         }            ///< @brief whether or not this read is the last read in a pair (or multiple pairs)
    bool Secondary() const { return m_body->core.flag & BAM_FSECONDARY;     }       ///< @brief whether or not this read is a secondary alignment (see definition in BAM spec)
    bool Fail() const { return m_body->core.flag & BAM_FQCFAIL;        }            ///< @brief whether or not this read is marked as failing vendor (sequencer) quality control
    bool Duplicate() const { return m_body->core.flag & BAM_FDUP;           }       ///< @brief whether or not this read is a duplicate
    bool Supplementary() const { return m_body->core.flag & BAM_FSUPPLEMENTARY; }   ///< @brief whether or not this read is a supplementary alignment (see definition in the BAM spec)

    // modify flags
    void set_paired()            { m_body->core.flag |= BAM_FPAIRED;         } 
    void set_not_paired()        { m_body->core.flag &= ~BAM_FPAIRED;        }
    void set_unmapped()          { m_body->core.flag |= BAM_FUNMAP;          }
    void set_not_unmapped()      { m_body->core.flag &= ~BAM_FUNMAP;         }
    void set_mate_unmapped()     { m_body->core.flag |= BAM_FMUNMAP;         }
    void set_not_mate_unmapped() { m_body->core.flag &= ~BAM_FMUNMAP;        }
    void set_reverse()           { m_body->core.flag |= BAM_FREVERSE;        }
    void set_not_reverse()       { m_body->core.flag &= ~BAM_FREVERSE;       }
    void set_mate_reverse()      { m_body->core.flag |= BAM_FMREVERSE;       }
    void set_not_mate_reverse()  { m_body->core.flag &= ~BAM_FMREVERSE;      }
    void set_first()             { m_body->core.flag |= BAM_FREAD1;          }
    void set_not_first()         { m_body->core.flag &= ~BAM_FREAD1;         }
    void set_last()              { m_body->core.flag |= BAM_FREAD2;          }
    void set_not_last()          { m_body->core.flag &= ~BAM_FREAD2;         }
    void set_secondary()         { m_body->core.flag |= BAM_FSECONDARY;      }
    void set_not_secondary()     { m_body->core.flag &= ~BAM_FSECONDARY;     }
    void set_fail()              { m_body->core.flag |= BAM_FQCFAIL;         }
    void set_not_fail()          { m_body->core.flag &= ~BAM_FQCFAIL;        }
    void set_duplicate()         { m_body->core.flag |= BAM_FDUP;            }
    void set_not_duplicate()     { m_body->core.flag &= ~BAM_FDUP;           }
    void set_supplementary()     { m_body->core.flag |= BAM_FSUPPLEMENTARY;  }
    void set_not_supplementary() { m_body->core.flag &= ~BAM_FSUPPLEMENTARY; }

    bool empty() const { return m_body == nullptr; } ///< @brief whether or not this Xam object is empty, meaning that the internal memory has not been initialized (i.e. a Xam object initialized with Xam()).

    int filenum;         ///the input filenumber from which the record came
    
    bam1_t * GetHtsLibBody() const { return m_body; }
 private:
    XamReader *m_reader;
    bam1_t *m_body;      ///< htslib pointer to the Xam body structure
};

#endif
