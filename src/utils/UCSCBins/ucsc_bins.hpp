#ifndef __UCSC_BINS_HPP__
#define __UCSC_BINS_HPP__
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <stdint.h> 


using namespace std;

typedef uint64_t CHR_POS;
typedef uint16_t BIN_LEVEL;
typedef uint32_t BIN;
typedef uint16_t USHORT;
typedef uint32_t UINT;

//{{{ template <class T> class UCSCElement
template <class T>
class UCSCElement
{
	public:
		unsigned long int id;
		CHR_POS start, end;
		char strand;
		string chr;
		T value;

		int overlap(UCSCElement<T> e);
		
		static bool sort_ucscelement_by_start(UCSCElement<T> e,
											  UCSCElement<T> f);

		static bool sort_ucscelement_by_value(UCSCElement<T> e,
											  UCSCElement<T> f);

		static bool compare_ucscelement_by_value(UCSCElement<T> e,
												 UCSCElement<T> f);
		static int overlap(CHR_POS as,
						   CHR_POS ae,
						   CHR_POS bs,
						   CHR_POS be);

		UCSCElement();
		UCSCElement(CHR_POS start_, CHR_POS end_);
		UCSCElement(CHR_POS start_, CHR_POS end_, char strand_);
		UCSCElement(CHR_POS start_, CHR_POS end_, char strand_, string chr_);
};
//}}}

//{{{ template <class T> class UCSCBins 
template <class T>
class UCSCBins 
{
	unsigned long int next_id;
	unsigned int size;

	typedef vector<UCSCElement<T> > element_vector;
	typedef map<BIN, element_vector, less<BIN> > bins;
	map<string, bins, less<string> > chrom_bins;

	public:
		const static BIN kNumBins;
		const static BIN_LEVEL kBinLevels;
		const static BIN kBinOffsets[];

		const static USHORT kBinFirstShift;
		const static USHORT kBinNextShift;

		static BIN getBin(CHR_POS start,
						  CHR_POS end);

		unsigned int num_bps();

		void add(string chr,
				 CHR_POS start,
				 CHR_POS end,
				 UCSCElement<T> element);

		void add(string chr,
				 CHR_POS start,
				 CHR_POS end,
				 char strand,
				 T value);

		int remove(UCSCElement<T>,
					bool remove_value,
					bool test_strand,
					bool test_id);

		vector< UCSCElement<T> > get(string chr,
									 CHR_POS start,
									 CHR_POS end,
									 char strand,
									 bool test_strand);

		vector< UCSCElement<T> >values();
		vector< UCSCElement<T> >values(string chr);
		vector< UCSCElement<T> >values(string chr,CHR_POS max_pos);

		UCSCBins();

};
//}}}

//{{{ template <class T> UCSCElement<T>::UCSCElement() { ; }
template <class T>
UCSCElement<T>::UCSCElement() { ; }
//}}}

//{{{ template <class T> UCSCElement<T>::UCSCElement(CHR_POS start_, CHR_POS
template <class T>
UCSCElement<T>::UCSCElement(CHR_POS start_, CHR_POS end_)
{
	start = start_;
	end = end_;
}
//}}}

//{{{ template <class T> UCSCElement<T>::UCSCElement(CHR_POS start_, CHR_POS
template <class T>
UCSCElement<T>::UCSCElement(CHR_POS start_, CHR_POS end_, char strand_)
{
	start = start_;
	end = end_;
	strand = strand_;
}
//}}}

//{{{ template <class T> UCSCElement<T>::UCSCElement(CHR_POS start_,
template <class T>
UCSCElement<T>::UCSCElement(CHR_POS start_,
							CHR_POS end_,
							char strand_,
							string chr_)
{
	start = start_;
	end = end_;
	strand = strand_;
	chr = chr_;
}
//}}}

//{{{ template <class T> int UCSCElement<T>::overlap(UCSCElement<T> e)
template <class T>
int
UCSCElement<T>::overlap(UCSCElement<T> e)
{
	return min(end, e.end) - max(start, e.start) + 1;
}      
//}}}

//{{{ template <class T> int UCSCElement<T>::overlap(CHR_POS as,
template <class T>
int
UCSCElement<T>::overlap(CHR_POS as,
						CHR_POS ae,
						CHR_POS bs,
						CHR_POS be)
{
	return min(ae, be) - max(as, bs) + 1;
}      
//}}}

/* Genome binning constants */
//{{{ 
template <class T>
const BIN UCSCBins<T>::kNumBins   = 37450;

template <class T>
const BIN_LEVEL UCSCBins<T>::kBinLevels = 7;

// bins range in size from 16kb to 512Mb
// Bin  0          spans 512Mbp,   # Level 1
// Bins 1-8        span 64Mbp,     # Level 2
// Bins 9-72       span 8Mbp,      # Level 3
// Bins 73-584     span 1Mbp       # Level 4
// Bins 585-4680   span 128Kbp     # Level 5
// Bins 4681-37449 span 16Kbp      # Level 6
template <class T>
const BIN UCSCBins<T>::kBinOffsets[] = {32678+4096+512+64+8+1,
										4096+512+64+8+1,
										512+64+8+1,
										64+8+1,
										8+1,
										1,
										0};

/* How much to shift to get to finest bin. */
template <class T>
const USHORT UCSCBins<T>::kBinFirstShift = 14;

/* How much to shift to get to next larger bin. */
template <class T>
const USHORT UCSCBins<T>::kBinNextShift = 3;
//}}}

//{{{ template <class T> UCSCBins<T>::UCSCBins()
template <class T>
UCSCBins<T>::UCSCBins()
{
	size = 0;
	next_id = 0;
	chrom_bins = map<string, bins, less<string> >();
}
//}}}

//{{{ template <class T> void UCSCBins<T>::add(string chr,
template <class T>
void
UCSCBins<T>::add(string chr,
			  CHR_POS start,
			  CHR_POS end,
			  UCSCElement<T> element) 
{
	++size;
	element.id = next_id++;
	BIN bin = getBin(start, end);
	chrom_bins[chr][bin].push_back(element);
}
//}}}

//{{{ template <class T> void UCSCBins<T>:: add(string chr,
template <class T>
void
UCSCBins<T>:: add(string chr,
				  CHR_POS start,
				  CHR_POS end,
				  char strand,
				  T value)
{
	UCSCElement<T> element = UCSCElement<T>(start, end, strand, chr);
	element.value = value;
	element.id = next_id++;
	++size;
	BIN bin = getBin(start, end);
	chrom_bins[chr][bin].push_back(element);
}
//}}}

//{{{ template <class T> BIN UCSCBins<T>::getBin(CHR_POS start,
template <class T>
BIN
UCSCBins<T>::getBin(CHR_POS start,
					CHR_POS end)
{
	--end;
	start >>= kBinFirstShift;
	end   >>= kBinFirstShift;

	for (register short i = 0; i < kBinLevels; ++i) {
		if (start == end)
			return kBinOffsets[i] + start;
		start >>= kBinNextShift;
		end   >>= kBinNextShift;
	}

	return 0;
}
//}}}

//{{{ template <class T> vector<UCSCElement<T> > UCSCBins<T>::get(string chr,
template <class T>
vector<UCSCElement<T> >
UCSCBins<T>::get(string chr,
				 CHR_POS start,
				 CHR_POS end,
				 char strand,
				 bool test_strand)
{
	vector<UCSCElement<T> > hits;

	BIN start_bin, end_bin;
	start_bin = (start >> kBinFirstShift);
	end_bin = ( (end - 1) >> kBinFirstShift);

	for (BIN_LEVEL i = 0; i < kBinLevels; ++i) {
		BIN offset = kBinOffsets[i];
		for (BIN j = start_bin + offset; j <= end_bin + offset; ++j) {

			typename vector< UCSCElement<T> >::const_iterator elementItr;

			for (	elementItr = chrom_bins[chr][j].begin();
					elementItr != chrom_bins[chr][j].end();
					++elementItr) {
				
				if ( UCSCElement<T>::overlap(start,
											 end,
											 elementItr->start,
											 elementItr->end) > 0 ) {
					bool strands_are_same = true;

					if ( test_strand )
						strands_are_same = (strand == elementItr->strand);

					if (strands_are_same)
						hits.push_back(*elementItr);
				}
			}
		}
		start_bin >>= kBinNextShift;
		end_bin   >>= kBinNextShift;
	}

	return hits;
}
//}}}

//{{{ template <class T> int UCSCBins<T>::remove(UCSCElement<T> e,
template <class T>
int
UCSCBins<T>::remove(UCSCElement<T> e,
					bool remove_value,
					bool test_strand,
					bool test_id)
{
	
	BIN start_bin, end_bin;
	start_bin = (e.start >> kBinFirstShift);
	end_bin = ( (e.end - 1) >> kBinFirstShift);

	for (BIN_LEVEL i = 0; i < kBinLevels; ++i) {
		BIN offset = kBinOffsets[i];
		for (BIN j = start_bin + offset; j <= end_bin + offset; ++j) {

			typename vector< UCSCElement<T> >::iterator elementItr;

			for (	elementItr = chrom_bins[e.chr][j].begin();
					elementItr != chrom_bins[e.chr][j].end();
					++elementItr) {
				if ( UCSCElement<T>::overlap(e.start,
											 e.end,
											 elementItr->start,
											 elementItr->end) > 0 ) {
					bool strands_are_same = true;
					bool ids_are_same = true;

					if ( test_strand )
						strands_are_same = (e.strand == elementItr->strand);

					if ( test_id )
						ids_are_same = (e.id == elementItr->id);

					if ( strands_are_same && ids_are_same )  {
						if (remove_value)
							free(elementItr->value);
						chrom_bins[e.chr][j].erase(elementItr);
						--size;
						return 0;
					} 
				}
			}
		}
		start_bin >>= kBinNextShift;
		end_bin   >>= kBinNextShift;
	}
	return 1;
}
//}}}

//{{{ template <class T> unsigned int UCSCBins<T>::num_bps()
template <class T>
unsigned int
UCSCBins<T>::num_bps()
{
	return size;
}
//}}}

//{{{ template <class T> vector<UCSCElement<T> > UCSCBins<T>::values(string
template <class T>
vector<UCSCElement<T> >
UCSCBins<T>::values(string target_chr)
{

	//cerr << "values:" << target_chr << "\t";
	vector< UCSCElement<T> > values;

	typename map<BIN, element_vector >::iterator bin_it;

	for (	bin_it = chrom_bins[target_chr].begin();
			bin_it != chrom_bins[target_chr].end();
			++bin_it ) 
		values.insert( values.end(),
					   (*bin_it).second.begin(),
					   (*bin_it).second.end() );

	return values;
}
//}}}

//{{{ template <class T> vector<UCSCElement<T> > UCSCBins<T>::values(string
template <class T>
vector<UCSCElement<T> >
UCSCBins<T>::values(string target_chr,
					CHR_POS max_pos)
{
	vector< UCSCElement<T> > values;

	typename map<BIN, element_vector >::iterator bin_it;

	for (	bin_it = chrom_bins[target_chr].begin();
			bin_it != chrom_bins[target_chr].end();
			++bin_it ) {

		typename vector< UCSCElement<T> >::iterator val_it;

		for (	val_it =    (*bin_it).second.begin();
				val_it !=   (*bin_it).second.end() ;
				++val_it ) 
			if (val_it->start < max_pos) 
				values.push_back( *val_it );

	}
	
	return values;
}
//}}}

//{{{ template <class T> vector<UCSCElement<T> > UCSCBins<T>::values()
template <class T>
vector<UCSCElement<T> >
UCSCBins<T>::values()
{
	vector< UCSCElement<T> > values;
	
	/*
	 * typedef vector<UCSCElement<T> > element_vector;
	 * typedef map<BIN, element_vector, less<BIN> > bins;
	 * map<string, bins, less<string> > chrom_bins;
	 */

	typename map<string,bins >::iterator chr_it;

	for (chr_it = chrom_bins.begin(); chr_it != chrom_bins.end(); ++chr_it) {

		typename map<BIN, element_vector >::iterator bin_it;

		for (	bin_it = (*chr_it).second.begin();
				bin_it != (*chr_it).second.end(); 
				++bin_it ) 

			values.insert( values.end(),
						   (*bin_it).second.begin(),
						   (*bin_it).second.end() );
				
	}

	return values;
}
//}}}

//{{{ template <class T> bool UCSCElement<T>::sort_ucscelement_by_value(
template <class T>
bool
UCSCElement<T>::sort_ucscelement_by_start(UCSCElement<T> e, UCSCElement<T> f)
{
	return e.start < f.start;
}      
//}}}

//{{{ template <class T> bool UCSCElement<T>::sort_ucscelement_by_start(
template <class T>
bool
UCSCElement<T>::sort_ucscelement_by_value(UCSCElement<T> e, UCSCElement<T> f)
{
	return e.value < f.value;
}      
//}}}

//{{{ template <class T> bool UCSCElement<T>::compare_ucscelement_by_value(
template <class T>
bool
UCSCElement<T>::compare_ucscelement_by_value(UCSCElement<T> e, UCSCElement<T> f)
{
	return e.value == f.value;
}      
//}}}
#endif //  __UCSC_BINS_HPP__
