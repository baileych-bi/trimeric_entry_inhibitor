
//Directions:
//To compile: $clang++ -std=c++23 -march=native -O3 -o extract_sequences main.cc
//To run:     $./extract_sequences /path/to/fastq/*.fastq
//            (hint, make sure .fastq files are unzipped first)
//Output will be found in /path/to/fastq/analysis/FASTQ_FILENAME/SAMPLE/ with
//one subfolder for each fastq input file. Sample names and barcodes are
//hardcoded in main(). Every input file will have a (1) summary.txt file
//inidcating the numbers of records parsed and how many were rejected and for
//what reasons, (2) a unique_sequences.fasta file which consists of each
//accepted unique nucleotide sequence, in descending order of occurrence; the
//fasta header row specifies N=XX where XX is the number of times that exact
//sequence occurred in the fastq data, (3) unique_translations.fasta which is
//analagous to unique_sequences.fasta but grouped according to translation
//rather than nucleotide sequence, and (4) a mutations_frequencies.csv file
//which contains the frequency with which each non-wild type amino acid is
//observed at every position of the wild type Env protein after alignment, and
//(5) a top_10_mutations.csv file that contains a list of the 10 most common
//amino acid substitutions.

//For questons, please contact Charles Bailey (baileych@broadinstitute.org)

#include <algorithm>                                                            //standard library headers
#include <functional>
#include <iostream>
#include <filesystem>
#include <format>
#include <fstream>
#include <list>
#include <ranges>
#include <regex>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include <fcntl.h>                                                              //Unix/Linux specific headers (for mmap, etc.)
#include <unistd.h>
#include <sys/mman.h>

/*******************************************************************************
 * Globals
 ******************************************************************************/
static const std::unordered_map<std::string_view, char> TTABLE = {              //simple/slow translation table
    {"AAA", 'K'}, {"AAC", 'N'}, {"AAG", 'K'}, {"AAT", 'N'},
    {"ACA", 'T'}, {"ACC", 'T'}, {"ACG", 'T'}, {"ACT", 'T'},
    {"AGA", 'R'}, {"AGC", 'S'}, {"AGG", 'R'}, {"AGT", 'S'},
    {"ATA", 'I'}, {"ATC", 'I'}, {"ATG", 'M'}, {"ATT", 'I'},
    {"CAA", 'Q'}, {"CAC", 'H'}, {"CAG", 'Q'}, {"CAT", 'H'},
    {"CCA", 'P'}, {"CCC", 'P'}, {"CCG", 'P'}, {"CCT", 'P'},
    {"CGA", 'R'}, {"CGC", 'R'}, {"CGG", 'R'}, {"CGT", 'R'},
    {"CTA", 'L'}, {"CTC", 'L'}, {"CTG", 'L'}, {"CTT", 'L'},
    {"GAA", 'E'}, {"GAC", 'D'}, {"GAG", 'E'}, {"GAT", 'D'},
    {"GCA", 'A'}, {"GCC", 'A'}, {"GCG", 'A'}, {"GCT", 'A'},
    {"GGA", 'G'}, {"GGC", 'G'}, {"GGG", 'G'}, {"GGT", 'G'},
    {"GTA", 'V'}, {"GTC", 'V'}, {"GTG", 'V'}, {"GTT", 'V'},
    {"TAA", '*'}, {"TAC", 'Y'}, {"TAG", '*'}, {"TAT", 'Y'},
    {"TCA", 'S'}, {"TCC", 'S'}, {"TCG", 'S'}, {"TCT", 'S'},
    {"TGA", '*'}, {"TGC", 'C'}, {"TGG", 'W'}, {"TGT", 'C'}, 
    {"TTA", 'L'}, {"TTC", 'F'}, {"TTG", 'L'}, {"TTT", 'F'}
};

static const int32_t BLOSUM62[22][22] = {                                       //BLOSUM62 expanded to include * and - (although we should never align - with anything)
//	     A	 C	 D	 E	 F	 G	 H	 I	 K	 L	 M	 N	 P	 Q	 R	 S	 T	 V	 W	 Y	 *	 -
/*A*/{	 4,	 0,	-2,	-1,	-2,	 0,	-2,	-1,	-1,	-1,	-1,	-2,	-1,	-1,	-1,	 1,	 0,	 0,	-3,	-2,	-4,	-4},
/*C*/{	 0,	 9,	-3,	-4,	-2,	-3,	-3,	-1,	-3,	-1,	-1,	-3,	-3,	-3,	-3,	-1,	-1,	-1,	-2,	-2,	-4,	-4},
/*D*/{	-2,	-3,	 6,	 2,	-3,	-1,	-1,	-3,	-1,	-4,	-3,	 1,	-1,	 0,	-2,	 0,	-1,	-3,	-4,	-3,	-4,	-4},
/*E*/{	-1,	-4,	 2,	 5,	-3,	-2,	 0,	-3,	 1,	-3,	-2,	 0,	-1,	 2,	 0,	 0,	-1,	-2,	-3,	-2,	-4,	-4},
/*F*/{	-2,	-2,	-3,	-3,	 6,	-3,	-1,	 0,	-3,	 0,	 0,	-3,	-4,	-3,	-3,	-2,	-2,	-1,	 1,	 3,	-4,	-4},
/*G*/{	 0,	-3,	-1,	-2,	-3,	 6,	-2,	-4,	-2,	-4,	-3,	 0,	-2,	-2,	-2,	 0,	-2,	-3,	-2,	-3,	-4,	-4},
/*H*/{	-2,	-3,	-1,	 0,	-1,	-2,	 8,	-3,	-1,	-3,	-2,	 1,	-2,	 0,	 0,	-1,	-2,	-3,	-2,	 2,	-4,	-4},
/*I*/{	-1,	-1,	-3,	-3,	 0,	-4,	-3,	 4,	-3,	 2,	 1,	-3,	-3,	-3,	-3,	-2,	-1,	 3,	-3,	-1,	-4,	-4},
/*K*/{	-1,	-3,	-1,	 1,	-3,	-2,	-1,	-3,	 5,	-2,	-1,	 0,	-1,	 1,	 2,	 0,	-1,	-2,	-3,	-2,	-4,	-4},
/*L*/{	-1,	-1,	-4,	-3,	 0,	-4,	-3,	 2,	-2,	 4,	 2,	-3,	-3,	-2,	-2,	-2,	-1,	 1,	-2,	-1,	-4,	-4},
/*M*/{	-1,	-1,	-3,	-2,	 0,	-3,	-2,	 1,	-1,	 2,	 5,	-2,	-2,	 0,	-1,	-1,	-1,	 1,	-1,	-1,	-4,	-4},
/*N*/{	-2,	-3,	 1,	 0,	-3,	 0,	 1,	-3,	 0,	-3,	-2,	 6,	-2,	 0,	 0,	 1,	 0,	-3,	-4,	-2,	-4,	-4},
/*P*/{	-1,	-3,	-1,	-1,	-4,	-2,	-2,	-3,	-1,	-3,	-2,	-2,	 7,	-1,	-2,	-1,	-1,	-2,	-4,	-3,	-4,	-4},
/*Q*/{	-1,	-3,	 0,	 2,	-3,	-2,	 0,	-3,	 1,	-2,	 0,	 0,	-1,	 5,	 1,	 0,	-1,	-2,	-2,	-1,	-4,	-4},
/*R*/{	-1,	-3,	-2,	 0,	-3,	-2,	 0,	-3,	 2,	-2,	-1,	 0,	-2,	 1,	 5,	-1,	-1,	-3,	-3,	-2,	-4,	-4},
/*S*/{	 1,	-1,	 0,	 0,	-2,	 0,	-1,	-2,	 0,	-2,	-1,	 1,	-1,	 0,	-1,	 4,	 1,	-2,	-3,	-2,	-4,	-4},
/*T*/{	 0,	-1,	-1,	-1,	-2,	-2,	-2,	-1,	-1,	-1,	-1,	 0,	-1,	-1,	-1,	 1,	 5,	 0,	-2,	-2,	-4,	-4},
/*V*/{	 0,	-1,	-3,	-2,	-1,	-3,	-3,	 3,	-2,	 1,	 1,	-3,	-2,	-2,	-3,	-2,	 0,	 4,	-3,	-1,	-4,	-4},
/*W*/{	-3,	-2,	-4,	-3,	 1,	-2,	-2,	-3,	-3,	-2,	-1,	-4,	-4,	-2,	-3,	-3,	-2,	-3,	11,	 2,	-4,	-4},
/*Y*/{	-2,	-2,	-3,	-2,	 3,	-3,	 2,	-1,	-2,	-1,	-1,	-2,	-3,	-1,	-2,	-2,	-2,	-1,	 2,	 7,	-4,	-4},
/***/{	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	 1,	-4},
/*-*/{	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	-4,	 1}
};

static const std::string XBB_S = "FVFLVLLPLVSSQCVNLITRTQSYTNSFTRGVYY"
    "PDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPALPFNDGVYFASTEKSNIIRGWIFGTTLDS"
    "KTQSLLIVNNATNVVIKVCEFQFCNDPFLDVYQKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKEGNFK"
    "NLREFVFKNIDGYFKIYSKHTPINLERDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPVDSSSGWTAG"
    "AAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNL"
    "CPFHEVFNATTFASVYAWNRKRISNCVADYSVIYNFAPFFAFKCYGVSPTKLNDLCFTNVYADSFVIRGNEVSQ"
    "IAPGQTGNIADYNYKLPDDFTGCVIAWNSNKLDSKPSGNYNYLYRLFRKSKLKPFERDISTEIYQAGNKPCNGV"
    "AGPNCYSPLQSYGFRPTYGVGHQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNK"
    "KFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTP"
    "TWRVYSTGSNVFQTRAGCLIGAEYVNNSYECDIPIGAGICASYQTQTKSHRRARSVASQSIIAYTMSLGAENSV"
    "AYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLKRALTGIAVEQDKNTQ"
    "EVFAQVKQIYKTPPIKYFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQ"
    "KFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQF"
    "NSAIGKIQDSLSSTASALGKLQDVVNHNAQALNTLVKQLSSKFGAISSVLNDILSRLDKVEAEVQIDRLITGRL"
    "QSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNF"
    "TTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFK"
    "EELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAG"
    "LIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT*";

static const std::string FW_REFERENCE = "bbbbbbNNNNNNCTTGTTAACAACTAAACGAACAATG";//the forward reference abuts the start codon
static const std::string RV_REFERENCE = "ACGAACTTATGGATTTGTTT";                 //the reverse reference abuts the stop codon

/*******************************************************************************
 * Sequence & Alignment Utilities
*******************************************************************************/
/** Map codons onto amino acids. Will throw on invalid codon. */
std::string
translate_dna(std::string_view dna) {
    std::string aas; aas.reserve(dna.size() / 3);
    for (size_t i=0; i + 2 < dna.size(); i += 3)
        aas.push_back(TTABLE.at(dna.substr(i, 3)));
    return aas;
}

/** Align sequences to a subject and record the identities of the amino acids
  * that align to each residue in the subject.
 */
struct Aligner {
    Aligner(std::string_view sbj_aas, uint32_t gapp=4)
    : sbj_(sbj_aas.size())
    , counts_(22 * sbj_aas.size(), 0)
    , gapp_(gapp) {
        std::transform(sbj_aas.begin(), sbj_aas.end(), sbj_.begin(), encode_aa);
    }

    /** Perform Needleman-Wunsch alignment of qry_aas to sbj_aas supplied during
        construction. Uses to global BLOSUM62 scoring matrix. No penalty for 
        leading or trailing gaps is applied. By default, it only records which
        amino acids align to each position in the subject with each call,
        data which can be accessed by calling aa_count(). Each time an amino
        is recorded, it will be counted 'multiplicity' times. If build_strings
        is true, then human-readable gapped and aligned query and subject 
        strings can be recovered by calling strings().
     */
    void
    operator()(std::string_view qry_aas, 
        uint32_t multiplicity=1, 
        bool build_strings=false) {

        qstring_.clear();                                                       
        sstring_.clear();

        if (qry_aas.size() == 0 || sbj_.size() == 0)                            //nothing to align in case of an empty string
            return;

        struct Cell {                                                           //element of traceback matrix; has score and pointer to previous cell
            int32_t score = 0;
            const Cell *prv = nullptr;
        };

        alignments_done_ += multiplicity;

        thread_local std::vector<uint32_t> qry;
        qry.clear();
        for (char aa : qry_aas)
            qry.push_back(encode_aa(aa));

        thread_local std::vector<Cell> cells;                                   //init traceback matrix
        cells.clear();
        const size_t qrows = qry .size() + 1;
        const size_t scols = sbj_.size() + 1;
        cells.resize(qrows * scols);
        for (size_t i=1; i < qrows; ++i) {
            cells[i * scols].score = -i * gapp_;
            cells[i * scols].prv = &cells[(i - 1) * scols];                     //no penalty for leading gaps
        }

        for (size_t j=1; j < scols; ++j) {
            cells[j].score = -j * gapp_;
            cells[j].prv = &cells[j - 1];                                       //no penalty for leading gaps
        }

        Cell *cur = nullptr;                                              
        for (size_t i = 1; i < qrows; ++i) {                                    
            uint32_t qc = qry[i - 1];                                           
            for (size_t j = 1; j < scols; ++j) {                                //incr cur from cells[i, j] to cells[i, j+1]
                uint32_t sc = sbj_[j - 1];
                cur = &cells[i * scols + j];
                                                                                //our three possible traceback origins:
                const Cell *ul = cur - scols - 1;                               //up + left
                const Cell *u  = cur - scols;                                   //up
                const Cell *l  = cur - 1;                                       //left

                cur->score = ul->score + BLOSUM62[qc][sc];                      //we start by assuming that qry and sbj should match
                cur->prv = ul;                                                  

                int32_t alt_score;                                              //and then consider whether insertion/deletion is better

                alt_score = u->score - gapp_;                                   //no penalty for trailing gaps
                if (alt_score > cur->score) {
                    cur->score = alt_score;           
                    cur->prv = u;
                }

                alt_score = l->score - gapp_;                                   //no penalty for trailing gaps
                if (alt_score > cur->score) {
                    cur->score = alt_score;                                     
                    cur->prv = l;
                }
            }
        }

        const uint32_t *qc = qry .data() + qry .size();
        const uint32_t *sc = sbj_.data() + sbj_.size();
        size_t aln_pos = sbj_.size();

        const uint32_t gap_code = encode_aa('-');
        const Cell *last = &cells[cells.size() - 1];
        for (; last->prv; last = last->prv) {
            if (last - last->prv == scols + 1) {                                //came from up/left => match
                counts_[22 * (--aln_pos) + *--qc] += multiplicity;
                --sc;
            } else if (last - last->prv == scols) {                             //came from up      => insertion
                //we don't record insertions yet
                --qc;
            } else {                                                            //came from left    => deletion
                counts_[22 * (--aln_pos) + gap_code] += multiplicity;
                --sc;
            }
        }

        if (build_strings) {
            qc = qry .data() + qry .size();
            sc = sbj_.data() + sbj_.size();

            qstring_.reserve(std::max(qry.size(), sbj_.size()));
            sstring_.reserve(std::max(qry.size(), sbj_.size()));

            const Cell *last = &cells[cells.size() - 1];
            for (; last->prv; last = last->prv) {
                if (last->prv + scols + 1 == last) {                            //came from up/left => match
                    qstring_.push_back(decode_aa(*--qc));
                    sstring_.push_back(decode_aa(*--sc));
                } else if (last->prv + scols == last) {                         //came from up      => insertion
                    qstring_.push_back(decode_aa(*--qc));
                    sstring_.push_back('-'             );
                } else {                                                        //came from left    => deletion
                    qstring_.push_back('-'             );
                    sstring_.push_back(decode_aa(*--sc));
                }
            }

            std::reverse(qstring_.begin(), qstring_.end());
            std::reverse(sstring_.begin(), sstring_.end());
        }        
    }

    uint32_t aa_count(size_t pos, char aa) const {
        return counts_[22 * pos + encode_aa(aa)];
    }

    size_t alignments_done() const {
        return alignments_done_;
    }

    std::pair<std::string_view, std::string_view>
    strings() const { 
        return std::make_pair(
            std::string_view(qstring_),
            std::string_view(sstring_)
        );
    }

    //map amino acids onto unsigned [0, 19], '*' onto 20, and '-' onto 21
    static uint32_t encode_aa(char aa) {
        static const uint32_t lut[26] = {
        //     A     ?     C     D     E     F     G
            0x00, 0xFF, 0x01, 0x02, 0x03, 0x04, 0x05,
        //     H     I     ?     K     L     M     N
            0x06, 0x07, 0xFF, 0x08, 0x09, 0x0A, 0x0B,
        //     O     P     Q     R     S     T     U
            0xFF, 0x0C, 0x0D, 0x0E, 0x0F, 0x10, 0xFF,
        //     V     W     X     Y     Z
            0x11, 0x12, 0xFF, 0x13, 0xFF 
        };
        uint32_t code = 0;
        code += ('A' <= aa) * lut[aa - 'A'];
        code += (aa == '*') * 20;
        code += (aa == '-') * 21;
        return code;
    }

    //map number from [0..21] onto amino acid, stop, or gap char
    static char
    decode_aa(uint32_t code) {
        static const char *aas = "ACDEFGHIKLMNPQRSTVWY*-";
        return aas[code];
    }

    Aligner
    operator+=(Aligner &rhs) {
        alignments_done_ += rhs.alignments_done_;
        for (size_t i=0; i != counts_.size(); ++i)
            counts_[i] += rhs.counts_[i];
        return *this;
    }

private:
    size_t alignments_done_ = 0;
    std::vector<uint32_t> sbj_;                                                 //the subject amino acid sequence, encoded
    std::vector<uint32_t> counts_;                                              //counts_[i, j] := #times aa w/code j aligns with aa at sbj_[i]
    int32_t gapp_ = 4;                                                          //the NW gap penalty (greater numbers penalize more)
    std::string qstring_, sstring_;                                             //sotres the human-readable alignment from the latest call, if build_strings = true
};

/*******************************************************************************
 * Fastq parsing structures and functions.
*******************************************************************************/
/** Holds two string_views for the sequence and quality lines of a fastq file */
struct Record {
    Record(std::string_view sequence, std::string_view quality)
    : seq_(sequence)
    , qly_(quality) { }

    std::string_view sequence() const { return seq_;     }
    std::string_view quality()  const { return qly_;     }

    std::string &umi()                { return umi_;     }
    std::string &barcode()            { return barcode_; }
    std::string_view umi()      const { return umi_;     }
    std::string_view barcode()  const { return barcode_; }

    size_t size() const { return seq_.size(); }

    void trim(size_t pos, size_t len=std::string::npos) {
        seq_ = seq_.substr(pos, len);
        qly_ = qly_.substr(pos, len);
    }

private:
    std::string_view seq_;
    std::string_view qly_;
    std::string umi_;
    std::string barcode_;
};

/** Memory maps the file supplied during construction and automatically unmaps
  * it when destroyed. 
  */
struct MappedFile {
    MappedFile(const std::filesystem::path &path) { map(path); }
    ~MappedFile() { unmap(); };

    void
    unmap() {
        if (begin_) {
            munmap((void *)begin_, end_ - begin_);
            begin_ = end_ = nullptr;
            path_ = "";
        }
    }

    bool
    map(const std::filesystem::path &path) {
        if (*this) unmap();

        int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) return false;
        const size_t sz = lseek(fd, 0, SEEK_END);

        begin_ = reinterpret_cast<char *>(
            mmap(NULL, sz, PROT_READ, MAP_SHARED, fd, 0)
        );

        close(fd);

        if (begin_ == (void *)(-1)) {
            begin_ = nullptr;
            return false;
        } else {
            madvise((void *)begin_, sz, MADV_SEQUENTIAL | MADV_HUGEPAGE);
            end_ = begin_ + sz;
            path_ = path;
            return true;
        }
    }

    operator bool() { return bool(begin_); }

    const char *begin() const { return begin_; }
    const char *end()   const { return end_;   }

    size_t size() const { return end_ - begin_; }
private:
    std::filesystem::path path_;
    const char *begin_ = nullptr, *end_ = nullptr;
};

/** Tasks are any functors bindable by a std::function<bool(Record &)>. They are
  * called on a Record and may possibly modify it. They return true if the
  * Record is deemed acceptable and/or successfully modified or false otherwise.
 */
namespace tasks {
/** Ensure .fastq sequence data is capital ATGC only. */
struct RejectInvalidChars {
    bool operator()(Record &record) const {
        return record.sequence().find_first_not_of("ACGT") == std::string::npos;
    }
};

/** Is the nucleic acid sequence length in [min, max]? */
struct RequireLength {
    RequireLength(size_t len) : min_length_(len), max_length_(len) { }
    RequireLength(size_t min, size_t max) : min_length_(min), max_length_(max) {
        if (max < min) {
            max_length_ = min_length_ = 0;
            throw std::invalid_argument("RequireLength called with max < min");
        }
    }

    bool operator()(Record &record) const { 
        return min_length_ <= record.size() && record.size() <= max_length_;
    }

private:
    size_t min_length_ = 0;
    size_t max_length_ = std::numeric_limits<size_t>::max();
};

/** 
 * Ensure presence of reference sequence and trim, removing the reference and
 * its 5' or 3' flanking sequence as determined by is_reverse.
 */
struct TrimReferenceExcl {
    TrimReferenceExcl(std::string_view reference, bool is_reverse)
    : is_reverse_(is_reverse) {
        if (reference.find_first_not_of("ACGTNub") != std::string::npos) {
            std::cerr << reference << std::endl;
            throw std::runtime_error("invalid reference sequence");
        }

        umi_msk_.resize(reference.size(), false);
        bar_msk_.resize(reference.size(), false);
        re_str_.reserve(reference.size());
        for (size_t i=0; i != reference.size(); ++i) {
            char c = reference[i];
            switch (c) {
                case 'u':
                    re_str_.push_back('.');
                    umi_msk_[i] = true;
                break;
                case 'b':
                    re_str_.push_back('.');
                    bar_msk_[i] = true;
                break;
                case 'N':
                    re_str_.push_back('.');
                break;
                default:
                    re_str_.push_back(c);
            }
        }
        re_ = std::regex(re_str_);
    }

    bool operator()(Record &record) const {
        std::string_view sbj = record.sequence();
        std::cmatch match;
        bool found = std::regex_search(sbj.data(), sbj.data() + sbj.size(), 
            match, re_);

        if (found) {
            size_t pos = match.position(0);
            for (size_t i=0; i != re_str_.size(); ++i) {
                char c = sbj[pos + i];
                if (umi_msk_[i])
                    record.umi().push_back(c);
                else if (bar_msk_[i])
                    record.barcode().push_back(c);
            }
            
            if (is_reverse_)
                record.trim(0, match.position());
            else 
                record.trim(match.position() + match.length());
            return true;
        } else {
            return false;
        }
    }
private:
    bool is_reverse_  = false;
    std::vector<bool> umi_msk_;
    std::vector<bool> bar_msk_;
    std::string re_str_;
    std::regex re_;
};

/**
 */
 struct FilterInvalidBarcodes {
    FilterInvalidBarcodes(
        const std::vector<std::pair<std::string, std::string>> &samples) {
            for (const auto &[sample_name, barcode] : samples)
                valid_barcodes_.insert(barcode);
         }

    bool operator()(Record &record) {
        return valid_barcodes_.find(record.barcode()) != valid_barcodes_.end();
    }
private:
    std::unordered_set<std::string> valid_barcodes_;
};

/** Require that a DNA sequence be a multiple of 3, start with a Met codon
  * and end with a stop codon.
 */
struct RequireORFLike {
    bool is_stop(const char *a) {                                               //a slow, branch-y stop codon classifier
        if (*a != 'T')
            return false;

        const char *b = a + 1, *c = a + 2;
        if (*b == 'A' && (*c == 'A' || *c == 'G'))
            return true;
        if (*b == 'G' && (*c == 'A'))
            return true;
        return false;
    }

    bool operator()(Record &record) {
        if (record.sequence().size() % 3 != 0 || record.sequence().empty())     //ensure length is multiple of 3
            return false;

        const char *first = record.sequence().data();
        const char *last  = first + record.sequence().size();
        for (; first + 3 != last; first +=3)                                    //require no stop codons before last codon
            if (is_stop(first)) {
                //std::cout << "not orf #1" << std::endl;
                //std::cout << record.sequence() << std::endl;
                return false;
            }

        if (!is_stop(first)) {                                                  //also require last codon IS stop codon
            //std::cout << "not orf #2" << std::endl;
            return false;
        }

        return true;
    }
};

};

struct Count {
    bool operator==(const Count &rhs) const { return n == rhs.n; }
    bool operator< (const Count &rhs) const { return n <  rhs.n; }
    Count &operator+=(const Count &rhs) { n += rhs.n; return *this; }
    operator       size_t &  () { return n; }
    operator const size_t &  () { return n; }
    operator       size_t && () { return std::move(n); }
    size_t n = 0;
};

template<typename K, typename Hash=std::hash<K>, typename Equal=std::equal_to<K>>
struct Counter {
    auto begin()        { return counts_.begin(); }
    auto end()          { return counts_.end();   }
    auto begin() const  { return counts_.begin(); }
    auto end()   const  { return counts_.end();   }

    size_t size() const { return counts_.size(); }

    Counter &merge(Counter &&ctr) {
        for (auto &[k, ct] : ctr) counts_[k] += ct;
        return *this;
    }

    template<typename T>
    Count &operator[](const T &key) {
        auto ii = counts_.find(key);
        if (ii == counts_.end()) {
            auto [jj, ins] = counts_.insert(std::pair<K, Count>{key, Count(0)});
            return jj->second;
        } else {
            return ii->second;
        }
    }

    template<typename T>
    Count operator[](const T &key) const {
        auto ii = counts_.find(key);
        if (ii == counts_.end())
            return Count(0);
        return ii->second;
    }

private:
    std::unordered_map<K, Count, Hash, Equal> counts_; 
};

struct StringHash {                                                             //we need this special hash to search for string keys with string_view queries
    using is_transparent = void;

    size_t operator()(std::string_view sv) const { 
        return std::hash<std::string_view>{}(sv);
    }

    size_t operator()(const std::string &st) const {
        return std::hash<std::string>{}(st);
    }
};

/** Analyze records in a chunk of fastq. Collect those that pass the filters
  * and record those that fail, and which filters they fail.
 */
struct WorkOrder {
    using Task = std::function<bool(Record &)>;                                 //a Task filters or modifies a Record; returns false if Record should be discarded 

    WorkOrder(std::span<const Task> tasks,                                      //[begin, end) refers to 1st and sentinel byte in a mmap'd file
        const char *begin, const char *end, 
        const char *lo,    const char *hi)                                      //[lo, hi) is the chunk of said file this WorkOrder will process
    : tasks_(tasks), failures_(tasks_.size(), 0)
    , begin_(begin), end_(end)
    , lo_(lo), hi_(hi) { 
        lo_ = snap_to_header(lo_);                                              //align lo and hi to 1st byte of the next header line
        hi_ = snap_to_header(hi_);                                              //note: our lo_ is the same as the previous WorkOrder's hi_
    }

    void
    operator()() {
        const char *delim = nullptr, *ptr = lo_;
        while (ptr != hi_) {
            delim = std::find(ptr, hi_, '\n');                                  //skip header
            ptr = delim + (delim != hi_);

            delim = std::find(ptr, hi_, '\n');                                  //extract sequence
            std::string_view sequence(ptr, delim - ptr);
            ptr = delim + (delim != hi_);

            delim = std::find(ptr, hi_, '\n');                                  //skip +
            ptr = delim + (delim != hi_);

            delim = std::find(ptr, hi_, '\n');                                  //extract quality
            std::string_view quality(ptr, delim - ptr);
            ptr = delim + (delim != hi_);

            if (sequence.size() == 0)
                throw std::runtime_error("empty sequence error!");

            if (sequence.size() != quality.size())
                throw std::runtime_error("sequence and quality length mismatch");

            Record record(sequence, quality);

            size_t failure = size_t(-1);
            for (size_t i = 0; i != tasks_.size(); ++i)                         //try to run all tasks in order, stopping if a record 'fails' a Task
                if (!tasks_[i](record)) {                                       //record which Task, if any, the Record fails
                    failure = i;
                    break;
                }

            if (failure < failures_.size()) {                                   //on failure record which task failed
                ++failures_[failure];
            } else {
                sequences_[record.barcode()][record.sequence()] += 1;
            }
        }
    }

    const char *
    snap_to_header(const char *ptr) {                                           //align ptr to first byte of next header, or EOF sentinel if there is no next header
        if (ptr == begin_) return ptr;                                          //first byte is always start of header line

        const char *plus_line = "\n+\n";                                        //otherwise look for the + on a line by itself
        const char *landmark = std::search(ptr, end_, 
            plus_line, plus_line + 3);                                          

        if (landmark == end_) return landmark;                                  //reached EOF

        ptr = std::find(landmark + 3, end_, '\n');                              //look for the newline of the proceeding quality line
        if (ptr == end_)
            throw std::runtime_error("malformed fastq file");                   //said newline should always exist
        return ptr + 1;                                                         //header (or EOF) should come after it
    }
    
    static auto
    into_results(WorkOrder &&order) {                                           //break the WorkOrder into its failures vector and sequences map
        return std::make_pair(std::move(order.failures_),
                              std::move(order.sequences_)
        );
    }

    size_t failure_count(size_t n) const { return failures_[n]; }               //number of records rejected by task n

private:
    const char *begin_ = nullptr, *end_ = nullptr;                              //begin and end are the pointers to our file mapping
    const char *lo_    = nullptr, *hi_  = nullptr;                              //lo and hi are the pointers to the chunk we're assigned

    std::span<const Task> tasks_;                                               //our ordered list of tasks
    std::vector<size_t> failures_;                                              //keep count of which tasks were failed
    std::unordered_map<std::string,
        Counter<std::string, StringHash, std::equal_to<>>,
        StringHash, std::equal_to<>
    > sequences_;                                                               //all the unique sequences we eventually accept and the #times they occur
};

std::vector<WorkOrder>
create_work_orders(const MappedFile &map,                                       //break a mmap'd file into contiguous chunks
                   std::span<const WorkOrder::Task> tasks) {                    //and construct a WorkOrder for each chunk
    size_t n_threads = std::thread::hardware_concurrency();
    std::vector<WorkOrder> work; work.reserve(n_threads);

    const char *begin = map.begin(), *end = map.end();
    const char *lo = begin, *hi = begin;
    const size_t chunk = map.size() / n_threads;
    for (size_t i=0; i != n_threads; ++i) {
        lo = begin + i * chunk;
        hi = end - lo < chunk ? end : lo + chunk;
        work.emplace_back(tasks, begin, end, lo, hi);
    }
    return work;
}

/** Print a summary of the run. */
void
summarize(std::ostream &os, 
    const std::filesystem::path &fastq_path,
    size_t n_unique_sequences,
    size_t n_unique_translations,
    std::span<const size_t> failures,
    std::span<const std::string> excuses) {
    size_t total_failures =                                                     //sum of all rejected reads
        std::ranges::fold_left(failures, size_t(0), std::plus<>());
    os << "Summary of " << fastq_path.string() << ':' << std::endl;
    os << "Records Parsed\t" 
              << (n_unique_sequences + total_failures) << std::endl;
    for (size_t i = 0; i != failures.size(); ++i)
        os << '(' << i << ") " << excuses[i] << "\t" << failures[i] << std::endl;
        os << "Unique nucleotide sequences\t" << n_unique_sequences    << std::endl;
        os << "Unique amino acid sequences\t" << n_unique_translations << std::endl;
}

void
analyze_file(const std::filesystem::path &fastq_path,
             const std::vector<std::pair<std::string, std::string>> &samples) {
    using Task = WorkOrder::Task;
                                                                                //here we define our list of processing/filtering steps ('tasks')
    Task t1 = tasks::RejectInvalidChars();                                      //forbid any sequence with non-ACGTN characters
    Task t2 = tasks::RequireLength(2610, size_t(-1));                           //require untrimmed reads with at least 2652 base pairs
    Task t3 = tasks::TrimReferenceExcl(FW_REFERENCE, false);                    //trim read to forward reference
    Task t4 = tasks::TrimReferenceExcl(RV_REFERENCE, true);                     //trim read to reverse reference
    Task t5 = tasks::FilterInvalidBarcodes(samples);                            //drop all sequences associated with unrecognized barcodes
    Task t6 = tasks::RequireLength(2610, size_t(-1));                           //require trimmed read to still have at least 2652 bp
    Task t7 = tasks::RequireORFLike();                                          //discard reads containing PTCs except those in the final 3 bp

    std::vector<std::function<bool(Record &)>> tasks({t1, t2, t3, t4, t5, t6, t7}); 

    std::vector<size_t> failures(tasks.size(), 0);                              //record #reads failing task n
    std::vector<std::string> excuses = {                                        //description of failure on task n
        "contained non-ACGT characters",
        "untrimmed read less than 2607 bp",
        "forward primer sequence not found",
        "reverse primer sequence not found",
        "read had unrecognized barcode",
        "trimmed read less than 2607 bp",
        "is not a valid ORF"
    };

    std::ios_base::sync_with_stdio(false);                                      //hopefully improve std::ostream performance
    std::unordered_map<std::string,
        Counter<std::string, StringHash, std::equal_to<>>
    > unique_sequences;

    MappedFile mapping(fastq_path.c_str());
    if (!mapping)
        throw std::runtime_error(
            std::format("Failed to mmap file '{}'", fastq_path.c_str())
        );

    std::vector<WorkOrder> work = create_work_orders(mapping, tasks);           //mapping gets brokwn up and processed in parallel

    std::list<std::thread> threads;                                             //invoke WorkOrders in parallel here
    for (size_t i=0; i != work.size(); ++i)
        threads.emplace_back(std::ref(work[i]));
    for (std::thread &th : threads)
        th.join();
    threads.clear();

    size_t new_index = 0;

    for (WorkOrder &order : work) {                                             //collate all completed WorkOrders
        auto [failed, passed] = WorkOrder::into_results(std::move(order));
        for (auto &[barcode, counts] : passed)
            unique_sequences[barcode].merge(std::move(counts));

        for (size_t i=0; i != failures.size(); ++i)                             //and tally the failures
            failures[i] += failed[i];
    }

    mapping.unmap();                                                            //clean up a bit
    work.clear();

    std::cout << "Found " << unique_sequences.size() << " barcodes" << std::endl;
    std::cout << "Failures:" << std::endl;
    for (size_t i=0; i != failures.size(); ++i)
        std::cout << excuses[i] << " - " << failures[i] << std::endl;

    std::unordered_map<std::string, std::string> bc_sample_map;
    for (const auto &[sample, barcode] : samples)
        bc_sample_map[barcode] = sample;

    for (const auto &[barcode, counts] : unique_sequences) {
        std::vector<std::pair<std::string_view, Count>> sorted_dna(                //sort by most common
            counts.begin(), counts.end()
        );

        std::sort(sorted_dna.begin(), sorted_dna.end(), 
            [](const std::pair<std::string_view, Count> &a,
               const std::pair<std::string_view, Count> &b) {
                return b.second < a.second;
            }
        );

        std::filesystem::path output_dir =                                          //create, e.g., /path/to/fastq/analysis/sample01/
            fastq_path.parent_path() / "analysis";
        std::filesystem::create_directory(output_dir);
        output_dir /= fastq_path.stem();
        std::filesystem::create_directory(output_dir);
        output_dir /= bc_sample_map[barcode];
        std::filesystem::create_directory(output_dir);

        std::cout << "writing results for barcode {" << barcode << "} to "
                  << (output_dir).string() << std::endl;

        std::ofstream accepted_file(output_dir / "unique_sequences.fasta");         //write reads to file that passed all filters
        for (const auto &[sequence, count] : sorted_dna)
            accepted_file << ">N=" << (const size_t &)(count) << '\n' << sequence << '\n';
        accepted_file.close();

        std::unordered_map<std::string, size_t> unique_translations;                //now collate unique translations

        for (const auto &[dna, count] : sorted_dna) {                                   
            auto [ii, _] = unique_translations.insert({translate_dna(dna), 0});
            ii->second += (const size_t &)(count);
        }

        std::vector<std::pair<std::string, size_t>> sorted_aas(                     //and sort by occurrence
            unique_translations.begin(), unique_translations.end()
        );

        std::sort(
            sorted_aas.begin(), sorted_aas.end(),
            [](const std::pair<std::string, size_t> &a,
            const std::pair<std::string, size_t> &b) {
                return b.second < a.second;
            }
        );

        std::ofstream translations_file(output_dir / "unique_translations.fasta");  //write translations to file          
        for (const auto &[sequence, count] : sorted_aas) {
            translations_file << ">N=" << count << '\n' << sequence << '\n';
        }
        translations_file.close();

        auto batch_align = [](Aligner &aln,                                         //the next step is to align each unique translation to wt Env
            std::span<std::pair<std::string, size_t>> sbjs) {                       //this lambda will be called in a worker thread on all sequences in sbjs
            for (const auto &[aas, count] : sbjs)                                   
                aln(aas, count);                                                    //note the count paramter ensures that we account for multiple occurrences of each sequence
        };

        std::vector<Aligner> aligners(std::thread::hardware_concurrency(),          //do alignments in parallel with results accumulating in the aligners vector
            Aligner(XBB_S));                                             
        const size_t alignment_batch_size =                                         
            (sorted_aas.size() + aligners.size() - 1) / aligners.size();            //divvy up the work into batches; the last actual batch will usually be a bit smaller
        std::span<std::pair<std::string, size_t>> all_subjects(sorted_aas);
        threads.clear();
        for (size_t i=0; i != std::thread::hardware_concurrency(); ++i) {
            const size_t lo = i * alignment_batch_size;
            const size_t sz = std::min(alignment_batch_size,                        //account for last batch being smaller than the others
                all_subjects.size() - lo);
            threads.emplace_back(batch_align, std::ref(aligners[i]), 
                all_subjects.subspan(lo, sz)
            );
        }
        for (auto &th : threads)
            th.join();

        for (size_t i=1; i < aligners.size(); ++i)                                  //collate results in the first aligner structure
            aligners[0] += aligners[i];

        std::map<float, std::string> top10;

        Aligner &aln = aligners[0];                                                 
        std::ofstream subs_file(output_dir / "mutation_frequencies.csv");
        const std::string symbols = "ACDEFGHIKLMNPQRSTVWY*-";                       //our row headers in proper order
        for (const auto [i, aa] : std::views::enumerate(XBB_S))          //print the wild type residues/column headers (M1, G2, C3... etc)
            subs_file << '\t' << aa << (i + 1);
        subs_file << '\n'; 
        for (char aa : symbols) {                                                   //now go 1-by-1 through the amino acids
            subs_file << aa;                                                        //writing a row for each, with the aa as the row header
            for (const auto [i, wt] : std::views::enumerate(XBB_S)) {
                float freq = aa == wt
                        ? 0.0f                                                   //write the frequency of the wt amino acid as 0
                        : float(aln.aa_count(i, aa)) /                           
                            std::max(aln.alignments_done(), size_t(1));            //prevent divide by 0 if, for some reason, the aligner did nothing
                top10.insert({freq, std::format("{}{}{}", wt, i, aa)});
                subs_file << '\t' << freq;
            }
            subs_file << '\n';
        }
        subs_file.close();

        std::ofstream top10_file(output_dir / "top_10_mutations.csv");
        size_t limit = std::min(size_t(10), top10.size());
        auto ri = top10.rbegin();
        for (size_t i=0; i<limit; ++i, ++ri)
            top10_file << ri->second << '\t' << ri->first << '\n';
        top10_file.close();

        std::ofstream summary_file(output_dir / "summary.txt");                     //finally, write the summary to file and stdout
        summarize(summary_file, fastq_path, counts.size(), 
                unique_translations.size(), failures, excuses);
        summary_file.close();

        summarize(std::cout, fastq_path, counts.size(),
                unique_translations.size(), failures, excuses);
    }
}

int 
main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " [file.fastq...]" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::filesystem::path> fastq_files;
    for (int argn=1; argn < argc; ++argn) {
        fastq_files.push_back(argv[argn]);
        std::ifstream test(fastq_files.back());
        if (!test) {
            std::cerr << "Could not open file '" << fastq_files.back().string() 
                      << "' for reading" << std::endl;
            return EXIT_FAILURE;
        } 
    }

    std::vector<std::pair<std::string, std::string>> samples = {
        {"PBS-day2-M1", "CGTGAT"},
        {"PBS-day2-M2", "ACATCG"},
        {"PBS-day2-M3", "GCCTAA"},
        {"PBS-day2-M4", "TGGTCA"},
        {"PBS-day2-M5", "CACTGT"},
        {"PBS-day5-M1", "ATTGGC"},
        {"PBS-day5-M2", "GATCTG"},
        {"PBS-day5-M3", "TCAAGT"},
        {"PBS-day5-M4", "AGCGAG"},
        {"PBS-day5-M5", "AAGCTA"},
        {"S22-day2-M1", "GTAGCC"},
        {"S22-day2-M2", "TACAAG"},
        {"S22-day2-M3", "GGACGG"},
        {"S22-day2-M4", "GCGGAC"},
        {"S22-day2-M5", "TTTCAC"},
        {"S22-day5-M1", "GGCCAC"},
        {"S22-day5-M2", "CGAAAC"},
        {"S22-day5-M3", "CGTACG"},
        {"S22-day5-M4", "AGGAAT"},
        {"S22-day5-M5", "AGTTTC"},
        {"HCQ-day2-M1", "GAACCT"},
        {"HCQ-day2-M2", "GCCCAG"},
        {"HCQ-day2-M3", "TGACAG"},
        {"HCQ-day2-M4", "CATCAC"},
        {"HCQ-day2-M5", "CTGGAG"},
        {"HCQ-day5-M1", "GATCCG"},
        {"HCQ-day5-M2", "AACACC"},
        {"HCQ-day5-M3", "GTGACG"},
        {"HCQ-day5-M4", "ACAGGA"},
        {"HCQ-day5-M5", "TTACCG"},
        {"Control"    , "TCGTCT"},
    };

    for (const auto &path : fastq_files) {
        std::cout << "Parsing " << path << "..." << std::endl;
        analyze_file(path, samples);
    }

    return EXIT_SUCCESS;
}
