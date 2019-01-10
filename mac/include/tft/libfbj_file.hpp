// body file: libfbj_file.cpp

#ifndef LIBFBJ_FILE
#define LIBFBJ_FILE

#include <string>
#include <vector>
#include <deque>
#include "libfbj_fn.hpp"
#include "libfbj_rowsbuff.hpp"

enum _DF_OPT_NAME {
	SUCCESSIVE_DELIMITERS_AS_ONE,\
	QUOTED,\
	SKIP_NOTES,\
	KEEP_NOTES,\
	SKIP_BLANKS,\
	READ_BLANKS,\
	TRIM_LEADING_WHITESPACES,\
	TRIM_LAST_EMPTY_FIELDS,\
	FORMAT_IRREGULAR,\
	FORMAT_CSV
};

typedef  long long RowNum_t;

class tfile_format {

	typedef std::string string;
	typedef std::vector<string> RowData;
	friend class tabular_file;
	
public:
	enum action_code { None,Expand,Wr_Cont,Continue,Break,Exit };
	// Sorted by insreasing severity and actions, suggest to start from the least.
	//   None     = do nothing
	//   Expand   = expand to the required number
	//   Wr_Cont  = write_r(program.outf,true) then continue to next line
	//   Continue = continue to next line
	//   Break    = stop reading the file
	//   Exit     = exit program
	
	// --------------- setup ---------------
	~tfile_format();
	tfile_format();
	tfile_format(bool is_default);
	tfile_format(const tfile_format& othr);
	tfile_format& operator=(const tfile_format& orig);	
	void reset();									// set to original settings defined in the source code
	void set_delimiters_to_default();				// default="\t "
	void set_delimiters(const string& dels);		// multiple delimiters / no delimiter is allowed
	void set_option(_DF_OPT_NAME opt, bool value);	// set option
	void set_titlelines(RowNum_t l);				// number of title lines
	void set_field_nums(field_numbers& f, const string& message_if_short="", action_code act_if_short=None);
	void set_field_nums(int n, const string& message_if_short="", action_code act_if_short=None);
	void clear_field_nums();						// clear all FN set by set_field_nums()
	void forbid_option(_DF_OPT_NAME opt_name);		// opt_name = _DF_OPT_NAME
	void forbid_option(const string& opt_name);		// opt_name = -d -t --notes-startw delimiting_related
	void forbid_sv_rpt();							// forbid reporting Size_Validity on field_nums
	void enable_sv_rpt();							// enable reporting Size_Validity on field_nums
	void forbid_nf_rpt();							// forbid reporting Num_of_Fields inconsistency
	void enable_nf_rpt();							// enable reporting Num_of_Fields inconsistency
	void read_arguments(std::vector<string>& arg, size_t start=1, bool expect_at_start=false, bool stop_at_unknown=false); // start is 0-based
	int  ReadOpt(std::vector<string>& opt, size_t& argi);// stop at the 1st unknown option; new argi point to its left
	string&	output_del();							// delimiter for output, =input_dels_[0] or an empty string
	string&	comment_sw();							// comment lines start with this string, can't be empty.
	
	// These 2 functions decide how to store data. If none is called, store data internally, where only the last row is stored
	// To use expression_wrapper,  you must set_storage_to(program.main_data()) because e_wrp takes data from program.main_data()
	// To use multiple_conditions, you must set_storage_to(program.main_data()) because MC uses string/math expressions.
	void set_storage_to(RowsBuffer&);				// store data in an external RowsBuffer&
	void set_storage_to(std::deque<RowData>&);		// store data in an external deque of vector<string>
	
	// --------------- access ---------------
	RowNum_t	titlelines();						// number of title lines
	int			min_NumFields();					// previously min_req()
	bool		opt(_DF_OPT_NAME o);				// value of option named o
	bool		is_header(RowNum_t rowID)			{ return rowID < titlelines(); }
	string		help_text();						// help text depending on forbiddance
	
	// --------------- test ---------------
	bool is_content(char c);						// test whether c is content,	mostly internal use
	bool is_delimiter(char c);						// test whether c is delimiter,	mostly internal use
	int  act_by_size(int size, bool& validity);		// return overall action_code based on size, validity is the overall result

private:
	struct tfile_format_data;
	tfile_format_data *fdta;
};

int ReadArg(std::vector<std::string>& args, size_t& argi, tfile_format& f, int ErrCode=-1); // call f.ReadOpt(args,argi); return num_read;
extern std::string DLMTR;									// delimiter for output, default to "\t"
extern tfile_format default_tfile_format;					// only this object's 1st delimiter is DLMTR

//===================== tablular file classes =====================

// features: 1) find and open <filename> <filename>.gz <filename>.bz2, even filename is not ended with .gz .bz2
//			 2) controllable skip-notes skip-blank read-blank ignore-leading-whitespace ignore-trailing-empty-column 
//			 3) robust to all newline characters (Unix, Mac, DOS)
//			 4) robust to no newline charater at the end of the file
//			 5) support >1 delimiters, no delimiters, or treating successive delimiters as one
//			 6) read_r => read a line as a row of a table, fields separated by delimiters;
//			 7)	read_l => read through a whole line without breaking.
//			 8) User specified number of header lines
//			 9) Automatically check for the number of fields across rows, and decide what to do if they varies
//			10) Store data internally / externally, for only the last line / the fast few lines / all lines.
//			11) Use default_tfile_format or fmt if fmt!=NULL; good for unanimous format through out the program.
//			12) Can show progress; good for running on a remote machine in interactive mode (no more broken pipe).

#define Rows_in_File(obj, name, fmt)					tabular_file (obj)(name,fmt);							(obj).read_r(); (obj).next()
#define Rows_after_N_rows_in_File(obj, n, name, fmt)	tabular_file (obj)(name,fmt,n);							(obj).read_r(); (obj).next()
#define every_M_Rows_in_File(m, obj, name, fmt)			tabular_file (obj)(name,fmt,0,m-1);						(obj).read_r(); (obj).next()
#define every_M_Rows_after_N_in_File(m,obj,n,name,fmt)	tabular_file (obj)(name,fmt,n,m-1);						(obj).read_r(); (obj).next()
#define Lines_in_File(obj, name, fmt)					tabular_file (obj)(name,fmt);							(obj).read_l(); (obj).next()
#define Lines_after_N_lines_in_File(obj, n, name,fmt)	tabular_file (obj)(name,fmt,n);							(obj).read_l(); (obj).next()
#define every_M_Lines_in_File(m, obj, name, fmt)		tabular_file (obj)(name,fmt,0,m-1);						(obj).read_l(); (obj).next()
#define every_M_Lines_after_N_in_File(m,obj,n,name,fmt)	tabular_file (obj)(name,fmt,n,m-1);						(obj).read_l(); (obj).next()
#define first_N_Rows_in_File(n, obj, name, fmt)			tabular_file (obj)(name,fmt);(obj).RowNumber()+1<(n) && (obj).read_r(); (obj).next()
#define first_N_Lines_in_File(n, obj, name, fmt)		tabular_file (obj)(name,fmt);(obj).RowNumber()+1<(n) && (obj).read_l(); (obj).next()
#define Data_in_File(obj, name, fmt)					tabular_file (obj)(name,fmt,fmt.titlelines());			(obj).read_r(); (obj).next()
// Upper case letters indicate arguments in the same sequence
// The last one requires that fmt is a tfile_format object, not a non-negative integer / NULL that is allowed by the others
// For M: 1) if line l should be skipped, no matter it satisfies fmt or not, it is counted and skipped
//        2) if line l should be read, if fmt says Wr_Cont/Continue, skip another M-1
//        3) about 2), if you want to go to immediate next content line without skipping M-1, obj.next(m-1) <= need some programming
//        4) RowNumber is still 0,1,2..
//        5) no skipping M-1 before the 1st content row/line

class tabular_file
{
	typedef std::string string;
	typedef std::vector<string> RowData;
public:
	// --------------- setup ---------------
	~tabular_file();
	tabular_file();
									tabular_file(const std::string& name, int min_NumOf_Fields, int skip_begin=0, int skip_every=0);
									void    open(const std::string& name, int min_NumOf_Fields, int skip_begin=0, int skip_every=0);
	template<typename StrContainr>	tabular_file(const StrContainr& name, int min_NumOf_Fields, int skip_begin=0, int skip_every=0);
	template<typename StrContainr>	void    open(const StrContainr& name, int min_NumOf_Fields, int skip_begin=0, int skip_every=0);
									tabular_file(const std::string& name, tfile_format* f=NULL, int skip_begin=0, int skip_every=0);
									void    open(const std::string& name, tfile_format* f=NULL, int skip_begin=0, int skip_every=0);
	template<typename StrContainr>	tabular_file(const StrContainr& name, tfile_format* f=NULL, int skip_begin=0, int skip_every=0);
	template<typename StrContainr>	void    open(const StrContainr& name, tfile_format* f=NULL, int skip_begin=0, int skip_every=0);
	// min_NumOf_Fields can be a negative or positive number; while 0 has no effect.
	// With min_NumOf_Fields, action(lack_NF)=continue; with tfile_format*, action depends on f.
	// StrContainr for the constructor can be vector / deque / set of strings.
	// You can't do tabular_file("",0); but can do tabular_file(string(),0).
	// if (skip_begin>0) tabular_file will skip the lines at the beginning of each file. RowNumber() doesn't include these lines.
	
	bool is_open();						// file is opened
	void   close();						// close the currently opened file and remove all others in the queue
	void add_file(const string& file);	// put file at the end of the queue for reading
	
	// --------------- data ---------------
	string FileName();					// currently opened file name
	int	 FileNumber();					// 0-based; 
	RowNum_t RowNumber();				// 0-based; skipped_lines(comments/blanks) does not count; but if --skip-comment=no, comments still count; init=-1
	int  NumFields();					// number of fields
	bool SizeValid();					// number of fields is valid
	bool is_blank();					// the current line is a blank   line (contents().size=1 && contents()[0].empty())
	bool is_header();					// the current line is a header  line (must use external tfile_format to take effect)
	bool is_comment();					// the current line is a comment line
	RowData& contents();				// get the last row , user can access and modify the contents
	string& operator[]( int i );		// get contents()[i], user can access and modify the contents
	void clear_nf();					// clear number of fields counts
	
	// --------------- read ---------------
	void seekg(long offset);			// goto offset from beginning. not work for gz/bz2
	bool next();						// skip m lines, return !eof().
	bool eof();							// file is at EOF
	bool read_r();						// skip comments and blank, then read a row (use delimiters).
	bool read_l();						// skip comments and blank, then read a line (no delimiters).
	
	// --------------- write ---------------
	void write_r(std::ostream& o,					   bool write_endl=false);	// write all fields sequentially
	void row2str(std::string& o,					   bool write_endl=false);	// write all fields sequentially
	void write_r(std::ostream& o, std::vector<int>& f, bool write_endl=false);	// write selected fields, missing = ""
	void write_r(std::ostream& o, field_numbers&    f, bool write_endl=false);	// write selected fields, missing = ""

private:
	struct tabular_file_data;
	tabular_file_data*	tdta;
	tabular_file& operator=(const tabular_file& orig);	// forbid this operation
	tabular_file(const tabular_file& othr);				// forbid this operation
};

template<typename StrContainr>
void tabular_file::open(const StrContainr& filenames, int min_NumOf_Fields, int skip_begin, int skip_every)
{
	if (filenames.empty())
		open(string(),min_NumOf_Fields,skip_begin,skip_every);
	else
	{
		typename StrContainr::const_iterator it(filenames.begin());
		open(*it,min_NumOf_Fields,skip_begin,skip_every);
		for (++it;it!=filenames.end();++it) add_file(*it);
	}
}

template<typename StrContainr>
void tabular_file::open(const StrContainr& filenames, tfile_format* fmtPtr, int skip_begin, int skip_every)
{
	if (filenames.empty())
		open(string(),fmtPtr,skip_begin,skip_every);
	else
	{
		typename StrContainr::const_iterator it(filenames.begin());
		open(*it,fmtPtr,skip_begin,skip_every);
		for (++it;it!=filenames.end();++it) add_file(*it);
	}
}

#endif
