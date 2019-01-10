// body file: libfbj_program.cpp

#ifndef LIBFBJ_DFEDITORS
#define LIBFBJ_DFEDITORS

#include <string>
#include <vector>
#include <boost/regex.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include "libfbj_base.hpp"		// for extract_set()
#include "libfbj_rowsbuff.hpp"	// for RowsBuffer
#include "libfbj_log.hpp"		// for lns
#include "libfbj_elog.hpp"		// for elog
#include "libfbj_int.hpp"		// for IntRanges

class ProgramHandle {

public:
	bool	quiet;										// default is 0, run quietly
	int		nt;											// default is 1, number of threads
	int		nudge;										// default is 0, internal use only
	boost::iostreams::filtering_ostream outf;			// default is stdout
	boost::regex::flag_type				regex_syntax;	// default is POSIX basic
	std::string							trademark;		// default is empty, show in help text header
	std::string							manual;
	
	~ProgramHandle();
	ProgramHandle();
	std::string		name();			// program's original name, not the executable file's name
	std::string		prefix();		// prefix for input and/or output
	std::string		commands();		// one option per line
	std::string		out_name();		// current output file name
	std::string		exe_name();		// name of the executable file
	std::string		exe_path();		// path to the executable file
	std::string		init_path();	// path to the initial directory
	std::string		conf_path();	// path to the configure files
	RowsBuffer&		main_data();	// data storage
	bool			no_web();
	void check_help_request();
	void check_help_at_arg();
	void set_prefix			(const std::string& prefix);
	void set_check_url		(const std::string& URL1, const std::string& URL2, const std::string& version);
	void set_no_web			(const bool noweb);
	void read_arguments		(int argc, char * const argv[], bool rpl_esc=true, bool MustHaveArg=false); // for main(), call once
	void process_arguments	(bool rpl_esc=true);							// internal use only
	void print_messages		(std::ostream& os,std::string message);			// internal use only
	void print_help_text	(std::ostream& os);								// print help text to cout
	void show_version		(std::ostream& os);								// print version number and web check the latest version
	void push_back_help		(std::string S, bool up_front=false);			// add S to help text, up front / at the end
	void help_text_var		(const std::string& Fr, const std::string& To);	// replace Fr to To in help text. Multiplicity allowed.
	void re_open_output		(const std::string& output_filename);			// by default, output to standard output
	void forbid_option		(const std::string& option);					// by default, -h -o -q --tft-xxxx are  enabled.
	void enable_option		(const std::string& option);					// by default, -E -P --nt --prefix --spin-per --no-web are disabled.
	std::vector<std::string>& arg();										// remaining arguments after read_arguments()
	
private:
	struct ProgramData;
	ProgramData * d;
	ProgramHandle(const ProgramHandle& othr);
	ProgramHandle& operator=(const ProgramHandle& orig);	
};

extern ErrorLogger elog;
extern logger lns;
extern ProgramHandle program;

void NeedArg(const std::vector<std::string>& args, int need_how_many_arguments_after_argi, int argi); // determined by args.size()

template <typename T> // T: char float double long_double [u]int8/16/32/64_t (including short,int,long,unsigned,etc)
int ReadStr(const std::string& input,		  T& output, int ErrCode=-1);// return: 1=good 0=bad. ErrCode: >0=elog <0=exit 0=do_nothing
int ReadStr(const std::string& input,	   bool& output, int ErrCode=-1);// output=true if input=1/true/t/yes/y/on/positive (case-insensitive)
int ReadStr(const std::string& input,	 double& output, int ErrCode=-1);// success if read a number and !isnan()
int ReadStr(const std::string& input, IntRanges& output, int ErrCode=-1);// success if read_something && nothing_left
// Beware: If T is double, it reads nan (case insensitive) and return good.

template <typename T>
inline void ReadArg(const std::vector<std::string>& args, size_t& argi, T& value, int ErrCode=-1)
{
	if (argi>=args.size()) exit_error("argi is out of boundary");
	std::string input;
	std::size_t found = args[argi].find('=');
	if (found!=std::string::npos)	{ NeedArg(args,0,argi); input = args[argi].substr(found+1); }
	else							{ NeedArg(args,1,argi); input = args[++argi]; }
	ReadStr(input,value,ErrCode);
}
void ReadArg(const std::vector<std::string>& args, size_t& argi, bool& value, int ErrCode=-1); // read: --opt=Y | --opt Y | --opt

template <typename Container>
inline int ReadSet_add(const std::string& input, Container& values, bool read_empty, std::string d=",:;|\\ \t\n", int ErrCode=-1)
{
	typedef typename Container::value_type T;
	std::string s(input);
	extract_set_add(s,values,read_empty,d);
	if (s.empty())		{ return 1; }
	else if	(ErrCode>0)	{ elog.add(ErrCode); return 0; }
	else if (ErrCode<0)	{ T v; ReadStr(s,v,-1); }
	else				{ return 0; }
	return -1; // never happen
}
/*template <typename Container>
inline int ReadSet(const std::string& input, Container& values, bool read_empty, std::string d=",:;|\\ \t\n", int ErrCode=-1)
{
	values.clear();
	return ReadSet_add(input,values,read_empty,d,ErrCode);
}*/

template <typename Container>
inline void ReadSet(const std::vector<std::string>& args, size_t& argi, Container& values, std::string d=",:;|\\ \t\n", int ErrCode=-1)
{
	std::string input;
	std::size_t found = args[argi].find('=');
	char b4='\0';
	if (found!=std::string::npos)	{ b4 = args[argi][found-1]; input = args[argi].substr(found+1); }
	else							{ NeedArg(args,1,argi);			input = args[++argi]; }
//	if		(b4=='-')	{                 extract_set_remove (input,values);           } // for -=
//	else if	(b4=='&')	{                 extract_set_overlap(input,values);           } // for &=
//	else if	(b4=='*')	{                 extract_set_overlap(input,values);           } // for *=
	if (b4=='+')		{                 ReadSet_add		 (input,values,true, d,ErrCode); } // for +=
	else				{ values.clear(); ReadSet_add		 (input,values,false,d,ErrCode); } // for  =
}

// Output container type C = vector|deque|set|multiset<string>. Return number matched.
template <typename C> int files_matched			(const std::string& regex_str,C& fullpaths,const std::string& dir="");
template <typename C> int files_matched_zipOrNot(const std::string& regex_str,C& fullpaths,const std::string& dir="");

// return 0:not_checked 1:up_to_date -1:not_up_to_date
int chk_net_vers(const std::string& url1, const std::string& url2, const std::string& this_version, std::string& read_version);

#endif
