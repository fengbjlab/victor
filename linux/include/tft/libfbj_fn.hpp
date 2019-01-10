// body file: libfbj_fn.cpp

#ifndef LIB_FN_FILE
#define LIB_FN_FILE

#include <iostream>
#include <string>
#include <vector>
#include <limits>       // std::numeric_limits

/*
FN (Field-number) starts from 1, last FN is -1, second last -2, and so on.
FN 0 denotes the field next to the last field (ie, not exist in input).

To input data: 1) parse_or_exit; 2) push_back
To setup size: 1) size_valid/size_solvable; 2) contents_to_xxx; 3) if (s >= .max_field_num())
To obtain FNs: 1) contents_to_xxx (automatic size check); 2) vector<int>::const_iterator (return 0-based coordinate)
For access by iterator, user is responsible for 1) the validity of size; 2) size_valid() to solve neg. num. 
 
 Characters for the begin of field names instead of field numbers:
 Do not use ()<>, always need quations.
 Do not use {}, because {hello,world} becomes hello and world (two arguments instead of one)
 Do not use [], because my single_condition parser relies on whether [] presents.
 Do not use @,, because I have used it as the begin of a serie.
 Do not use #$`&*?!|\;, because the system uses them
 Do not use _, becaue they may be part of the name
 Do not use ~, because ~user means $HOME of another user
 Possible: % ^ : /
 Do not use ., because . is a valid character for R variable names.
 Do not use :, because I want to reserve it for filenames.
 Do not use ^, because I want to reserve it for NOT.
 Do not use /, because I want to reserve it for OR.
 % looks good, because it looks like a variable name like "date +%m-%d-%Y"
 Need to have the closing character, because the single_condition parser needs to know when it stops. But should also allow not having one.
 */

class bad_query_fieldnumbers_absent  { public: const char *ShowReason() const { return "Querying absent fields."; } };
class bad_query_fieldnumbers_unsolved{ public: const char *ShowReason() const { return "Querying unsolved fields."; } };

class field_numbers {

public:
	// ---------- setup ----------
	~field_numbers();
	field_numbers(bool n,bool a);		// NewEachTime=n AllRequired=a
	field_numbers();					// NewEachTime=false AllRequired=true
	field_numbers(const field_numbers& othr);
	field_numbers& operator=(const field_numbers& orig);
	void clear();
	void set_NewEachTime(bool status);	// Solve negative FN for each row. default=false
	void set_AllRequired(bool status);	// All fields are required. default=true

	// ---------- input data ----------
	void parse_or_exit(std::istream& ss,bool one_input=false);				// Add data aggregately. Exit if nothing parsed.
	std::string parse_or_exit(const std::string& str,bool one_input=false);	// return remaining string
	void push_back(int i);								// Add i-i
	void push_back(int i,int j);						// Add i-j
	std::string print();								// print what has been inputed
	
	// ---------- test ----------
	bool size_valid(int s);				// setup size, return whether s is good (neg. solvable & required exist)
	bool size_solvable(int s);			// setup size, return whether negative numbers are solvable if data_vector_size=s
	int  min_required();				// after size_valid(), return required size if AllRequired=true, and 0 otherwise
	bool no_input();					// has no input yet
	int  num_inputs();					// number of inputs (number of times parse/push_back is called)
	bool is_double();					// has two single field inputs, e.g. 1,2; while 1-2 doesn't count
	bool is_multi();					// has multiple fields
	bool is_multi(int section);			// has multiple fields
	
	// ---------- get field contents, no need to solve size before ----------
	// IF there's an error (query an absent/unsolved field number), throw an exception.
	// To avoid absent field exception: (AllRequired=true && tfile_format::action=Expand/Wr_Cont/Continue/Break/Exit) || AllRequired=false
	// To write elog: AllRequired=true && tfile_format::action=None;
	// To insert empty and write elog: AllRequired=true && tfile_format::action=Expand;
	// To insert default: AllRequired=false && tfile_format::action=None && use the argument def_val;
	// For contents_to_a_string, delimiter is always inserted even when the field is missing
	template<typename T> // delimiter T = string or char
	void contents_to_a_string(const std::vector<std::string>& srce, std::string& dest,  const T& delimiter, bool wr_endl=false, bool quoted=false, const std::string& def_val="");
	template<typename T> // delimiter T = string or char
	void contents_to_ostream (const std::vector<std::string>& srce, std::ostream& dest, const T& delimiter, bool wr_endl=false, bool quoted=false, const std::string& def_val="");
	template<typename C> // Container = vector|deque|set|multiset<double>
	void contents_to_doubles (const std::vector<std::string>& srce, bool clr, C& dest, bool add_def=true, const double def_val=std::numeric_limits<double>::signaling_NaN());
	template<typename C> // Container = vector|deque|set|multiset<string>
	void contents_to_strings (const std::vector<std::string>& srce, bool clr, C& dest, bool add_def=true, const std::string& def_val="");

	// ---------- get field numbers, call only after size is solved, functions check for this except end() rend() ----------
	typedef std::vector<int>::const_iterator const_iterator;

	int									max_field_num();
	bool								contain(const int& n) const;
	std::vector<int>::const_iterator	find(const int& n) const;
	std::vector<int>::size_type			size() const;
	bool								empty() const;
	int operator[]( int pos );			// Beware: return by value, not reference
	std::vector<int>::const_iterator			begin()  const;
	std::vector<int>::const_iterator			end()    const;
	std::vector<int>::const_reverse_iterator	rbegin() const;
	std::vector<int>::const_reverse_iterator	rend()   const;
	std::vector<int>::const_iterator			cbegin() const {return begin();}
	std::vector<int>::const_iterator			cend()   const {return   end();}
	std::vector<int>::const_iterator			begin(int section) const;
	std::vector<int>::const_iterator			end(int section) const;
	int& front();
	int& front(int section);
	
	// ---------- operators between two field_numbers, call only after size is solved, functions check for this ----------
	friend field_numbers operator&(const field_numbers& left, const field_numbers& right);  // previously *
	friend field_numbers operator|(const field_numbers& left, const field_numbers& right);
	friend field_numbers operator+(const field_numbers& left, const field_numbers& right);
	friend field_numbers operator-(const field_numbers& left, const field_numbers& right);

private:
	struct field_numbers_data;
	field_numbers_data * d;
};

int  ReadFld(const std::string& input,	field_numbers& f, bool one_input=false, int ErrCode=-1);
void ReadArg(const std::vector<std::string>& args, size_t& argi, field_numbers& f, bool one_input=false, int ErrCode=-1);

#endif

