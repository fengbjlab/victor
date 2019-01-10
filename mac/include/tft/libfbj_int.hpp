// body file: libfbj_int.cpp

#ifndef LIBFBJ_IntRanges
#define LIBFBJ_IntRanges

#include <iostream>

// This class has three modes of input: 1) num_pairs; 2) {math_expression}; 3) :file: w/ num_pairs.
// If math expression is used and has infinity number, results must be monotonously increasing !!!
// next() has a check point for this, and parse() check for the first 100 iterations. However,
// it still could be wrong if result decrease after 100 it. min() looks for the first number only,
// and max() return INT_MAX if math_expr has infinity, which all assume monotonic increment.
// The reason I use {} and :: as signals for math_expr and file, respectively, is
// 1) [] is reserved for conditions, and unix command line doesn't take () literally
// 2) : is the only character after which one can use TAB to auto-fill a filename in the command line

class IntRanges {

public:
	// initialize
	~IntRanges();
	IntRanges();
	IntRanges(bool omit_inverse_input);	// default=false
	IntRanges(const IntRanges& othr);
	IntRanges&	operator=(const IntRanges& orig);
	void		operator=(const std::string& in);
	void set_math_signal(char c='{');	// default is {
	void set_file_signal(char c=':');	// default is :
	void clear();						// clear all data
	
	// setup data
	void push_back(int n);				// push_back n  to vector
	void push_back_n_to_max(int n);		// push_back n- to vector
	int  parse(std::string& str);		// return whether sth read, str become the rest
	int  parse(std::istream& es);		// return whether sth read
	
	// access data
	bool empty();						// no input yet
	bool has_infinity();				// have infinity like 1- or {[1-]<math_expr>}
	int  max();							// max possible number, return INT_MAX if is math_expr & has infinity
	int  min();							// min possible number, wrong if math_expr has_inf & not mono increasing
	bool include(int n);				// whether n is included
	bool is_n_to_max(int n);			// whether is [n,INT_MAX]
	
	// output
	void print(std::ostream& os) const;	// iterator will be missing if in math_expr
	friend std::ostream& operator<<(std::ostream& os, const IntRanges& ir);
	
private:
	struct IntRangesData;
	IntRangesData *idta;
};

inline std::ostream& operator<<(std::ostream& os, const IntRanges& ir)
{
	ir.print(os);
    return os;
}

#endif
