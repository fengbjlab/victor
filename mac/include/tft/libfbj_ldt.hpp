// body file: libfbj_ldt.cpp

#ifndef LIBLDT_FILE
#define LIBLDT_FILE

#include <iostream>
#include <vector>
#include <string>

/* possible problem:
1) single condition - no flag to show whether it has successfully read a condition before testing
*/

extern std::string cond_row_selection_help_text;

class single_condition {

private:
	struct SingleConditionData;
	SingleConditionData *d;
	
public:
	~single_condition();
	single_condition();
	single_condition(const single_condition& othr);
	single_condition& operator=(const single_condition& orig);
	void initialize();
	void suppress_elog();
	char type_of_test();
	std::string ref_string();
	int first_field();
	int  max_field_num();
	bool size_valid(int size);
	std::pair<bool,bool> evaluate(const std::vector<std::string>& sv);
	std::pair<bool,bool> eval_without_inverse(const std::vector<std::string>& sv);
	bool read (std::istream & in);
	void write_eval(std::ostream& out);
	std::string help_text() {return cond_row_selection_help_text; }
};

class multiple_conditions {

private:
	struct MultipleConditionData;
	MultipleConditionData *d;
	
public:
	~multiple_conditions();
	multiple_conditions();
	multiple_conditions(const multiple_conditions& othr);
	multiple_conditions& operator=(const multiple_conditions& orig);
	void read_arguments(std::vector<std::string>& srce_opt);
	void read(const std::string& s);
	std::pair<bool,bool> evaluate(std::vector<std::string>& v);
	bool is_defined();
	bool is_valid();
	bool satisfied();
	int  group_number();
	int  max_field_num();
	bool size_valid(int size);
	std::string help_text() {return cond_row_selection_help_text; }
//	bool is_true();
//	std::pair<bool,bool> result();
};

int  ReadCnd(const std::string& input,	multiple_conditions& c, int ErrCode=-1);
void ReadArg(const std::vector<std::string>& args, size_t& argi, multiple_conditions& c, int ErrCode=-1);

#endif
