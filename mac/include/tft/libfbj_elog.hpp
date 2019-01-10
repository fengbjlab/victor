// body file: libfbj_elog.cpp

#ifndef LIBELOG_FILE
#define LIBELOG_FILE

#include <iostream>
#include <string>

class ErrorLogger {

public:
	~ErrorLogger();
	ErrorLogger();
	ErrorLogger(const ErrorLogger& othr);				// default output is cerr
	ErrorLogger& operator=(const ErrorLogger& orig);
	void set_output(std::ostream& ostream_reference);	// you can do this: elog.set_output(cerr);
	void set_output(std::ostream* ostream_pointer);		// you can do this: elog.set_output(NULL);
	int  get_token(const std::string& error_message, bool to_show_error_count=true, bool to_show_message=true);
	void rewind(int token);
	void add(int token);
	void add(int token, const std::string& associated_str);
	int  add(const std::string& error_message, bool to_show_error_count=true);
	int  add(const std::string& error_message, const std::string& associated_str, bool to_show_error_count=true);
	void write(std::ostream& ostream_reference);
	void write(std::ostream& ostream_reference, int token);
	void write(int token);								// write to default output
	void write();										// write to default output
	void set_lock(const std::string& filename);
	bool lock_set();

private:
	struct ErrorLoggerData;
	ErrorLoggerData *d;
};

#endif

/* Error logging class - log errors with names, counts and associated texts, and print errors before exit (Default=cerr).
Usage 1: elog.add(name,[associated text],[to show counts]);
Usage 2: int token=get_token(name,toShow); elog.add(token); // faster <= doesn't check whether the token exists.
Token starts from 1, so that you can "if (token) .." and elog.add(0) does nothing quietly.
*/

