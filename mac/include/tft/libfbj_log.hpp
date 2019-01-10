// body file: libfbj_log.cpp
/* Hierarchical output onto screen and/or into a log file, count number of errors & warnings, log begin/end time of running.
 Category of message: line=line, err=error, wrn=warning.
 Category of output:  show=on_screen, write=on_screen+to_file
 Place the message type at the beginning, then the content, then optionally flush_logger.
   By lns<<"your message", logger simply stores the message without showing it, unless there's a flush_logger at the end.
   At showl/writel/..., logger will automatically flush the stored message and print a new line with a correct prefix (some spaces then Wrn/Err/None).
 To suppress showing on screen, use program.quiet
 To suppress writing to a file, do not use logger::open  */

#ifndef LIBFBJ_LOG
#define LIBFBJ_LOG

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#define _NUM_LOGGER_CAT 3

class logger : public std::stringstream 
{
public:
	enum counts_category  {err,wrn,line};
	logger();								// log file not opened
	logger(const std::string& filename);	// open log file
	logger(const logger& othr);				// forbidden operation
	logger& operator=(const logger& orig);	// forbidden operation
	~logger();								// close file
	void		open(std::string logfilename, const std::ios_base::openmode& mode);
	void		close(int is_fatal);		// close log file
	void		flush_log();
	void		test();						// only for debugging
	void		fatal();
	void		sub();						
	void		endsub();
	void		show(int type);
	void		write(int type);
	void		continue_show(int type);
	void		continue_write(int type);
	std::string logger_input_string();
	double		logger_input_double();
	char		logger_input_char();
	void		rewind_timer();
	void		set_lock(const std::string& filename);
	bool		lock_set();
private:
	struct LoggerData;
	LoggerData * d;
};

struct flush_logger_
{
	inline friend logger& operator<< (std::ostream& os, const flush_logger_& w)
	{
		dynamic_cast<logger&>(os).flush_log();
		return dynamic_cast<logger&>(os);
	}
} ;

struct fatal_
{
	inline friend logger& operator<< (std::ostream& os, const fatal_& w)
	{
		dynamic_cast<logger&>(os).fatal();
		return dynamic_cast<logger&>(os);
	}
} ;

struct showline_
{
	inline friend logger& operator<< (std::ostream& os, const showline_& w)
	{
		dynamic_cast<logger&>(os).show(logger::line);
		return dynamic_cast<logger&>(os);
	}
} ;

struct showerr_
{
	inline friend logger& operator<< (std::ostream& os, const showerr_& w)
	{
		dynamic_cast<logger&>(os).show(logger::err);
		return dynamic_cast<logger&>(os);
	}
} ;

struct showwarning_
{
	inline friend logger& operator<< (std::ostream& os, const showwarning_& w)
	{
		dynamic_cast<logger&>(os).show(logger::wrn);
		return dynamic_cast<logger&>(os);
	}
} ;

struct cshowline_
{
	inline friend logger& operator<< (std::ostream& os, const cshowline_& w)
	{
		dynamic_cast<logger&>(os).continue_show(logger::line);
		return dynamic_cast<logger&>(os);
	}
} ;

struct writeline_
{
	inline friend logger& operator<< (std::ostream& os, const writeline_& w)
	{
		dynamic_cast<logger&>(os).write(logger::line);
		return dynamic_cast<logger&>(os);
	}
} ;

struct writeerr_
{
	inline friend logger& operator<< (std::ostream& os, const writeerr_& w)
	{
		dynamic_cast<logger&>(os).write(logger::err);
		return dynamic_cast<logger&>(os);
	}
} ;

struct writewarning_
{
	inline friend logger& operator<< (std::ostream& os, const writewarning_& w)
	{
		dynamic_cast<logger&>(os).write(logger::wrn);
		return dynamic_cast<logger&>(os);
	}
} ;

extern flush_logger_	flush_logger;
extern fatal_			fatal;
extern showline_		showl;
extern showerr_			showe;
extern showwarning_		showw;
extern cshowline_		cshowl;
extern writeline_		writel;
extern writeerr_		writee;
extern writewarning_	writew;

#endif
