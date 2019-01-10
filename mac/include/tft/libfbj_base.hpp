// body file: libfbj_base.cpp

#ifndef LIBFBJ_BASE
#define LIBFBJ_BASE

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/algorithm/string.hpp>						// for boost::split 

//---------------- container-related macros --------------------

#define exist_element(c,e)		((c).find((e))!=(c).end())

#define DefTermWidth 150

template <typename Container1, typename value_type>
bool exist(Container1& c, value_type& e)
{
	for (auto &i:c)
		if (i==e) return true;
	return false;
}

template <typename Container1, typename Container2>
bool exist_any(Container1& c, Container2& e)	// c and e are interchangable
{
	for (auto &i:e)
		if (exist_element(c,i)) return true;
	return false;
}

template <typename Container1, typename Container2>
bool exist_all(Container1& c, Container2& e)
{
	for (auto &i:e)
		if (!exist_element(c,i)) return false;
	return true;
}

template <typename Container1, typename Container2>
bool exist_any_lower(Container1& c, Container2& e) // container c elements are all lower-case, e is case-insensitive
{
	for (auto &i:e)
		if (exist_element(c,boost::to_lower_copy(i))) return true;
	return false;
}

template <typename Container1, typename Container2>
bool exist_all_lower(Container1& c, Container2& e) // container c elements are all lower-case, e is case-insensitive
{
	for (auto &i:e)
		if (!exist_element(c,boost::to_lower_copy(i))) return false;
	return true;
}

template <typename Container1, typename Container2>
bool exist_any_upper(Container1& c, Container2& e) // container c elements are all upper-case, e is case-insensitive
{
	for (auto &i:e)
		if (exist_element(c,boost::to_upper_copy(i))) return true;
	return false;
}

template <typename Container1, typename Container2>
bool exist_all_upper(Container1& c, Container2& e) // container c elements are all upper-case, e is case-insensitive
{
	for (auto &i:e)
		if (!exist_element(c,boost::to_upper_copy(i))) return false;
	return true;
}

template<class T> T const& constant(T& v){ return v; }
#if __cplusplus >= 201103L
	#define each_element(c,i)			auto i(        (c).begin() ); i!=(c).end(); ++i
	#define each_element_rev(c,i)		auto i(        (c).rbegin()); i!=(c).rend();++i
	#define each_element_const(c,i)		auto i(constant(c).begin() ); i!=(c).end(); ++i
	#define each_interval(c,i1,i2)		auto i2=(c).begin(), i1=i2++; i2!=(c).end();  ++i1,++i2
	#define each_interval_rev(c,i1,i2)	auto i2=(c).rbegin(),i1=i2++; i2!=(c).rend(); ++i1,++i2
#elif defined __GNUG__
	#define each_element(c,i)		typeof((c).begin())                              i=(c).begin(); i!=(c).end(); ++i
	#define each_element_rev(c,i)	typeof((c).rbegin())                             i=(c).rbegin();i!=(c).rend();++i
	#define each_element_const(c,i)	typeof(const_cast<const typeof(c)&>(c).begin())  i=(c).begin(); i!=(c).end(); ++i
#else
	#error This library requires either a GNU C++ compiler or C++0x/C++11.
#endif

// c++11 has the new "for" syntax: for (auto &x:c) / for (auto x:c). x is not an iterator!
// each_element_const is different from both in that the compiler complains when the user try to modify content.

// --------------- print container data -----------------

// Container should have empty() begin() end() 
template <typename Container, typename DelType>
inline void print_container(const Container& val, std::ostream& file, const DelType& del, bool write_endl=false, bool quoted=false)
{
	if (!val.empty())
	{
		typename Container::const_iterator it(val.begin());
		if (quoted)
		{
			file << '\"' << *it << '\"' ;
			for (++it; it!=val.end(); ++it) file << del << '\"' << *it << '\"' ;
		}
		else
		{
			file << *it;
			for (++it; it!=val.end(); ++it) file << del << *it;
		}
	}
	if (write_endl) file<<std::endl;
}

template <typename Container, typename DelType>
inline void print_container_head(const Container& val, int num, std::ostream& file, const DelType& del, bool write_endl=false, bool quoted=false)
{
	if (!val.empty())
	{
		typename Container::const_iterator it(val.begin());
		int i(0);
		if (quoted)
		{
			for (; i<num && it!=val.end(); ++i,++it)
			{
				if (i) file << del;
				file << '\"' << *it << '\"' ;
			}
		}
		else
		{
			for (; i<num && it!=val.end(); ++i,++it)
			{
				if (i) file << del;
				file << *it;
			}
		}
	}
	if (write_endl) file<<std::endl;
}

template <typename Container, typename DelType>
inline void print_container_tail(const Container& val, int num, std::ostream& file, const DelType& del, bool write_endl=false, bool quoted=false)
{
	if (!val.empty())
	{
		typename Container::const_iterator it;
		if ((int)val.size()>num) {	it=val.end(); it-=num; }
		else						it=val.begin();
		if (quoted)
		{
			file << '\"' << *it << '\"' ;
			for (++it; it!=val.end(); ++it) file << del << '\"' << *it << '\"' ;
		}
		else
		{
			file << *it;
			for (++it; it!=val.end(); ++it) file << del << *it;
		}
	}
	if (write_endl) file<<std::endl;
}

template <typename Container, typename DelType>
inline void print_container_after(const Container& val, int num, std::ostream& file, const DelType& del, bool write_endl=false, bool quoted=false)
{
	if (!val.empty())
	{
		typename Container::const_iterator it(val.begin());
		for (int i=0; i<=num && it!=val.end(); ++i,++it) ;
		if (quoted)
		{
			for (int i=0; it!=val.end(); ++i,++it)
			{
				if (i) file << del;
				file << '\"' << *it << '\"' ;
			}
		}
		else
		{
			for (int i=0; it!=val.end(); ++i,++it)
			{
				if (i) file << del;
				file << *it;
			}
		}
	}
	if (write_endl) file<<std::endl;
}

template <typename T1, typename T2, typename DelType>
inline void print_map_first(const std::map<T1,T2>& val, std::ostream& file, const DelType& del, bool write_endl=false)
{
	if (!val.empty())
	{
		typename std::map<T1,T2>::const_iterator it(val.begin());
		file << it->first;
		for (++it; it!=val.end(); ++it) file << del << it->first;
	}
	if (write_endl) file<<std::endl;
}

template <typename T1, typename T2, typename DelType>
inline void print_map_second(const std::map<T1,T2>& val, std::ostream& file, const DelType& del, bool write_endl=false)
{
	if (!val.empty())
	{
		typename std::map<T1,T2>::const_iterator it(val.begin());
		file << it->second;
		for (++it; it!=val.end(); ++it) file << del << it->second;
	}
	if (write_endl) file<<std::endl;
}

template <typename T1, typename T2, typename DelType>
inline void print_map(const std::map<T1,T2>& val, std::ostream& file, const DelType& del_i, const DelType& del_e, bool write_endl=false)
{
	if (!val.empty())
	{
		typename std::map<T1,T2>::const_iterator it(val.begin());
		file << it->first << del_i << it->second;
		for (++it; it!=val.end(); ++it) file << del_e << it->first << del_i << it->second;
	}
	if (write_endl) file<<std::endl;
}

template <typename Container, typename DelType> // container item type must be string
inline void print_uniq(const Container& val, std::ostream& file, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::string previous="dummy_string_that_will_not_appear";
	bool is_first=true;
	for (auto &x:val)
	{
		if (x!=previous)
		{
			previous=x;
			if (is_first)	{ if (quoted) file<<'\"'<<x<<'\"'; 		else file<<x; 		is_first=false; }
			else			{ if (quoted) file<<del<<'\"'<<x<<'\"'; else file<<del<<x; }
		}
	}
	if (write_endl) file<<std::endl;
}

// --------------- print container to string ----------------

template <typename Container, typename DelType>
inline void container_to_str(const Container& val, std::string& str, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container(val,ostr,del,write_endl,quoted);
	str = ostr.str();
}

template <typename Container, typename DelType>
inline void container_head_to_str(const Container& val, int num, std::string& str, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_head(val,num,ostr,del,write_endl,quoted);
	str = ostr.str();
}

template <typename Container, typename DelType>
inline void container_tail_to_str(const Container& val, int num, std::string& str, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_tail(val,num,ostr,del,write_endl,quoted);
	str = ostr.str();
}

template <typename Container, typename DelType>
inline void container_after_to_str(const Container& val, int num, std::string& str, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_after(val,num,ostr,del,write_endl,quoted);
	str = ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline void map_first_to_str(const std::map<T1,T2>& val, std::string& str, const DelType& del, bool write_endl=false)
{
	std::stringstream ostr;
	print_map_first(val,ostr,del,write_endl);
	str = ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline void map_second_to_str(const std::map<T1,T2>& val, std::string& str, const DelType& del, bool write_endl=false)
{
	std::stringstream ostr;
	print_map_second(val,ostr,del,write_endl);
	str = ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline void map_to_str(const std::map<T1,T2>& val, std::string& str, const DelType& del_i, const DelType& del_e, bool write_endl=false)
{
	std::stringstream ostr;
	print_map(val,ostr,del_i,del_e,write_endl);
	str = ostr.str();
}

// --------------- print container to return string -----------------

template <typename Container, typename DelType>
inline std::string str_of_container(const Container& val, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container(val,ostr,del,write_endl,quoted);
	return ostr.str();
}

template <typename Container, typename DelType>
inline std::string str_of_container_head(const Container& val, int num, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_head(val,num,ostr,del,write_endl,quoted);
	return ostr.str();
}

template <typename Container, typename DelType>
inline std::string str_of_container_tail(const Container& val, int num, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_tail(val,num,ostr,del,write_endl,quoted);
	return ostr.str();
}

template <typename Container, typename DelType>
inline std::string str_of_container_after(const Container& val, int num, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_container_after(val,num,ostr,del,write_endl,quoted);
	return ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline std::string str_of_map_first(const std::map<T1,T2>& val, const DelType& del, bool write_endl=false)
{
	std::stringstream ostr;
	print_map_first(val,ostr,del,write_endl);
	return ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline std::string str_of_map_second(const std::map<T1,T2>& val, const DelType& del, bool write_endl=false)
{
	std::stringstream ostr;
	print_map_second(val,ostr,del,write_endl);
	return ostr.str();
}

template <typename T1, typename T2, typename DelType>
inline std::string str_of_map(const std::map<T1,T2>& val, const DelType& del_i, const DelType& del_e, bool write_endl=false)
{
	std::stringstream ostr;
	print_map(val,ostr,del_i,del_e,write_endl);
	return ostr.str();
}

template <typename Container, typename DelType>
inline std::string str_of_uniq(const Container& val, const DelType& del, bool write_endl=false, bool quoted=false)
{
	std::stringstream ostr;
	print_uniq(val,ostr,del,write_endl,quoted);
	return ostr.str();
}

// --------------- vec_deq functions -----------------
// These functions have the names vec_deq_xxx(). They are templated functions disgned for vectors or deques.
// The aim is to access/modify vector[i]/deque[i] even when "i" is bebyond the container's size().
// Negative i is deemed illegal, so no checkpoint is made and let the program crash.

template <typename Container>
inline typename Container::value_type vec_deq_get(const Container& vd, const int loc) // prv: get_vec_val
{
	typedef typename Container::value_type T;
	if (loc >= (int)vd.size())	return	T();
	else						return	vd[loc];
}

template <typename Container>
inline typename Container::value_type vec_deq_get(const Container& vd, const int loc, const typename Container::value_type & def) // prv: get_vec_val
{
	if (loc >= (int)vd.size())	return	def;
	else						return	vd[loc];
}

template <typename Container>
inline void vec_deq_set(Container& vd, const int loc, const typename Container::value_type & val)  // prv: set_vec_val
{
	typedef typename Container::value_type T;
	int i=loc-vd.size();
	if (i>=0) { vd.insert(vd.end(),i,T()); vd.push_back(val); }
	else vd[loc]=val;
}

template <typename Container>
inline void vec_deq_set(Container& vd, const int loc, const typename Container::value_type & val, const typename Container::value_type & def)  // prv: set_vec_val
{
	int i=loc-vd.size();
	if (i>=0) { vd.insert(vd.end(),i,def); vd.push_back(val); }
	else vd[loc]=val;
}

template <typename Container>
inline void vec_deq_ins(Container& vd, const int loc, const typename Container::value_type & val, const int num=1) // prv: insert_vec_val
{
	typedef typename Container::value_type T;
	int i=loc-vd.size();
	if (i>0) { vd.insert(vd.end(),i,T()); vd.insert(vd.end(),num,val); }
	else vd.insert(vd.begin()+loc,num,val);
}

template <typename Container>
inline void vec_deq_ins(Container& vd, const int loc, const typename Container::value_type & val, const int num, const typename Container::value_type & def) // prv: insert_vec_val
{
	int i=loc-vd.size();
	if (i>0) { vd.insert(vd.end(),i,def); vd.insert(vd.end(),num,val); }
	else vd.insert(vd.begin()+loc,num,val);
}

// --------- template that can use a std::vector or a std::set -----------
// http://stackoverflow.com/questions/15574499/writing-a-template-that-can-use-a-stdvector-or-a-stdset

#if __cplusplus >= 201103L
#include <utility>
template<typename C, typename T>
auto add_to_container(C& c, T&& t) ->
decltype(c.push_back(std::forward<T>(t)), void())
{
    c.push_back(std::forward<T>(t));
}

template<typename C, typename T>
auto add_to_container(C& c, T&& t) ->
decltype(c.insert(std::forward<T>(t)), void())
{
    c.insert(std::forward<T>(t));
}
#endif

//------------ system environment variables and commands --------------

struct EnvironmentVariable {
	bool exist;
	bool is_true;
	std::string value;
	EnvironmentVariable():exist(false),is_true(false){}
};

EnvironmentVariable		bool_env(const std::string& name);
EnvironmentVariable		str_env(const std::string& name);
void		get_term_dim(int& lines, int& columns); // columns will be real columns -1 so that emacs works well
std::string HOME(); // Always ends with /
std::string linux_command_which_dir(const std::string& prg, const std::string& var="PATH", bool check_exe=true);
// Different from linux "which" that returns /path/to/prg, this function returns /path/to/ or an empty string
std::string to_term(std::string msg);

//---------------- basic string functions --------------------

bool str_has		(const std::string& s, const std::string& subs);
bool str_startsw	(const std::string& s, const std::string& subs);
bool str_endsw		(const std::string& s, const std::string& subs);
bool is_numeric		(const std::string& s); // can be cast to double
bool is_integer		(const std::string& s); // can be cast to int
bool is_white_spaces(const std::string& s); // TABs or SPACEs, not counting \n \r
bool is_a_valid_name(const std::string& s); // contains [_a-zA-Z0-9], start w/ [_a-zA-Z]
bool IsBool			(const std::string& s); // 1,true,t,yes,y | 0,false,f,no,n (case-insensitive)
bool IsYes			(const std::string& s); // throw input_exception if not boolean
bool IsNo			(const std::string& s); // no throw input_exception if not boolean
inline std::string s(const char * cs)	{ return std::string(cs);  }
inline std::string s(const char c)		{ return std::string(1,c); }

// Even NaN or nan will be read successful: v will become nan and the function return true.
// If nan is not wanted, then it's better to use the read_val_xxx versions or write a checkpoint after calling this function.
template <typename T>
bool read_val(const std::string& s, T& v)
{
	try { v=boost::lexical_cast<T>(s); }
	catch (boost::bad_lexical_cast &) { return false; }
	return true;
}

inline bool read_val_noNaN(const std::string& s, double& t)
{
	double v;
	try { v=boost::lexical_cast<double>(s); }
	catch (boost::bad_lexical_cast &) { return false; }
	if (std::isnan(v)) return false;
	t=v;
	return true;
}

template <typename T> bool read_val_ge(const std::string& s, T& v, T threshold) { T t; if (!read_val(s,t)) return false; if (t>=threshold) { v=t; return true; } return false; }
template <typename T> bool read_val_gt(const std::string& s, T& v, T threshold) { T t; if (!read_val(s,t)) return false; if (t> threshold) { v=t; return true; } return false; }
template <typename T> bool read_val_le(const std::string& s, T& v, T threshold) { T t; if (!read_val(s,t)) return false; if (t<=threshold) { v=t; return true; } return false; }
template <typename T> bool read_val_lt(const std::string& s, T& v, T threshold) { T t; if (!read_val(s,t)) return false; if (t< threshold) { v=t; return true; } return false; }
template <typename T> bool read_val_ge_le(const std::string& s, T& v, T l, T r) { T t; if (!read_val(s,t)) return false; if (t>=l && t<=r) { v=t; return true; } return false; }
template <typename T> bool read_val_gt_lt(const std::string& s, T& v, T l, T r) { T t; if (!read_val(s,t)) return false; if (t> l && t< r) { v=t; return true; } return false; }
template <typename T> bool read_val_ge_lt(const std::string& s, T& v, T l, T r) { T t; if (!read_val(s,t)) return false; if (t>=l && t< r) { v=t; return true; } return false; }
template <typename T> bool read_val_gt_le(const std::string& s, T& v, T l, T r) { T t; if (!read_val(s,t)) return false; if (t> l && t<=r) { v=t; return true; } return false; }

double as_double_or_nan(const std::string& s);
int as_int_or_max(const std::string& s);
int as_int_or_min(const std::string& s);

// ---------- sub-string extraction from the beginning of a string -------------------

class bad_extracting { public:
	bad_extracting() {};
	~bad_extracting() {};
	const char *ShowReason() const { return "Failed to extract the right type of data from a string."; }
};

char		extract_char	 (std::string& s); // throw bad_extracting if nothing extracted
double		extract_double	 (std::string& s); // throw bad_extracting if nothing extracted
int			extract_int		 (std::string& s); // throw bad_extracting if nothing extracted
std::string extract_alphabets(std::string& s);
std::string extract_alnum	 (std::string& s);
std::string extract_name	 (std::string& s);
template <typename chr_t> bool extract_char  (std::string& s, chr_t& res) { try { res=extract_char(s);   return true; } catch (bad_extracting &e) { return false; } }
template <typename dbl_t> bool extract_double(std::string& s, dbl_t& res) { try { res=extract_double(s); return true; } catch (bad_extracting &e) { return false; } }
template <typename int_t> bool extract_int   (std::string& s, int_t& res) { try { res=extract_int(s);    return true; } catch (bad_extracting &e) { return false; } }

#if __cplusplus >= 201103L
// previous version didn't have read_empty, which equal to read_empty=false.
template <typename Container> // c not cleared before adding, s becomes remains like other extract_xxx()
inline void extract_set_add	 (std::string& s, Container& c, bool read_empty=false, std::string d=",:;|\\ \t\n") // no _ (part of a name) no / (part of a path)
{
	typedef typename Container::value_type ElementT;
	std::vector<std::string> strs;
	std::string remains;
	boost::split(strs, s, boost::is_any_of(d));
	if (!read_empty && strs.size()==1 && strs[0].empty()) strs.clear(); // empty is not read
	for (each_element(strs,it))
	{
		if (*it=="EmptyString") it->clear();
		try	{ add_to_container ( c, boost::lexical_cast<ElementT>(*it) ); }
		catch (boost::bad_lexical_cast &) { remains+=*it; remains+=d[0]; }
	}
	s=remains;
}/*
template <typename Container>
inline void extract_set_remove(std::string& s, Container& c, std::string d=",:;|\\ \t\n") // no _ (part of a name) no / (part of a path)
{
	typedef typename Container::value_type ElementT;
	std::set<std::string> strs;	boost::split(strs, s, boost::is_any_of(d));
	std::set<ElementT>    vals;	for (auto &i:strs) vals.insert(boost::lexical_cast<ElementT>(i));
	for (bool removed=true; removed;)
	{
		removed=false;
		typename Container::const_iterator it;
		for (it=c.begin(); it!=c.end(); ++it)
			if (exist_element(vals, *it)) { c.erase(it); removed=true; break; }
	}
//	for (auto &i : vals) c.erase( std::remove (c.begin(), c.end(), i), c.end() ); // doesn't work with Values<..>
}
template <typename Container>
inline void extract_set_overlap(std::string& s, Container& c, std::string d=",:;|\\ \t\n") // no _ (part of a name) no / (part of a path)
{
	typedef typename Container::value_type ElementT;
	std::set<std::string> strs;	boost::split(strs, s, boost::is_any_of(d));
	std::set<ElementT>    vals;	for (auto &i:strs) vals.insert(boost::lexical_cast<ElementT>(i));
	for (bool removed=true; removed;)
	{
		removed=false;
		typename Container::const_iterator it;
		for (it=c.begin(); it!=c.end(); ++it)
			if (!exist_element(vals, *it)) { c.erase(it); removed=true; break; }
	}
	//	for (auto &i : vals) c.erase( std::remove (c.begin(), c.end(), i), c.end() ); // doesn't work with Values<..>
}*/
template <typename Container> // c not cleared before adding, s becomes remains like other extract_xxx()
inline void extract_set		(std::string& s, Container& c, std::string d=",:;|\\ \t\n") // no _ (part of a name) no / (part of a path)
{
	c.clear();
	extract_set_add(s,c,false,d);
}
#else
  #error This library requires C++11.
#endif

// INFO has the format var=value
bool		get_int   (const std::vector<std::string>& INFO, const std::string& var, int& dest);	// return true if success
double		get_value (const std::vector<std::string>& INFO, const std::string& var);				// return nan if NaN / !exist
std::string get_string(const std::vector<std::string>& INFO, const std::string& var);				// return empty if not exist

bool		get_int_sw   (const std::vector<std::string>& INFO, const std::string& var, int& dest);	// return true if success
double		get_value_sw (const std::vector<std::string>& INFO, const std::string& var);			// return nan if NaN / !exist
std::string get_string_sw(const std::vector<std::string>& INFO, const std::string& var);			// return empty if not exist

void        put_string(      std::vector<std::string>& INFO, const std::string& var, const std::string& val); // return 1=add 0=update

//---------------- string conversion --------------------

int int_width(long long num,int base=10);	// number of digits of an integer
int sig_digits(const std::string& input);	// number of significant digits of a double, allow +/- sign and scientific annotation
char num2char(int num);	// c is [0-9a-z], num is [0-35], exit_error is out of range
int char2num(char c);	// c is [0-9a-z], num is [0-35], exit_error is out of range
int read_int(const std::string& str, std::string::size_type& loc);
char *      itoa (long long num, char* result, int base );
std::string itos (long long num);
char *      itoa_format (long long num, char align, int width, char fill, int base, char * output);
std::string itos_format (long long num, char align, int width, char fill, int base);
std::string itos_max	(long long num, long long maxnum);
std::string itos_26based(long long num);

std::string ftos(double value,const std::string& MissingStr="nan");
std::string ftos(double value, int num_sig_digits, const std::string& MissingStr="nan");
std::string ftos_FixWidth(double value, int width_limits=8, const std::string& MissingStr="nan");
std::string ftos_MaxWidth(double value, int width_limits=8, const std::string& MissingStr="nan");
std::string ftos_MaxDigit(double value, int numbr_digits=3, const std::string& MissingStr="nan");
inline std::string ftos_FixWidthNoTiny(double value, int width_limits=8, const std::string& MissingStr="nan") { std::string r=ftos_FixWidth(value,width_limits,MissingStr); if (str_has(r,"e-")) return "0"; else return r; }
inline std::string ftos_MaxWidthNoTiny(double value, int width_limits=8, const std::string& MissingStr="nan") { std::string r=ftos_MaxWidth(value,width_limits,MissingStr); if (str_has(r,"e-")) return "0"; else return r; }
// Above MissingStr doesn't affect inf.
//						f=0.1234567879	f=0.1		nan		[-]inf
// ftos(f)				0.123457		0.1			nan		[-]inf
// ftos(f,3)			0.123			0.100		nan		[-]inf
// ftos_FixWidth(f)		0.12346			0.10000		nan		[-]inf
// ftos_MaxWidth(f)		0.12346			0.1			nan		[-]inf
// width_limits=8 is recommended. If width_limits=7, the difference may be up to 5% depending on the range of the numbers.
// It's always a good to try and find the minimum width_limits before transforming data for each data set.

// several ways to do tolower/toupper
// 1) int toupper(const int character);
// 2) std::transform(str.begin(), str.end(), str2.begin(), std::tolower );
// 3) std::for_each(str.begin(), str.end(), std::tolower );
// 4) void boost::to_lower(string&); string boost::to_lower_copy(const string&) <boost/algorithm/string.hpp>
// 5) The following functions, which can be used as to_upper_copy(char*)
inline std::string to_upper_copy(const std::string str) { return boost::to_upper_copy(str); }
inline std::string to_lower_copy(const std::string str) { return boost::to_lower_copy(str); }
inline int convert_case(int c) {
	if		(isupper(c)) return tolower(c);
	else if (islower(c)) return toupper(c);
	else return c;
}

// This funciton returns a string that will NOT contain any 'from' sub-string,
// e.g., replacing // to /, if context=/// it returns /.
// This is different from boost::algorithm::replace_all, which returns //.
void		replace_all_R(			 std::string& context, const std::string& from, const std::string& to);
std::string replace_all_R_copy(const std::string& context, const std::string& from, const std::string& to);

// Similar to unix::echo except that it doesn't have \c \e \E and change \0num
std::string replace_escape_sequence_copy(const std::string& s);
std::string reverse_escape_sequence_copy(std::string s);// TAB to \t, etc.
std::string reverse_escape_char(const char c);			// TAB to \t, etc.

// input => one command argument. Therefore, double quote if there is any white space, escape sequence, special characters [#<>|;&], or starts with !.
// wild cards [*?] and environment variables [$] should be determined at runtime, so I don't quote them.
std::string to_cmd(const std::string& input);			// TAB to \t, etc. && double quote if contain SPACE / # / ESCAPE_SEQUENCEs

//---------------- string trimming ----------------

// Input is not a reference, because otherwise char[] is impossible. So the input will not be changed.
// If "sch" is not found, return the whole string.
std::string substr_after_find(std::string str, const std::string& sch, std::string::size_type index=0);
std::string substr_before_find(std::string str, const std::string& sch, std::string::size_type index=0);
std::string trim_before_find(std::string str, const std::string& sch, std::string::size_type index=0);
std::string trim_after_find(std::string str, const std::string& sch, std::string::size_type index=0);
std::string substr_after_rfind(std::string str, const std::string& sch);
std::string substr_before_rfind(std::string str, const std::string& sch);
std::string trim_before_rfind(std::string str, const std::string& sch);
std::string trim_after_rfind(std::string str, const std::string& sch);

//---------- file system -------------

void filename_change_home_path(std::string& filename); // replace ~ $HOME and other $EnvVar (even EnvVar does not exist)
bool FileExists(std::string strFilename);
bool DirExists(std::string strFilename);
void mkdir(const std::string& pathname); // pathname better not end with / but it's allowed.
bool find_file_zippedOrNot(std::string filename,std::string& trueName);
bool find_file_zippedOrNot(const std::string& filename);
bool is_stdin(const std::string& filename);
std::string label_stdin();
std::string self_path();
inline std::string self_name() { return substr_after_rfind(self_path(),"/"); }
std::string exec(const std::string& cmd, bool nothrow);

void exit_error (const std::string& s);
//inline void exit_error (const std::string& s)	{ std::cerr<<std::endl<<"###!!!### "<<self_name()<<" exit with an error: "<<s<<std::endl; exit(1); }
inline void exit_lackMemory(const std::string& s)	{ exit_error("Not enough memory for "+s);	}
inline void exit_cannotOpen(const std::string& s)	{ exit_error("Cannot open file: "+s);	}
inline void exit_cannotFind(const std::string& s)	{ exit_error("Cannot find file: "+s);	}
inline void exit_if        (bool condition)			{ exit_error("condition");	}

//---------------- string printing ----------------

// print_text automatically break lines between words and format each paragraph with regard to hanging_spaces, prefix_string, and line_width.
// print_text use a "priority" index to decide whether to print a paragraph. The priority of a paragraph is identified by __PRTXTFMT_PRIORITY__# at the begining of the paragraph.
// priority number is 0+, where 0 is the top priority. Only the text with priority <= priority_threshold will be printed. Therefore, by default (threshold=INT_MAX) all will be printed.
std::string print_text(std::string& text, const int hanging=0, const std::string prefix="", const int max_length=DefTermWidth, const int priority_threshold=INT_MAX);
std::string print_a_line(const std::string& pattern, int line_length=DefTermWidth);
std::string print_a_line(const std::pair<std::string,int>& pattern, int line_length=DefTermWidth);

template <typename T> // T = std::string / std::pair<std::string,int>
inline std::string print_left_aligned(const std::string& text, const T& pattern, int line_length=DefTermWidth) {
	std::string out_str;
	int n = line_length - text.length();
	out_str += text;
	out_str += print_a_line(pattern, n);
	return out_str;
}

template <typename T> // T = std::string / std::pair<std::string,int>
inline std::string print_right_aligned(const std::string& text, const T& pattern, int line_length=DefTermWidth) {
	std::string out_str;
	int n = line_length - text.length();
	out_str += print_a_line(pattern, n);
	out_str += text;
	return out_str;
}

template <typename T> // T = std::string / std::pair<std::string,int>
inline std::string print_in_middle(const std::string& text, const T& pattern, int line_length=DefTermWidth) {
	std::string out_str;
	int n2 = (line_length - text.length()) / 2;
	int n1 = line_length - text.length() - n2;
	out_str += print_a_line(pattern, n1);
	out_str += text;
	out_str += print_a_line(pattern, n2);
	return out_str;
}

namespace line_pattern {
	// controls and basic latin characters, can be used universally
	const std::pair<std::string,int> round_big      = std::pair<std::string,int>( "()"   , 2 );
	const std::pair<std::string,int> round_small    = std::pair<std::string,int>( "o"    , 1 );
	const std::pair<std::string,int> squared        = std::pair<std::string,int>( "[]"   , 2 );
	const std::pair<std::string,int> diamond        = std::pair<std::string,int>( "<>"   , 2 );
	const std::pair<std::string,int> zigzag_big     = std::pair<std::string,int>( "/\\"  , 2 );
	const std::pair<std::string,int> zigzag_zip     = std::pair<std::string,int>( "^v"   , 2 );
	const std::pair<std::string,int> zigzag_small   = std::pair<std::string,int>( "v"    , 1 );
	const std::pair<std::string,int> wave           = std::pair<std::string,int>( "~"    , 1 );
	const std::pair<std::string,int> dash           = std::pair<std::string,int>( "_ "   , 2 );
	const std::pair<std::string,int> long_dash      = std::pair<std::string,int>( "__  " , 4 );
	const std::pair<std::string,int> dash_dot       = std::pair<std::string,int>( "_."   , 2 );
	const std::pair<std::string,int> dash_dash_dot  = std::pair<std::string,int>( "__."  , 3 );
	const std::pair<std::string,int> dash_dot_sparse= std::pair<std::string,int>( "_ . " , 4 );
	const std::pair<std::string,int> dot            = std::pair<std::string,int>( "."    , 1 );
	const std::pair<std::string,int> dot_sparse     = std::pair<std::string,int>( ". "   , 2 );
	const std::pair<std::string,int> dot_double     = std::pair<std::string,int>( ":"    , 1 );
	const std::pair<std::string,int> single_line    = std::pair<std::string,int>( "-"    , 1 );
	const std::pair<std::string,int> double_line    = std::pair<std::string,int>( "="    , 1 );
	const std::pair<std::string,int> thick_line     = std::pair<std::string,int>( "#"    , 1 );
	const std::pair<std::string,int> asterisk       = std::pair<std::string,int>( "*"    , 1 );
	const std::pair<std::string,int> chain_like     = std::pair<std::string,int>( "-="   , 2 );
	// other characters, may have problem in some machines
	const std::pair<std::string,int> the_great_wall = std::pair<std::string,int>( "𐅭"    , 1 );
	const std::pair<std::string,int> black_4p_star  = std::pair<std::string,int>( "✦"    , 1 );
	const std::pair<std::string,int> white_4p_star  = std::pair<std::string,int>( "✧"    , 1 );
	const std::pair<std::string,int> black_circle   = std::pair<std::string,int>( "●"    , 1 );
	const std::pair<std::string,int> black_diamond  = std::pair<std::string,int>( "❖"    , 1 ); // previously black_diamond_minus_x
	const std::pair<std::string,int> tripple_line   = std::pair<std::string,int>( "≡"    , 1 );
	const std::pair<std::string,int> tripple_line_w = std::pair<std::string,int>( "Ξ"    , 1 );
	const std::pair<std::string,int> big_wave       = std::pair<std::string,int>( "∿"    , 1 );
	const std::pair<std::string,int> double_wave    = std::pair<std::string,int>( "≈"    , 1 );
	const std::pair<std::string,int> tripple_wave   = std::pair<std::string,int>( "≋"    , 1 );
	const std::pair<std::string,int> long_wave      = std::pair<std::string,int>( "◜◝◟◞" , 4 );
	const std::pair<std::string,int> chain          = std::pair<std::string,int>( "⋐⋑"   , 2 );
	const std::pair<std::string,int> open_ctr_cross = std::pair<std::string,int>( "✜"    , 1 );
	const std::pair<std::string,int> one_eighth     = std::pair<std::string,int>( "▬"    , 1 );
	const std::pair<std::string,int> one_quarter    = std::pair<std::string,int>( "▂"    , 1 );
	const std::pair<std::string,int> three_eighths  = std::pair<std::string,int>( "▃"    , 1 );
	const std::pair<std::string,int> half           = std::pair<std::string,int>( "▄"    , 1 );
	const std::pair<std::string,int> five_eighths   = std::pair<std::string,int>( "▅"    , 1 );
	const std::pair<std::string,int> three_quarters = std::pair<std::string,int>( "▆"    , 1 );
	const std::pair<std::string,int> seven_eighths  = std::pair<std::string,int>( "▇"    , 1 );
	const std::pair<std::string,int> full_block     = std::pair<std::string,int>( "█"    , 1 );
	const std::pair<std::string,int> light_shade    = std::pair<std::string,int>( "░"    , 1 );
	const std::pair<std::string,int> medium_shade   = std::pair<std::string,int>( "▒"    , 1 );
	const std::pair<std::string,int> quadrant       = std::pair<std::string,int>( "▚"    , 1 );
}

//-----file open/close--------

inline bool openfile_successfully(std::fstream& file, std::string filename, const std::ios_base::openmode& mode)
{
	filename_change_home_path(filename);
	if (str_has(filename,"/")) mkdir(substr_before_rfind(filename,"/"));
	file.open(filename.c_str(),mode);
	if		(file.is_open())	return true;
	else	{file.clear();		return false;}
}

#define openfile_or_exit(f,n,m)		if ( !openfile_successfully( (f),(n),(m) ) ) exit_cannotOpen(n)
#define openfile_or_break(f,n,m)	if ( !openfile_successfully( (f),(n),(m) ) ) break
#define openfile_or_continue(f,n,m)	if ( !openfile_successfully( (f),(n),(m) ) ) continue

inline void closefile(std::fstream & file)	{ file.close();	file.clear(); }

// ----- boost::iostreams file open/close -----

bool openInpFile(boost::iostreams::filtering_istream& f, std::string n);

// I have tried std::ios::app, but the result is weird: the file can be gunzip/bunzip2 in command line and is fine;
// but when I use boost::iostreams to read it again, the appended content is not accessible for .gz, while number of fields
// is wrong for .bz2. In this trial I used boost 1.46, zlib 1.2.5, bzip2 1.0.6, gcc 4.2
bool openOutFile(boost::iostreams::filtering_ostream& f, std::string n);

// remember 2 things: 
// a) A macro is 2 centenses; if the macro is conditional, use {}, for example, if (..) { openInpFile_or_exit(); .. }
// b) If {} is used, the f is used only within the {}.
#define openInpFile_or_exit(f,n)		boost::iostreams::filtering_istream f; if ( !openInpFile( (f),(n) ) ) exit_cannotOpen(n)
#define openInpFile_or_break(f,n)		boost::iostreams::filtering_istream f; if ( !openInpFile( (f),(n) ) ) break
#define openInpFile_or_continue(f,n)	boost::iostreams::filtering_istream f; if ( !openInpFile( (f),(n) ) ) continue
#define openOutFile_or_exit(f,n)		boost::iostreams::filtering_ostream f; if ( !openOutFile( (f),(n) ) ) exit_cannotOpen(n)
#define openOutFile_or_break(f,n)		boost::iostreams::filtering_ostream f; if ( !openOutFile( (f),(n) ) ) break
#define openOutFile_or_continue(f,n)	boost::iostreams::filtering_ostream f; if ( !openOutFile( (f),(n) ) ) continue

inline void closefile(boost::iostreams::filtering_ostream& f) {	f.flush(); f.reset(); }
inline void closefile(boost::iostreams::filtering_istream& f) {	f.reset(); } 

//----other read/write utilities----

// return whether the string is found. 
// If found, file loc=end of string; otherwise, file.eof()=true;
int read_until(std::istream& file,const std::string& s);

// to compile this function, -lboost_thread -lrt on Linux and -lboost_thread on Mac. Recommend to use the next function.
void write_file_try_lock(const std::string& text, const std::string& filename, const std::ios_base::openmode mode=std::fstream::out, int num_tries=86400, int sleep_sec=1, int incr=0);
void write_file_lock(const std::string& text, const std::string& filename, const std::ios_base::openmode mode=std::fstream::out);

//------------ data structure --------------
//previously T is int.
//read integer formatted like "1,5-8,11,15-"
//if the 2nd num is empty, then use std::numeric_limits<T>::max()
//lower bound must be specified, because -# is negative number
//although "#-" can be in the middle, it's highly recommended to not do it
//return whether something's read
template <typename T>
inline int parse_num_pairs(std::istream& ss, std::vector< std::pair<T,T> >& val, bool no_rev=false) // parse_int_pair
{
	size_t old_val_size=val.size();
	T i1,i2;
	char c;
	for (c=ss.peek(); (c>='0' && c<='9') || c=='-'; c=ss.peek())
	{
		if (!( ss>>i1 )) break;
		c=ss.peek();
		if (c==',') { ss.get(); val.push_back(std::pair<T,T>(i1,i1)); continue; }
		if (c=='-') 
		{
			ss.get();
			c=ss.peek();
			if ((c>='0' && c<='9') || c=='-')
			{
				if (!( ss>>i2 )) { val.push_back(std::pair<T,T>(i1,std::numeric_limits<T>::max())); break; } //prv i1,i1
				if (i2>=i1)			val.push_back(std::pair<T,T>(i1,i2));
				else if (!no_rev)	val.push_back(std::pair<T,T>(i1,i2));
					 else			exit_error("Input integer-pair reversed.");
				c=ss.peek();
				if (c==',') { ss.get(); continue; }
				else		break;
			}
			else
			{
				val.push_back(std::pair<T,T>(i1,std::numeric_limits<T>::max()));
				if (c==',') { ss.get(); continue; }
				else		break;
			}
		}
		else
		{
			val.push_back(std::pair<T,T>(i1,i1));
			break;
		}
	}
	return val.size()!=old_val_size;
}

//input string become the remaining string, return 
// consider to rename to extract_ for consistency
template <typename T>
inline int parse_num_pairs(std::string& str, std::vector< std::pair<T,T> >& val, bool no_rev=false) // parse_int_pair
{
	std::istringstream ss(str);
	int n=parse_num_pairs(ss,val,no_rev);
	str.clear();
	getline(ss,str,(char)EOF);
	return n;
}

template <typename T>
inline void print_num_pairs(const std::vector< std::pair<T,T> >& val, std::ostream& file) // print_int_pair
{
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin(); i!=val.end(); ++i) 
	{
		if (i!=val.begin()) file<<',';
		file<<i->first;
		if (i->first!=i->second) file<<'-'; else continue;
		if (i->second!=std::numeric_limits<T>::max())  file<<i->second;
	}
}

template <typename T>
inline void convert_num_pairs(std::vector< std::pair<T,T> >& val, T max_num) // convert_int_pair
{
	typedef typename std::vector< std::pair<T,T> >::iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->second==std::numeric_limits<T>::max()) i->second=max_num;
	}
}

template <typename T>
inline bool num_pairs_includes(const std::vector< std::pair<T,T> >& val, T num) // within_int_pair
{
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->first <= i->second)	{ if (num>=i->first && num<=i->second) return true; }
		else						{ if (num<=i->first && num>=i->second) return true; }
	}
	return false;
}

template <typename T>
inline bool num_pairs_overlaps(const std::vector< std::pair<T,T> >& val, T num1, T num2)
{
	if (num1>num2) { T c(num1); num1=num2; num2=c; }
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->first <= i->second)	{ if (num2>=i->first  && num1<=i->second) return true; }
		else						{ if (num2>=i->second && num1<=i->first)  return true; }
	}
	return false;
}

template <typename T>
inline void num_pairs_to_vector(std::vector< std::pair<T,T> >& val, std::deque<bool>& vec)
{
	typedef typename std::vector< std::pair<T,T> >::iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		int l = std::min(i->first,i->second);
		int r = std::max(i->first,i->second);
		for (int j=l; j<=r; ++j) vec_deq_set(vec,j,true);
	}
}

template <typename T>
inline bool reg_num_pairs_includes(const std::vector< std::pair<T,T> >& val, T num) // within_int_pair
{
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
		if (num>=i->first && num<=i->second) return true;
	return false;
}

template <typename T>
inline bool reg_num_pairs_overlaps(const std::vector< std::pair<T,T> >& val, T num1, T num2)
{
	if (num1>num2) { T c(num1); num1=num2; num2=c; }
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
		if (num2>=i->first  && num1<=i->second) return true;
	return false;
}

template <typename T>
inline bool num_pairs_includes_inf(const std::vector< std::pair<T,T> >& val) // inf_in_int_pair
{
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->second==std::numeric_limits<T>::max()) return true;
	}
	return false;
}

template <typename T>
inline T num_pairs_max(const std::vector< std::pair<T,T> >& val) // max_of_int_pair
{
	T m=std::numeric_limits<T>::min();
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->first > m) m=i->first;
		if (i->second > m) m=i->second;
	}
	return m;
}

template <typename T>
inline T num_pairs_min(const std::vector< std::pair<T,T> >& val) // min_of_int_pair
{
	T m=std::numeric_limits<T>::max();
	typedef typename std::vector< std::pair<T,T> >::const_iterator vIt;
	for (vIt i=val.begin();i!=val.end();++i)
	{
		if (i->first < m) m=i->first;
		if (i->second < m) m=i->second;
	}
	return m;
}

// --------- date and time ---------

std::string macro_date_to_boost(const std::string& in); // macro_date_to_boost(string(__DATE__)); "Nov 20 2011" => 2011-11-20. same as bash command:  date +"%Y-%m-%d"
bool date_incorrect();		// whether system date is correct by checking the internet
inline boost::gregorian::date	today_s_date()		{ return boost::gregorian::day_clock::local_day(); }
inline std::string				today_s_date_str()	{ return boost::gregorian::to_simple_string(today_s_date()); }

// usage: auto ds=duration_start(); cerr<<"Duration="<<duration_end(ds)<<std::endl;
boost::posix_time::ptime			duration_start();
boost::posix_time::time_duration	duration_end(const boost::posix_time::ptime& time_begin);

// =================== file processing, OS safe ===================

// for reading, EOF is defined (as is usual) as -1, or the unsigned decimal value 255 ('\377' or '\xff')
// for writing, EOF in DOS/Win is Cntrl-Z, 0x1A or 26. For most Unix style machines, it's Cntrl-D, 0x04, or 4.
// you can write "int my_it=0; for (each_line(file),++my_it)". But nothing before each_xxx.
#define each_row(f,v)		 int row_num=0; next_row(f) && get_fields(f,v); ++row_num, skip_one_line(f)
#define each_line(f)		 int row_num=0; next_row(f);					++row_num, skip_one_line(f)
#define each_line_sync2(f,g) int row_num=0; next_row(f) && next_row(g);		++row_num, skip_one_line(f),skip_one_line(g)
#define all_lines(f)		 int row_num=0; f.peek()!=EOF;					++row_num, skip_one_line(f)
#define all_lines_sync2(f,g) int row_num=0; f.peek()!=EOF && g.peek()!=EOF;	++row_num, skip_one_line(f),skip_one_line(g)

bool is_EndOfLine(char c);						// returns whether c is an EOL/EOF
bool is_blank_row(std::istream & file);			// returns whether file is at EOL/EOF
bool next_row(std::istream & file);				// returns whether has a next row (skip comments, skip whitespaces, stop at blank row)
bool skip_EndOfLine(std::istream & file);		// returns whether EOL extracted, skip the EOL if right at it
bool skip_one_line(std::istream & file);		// returns whether EOL extracted, skip to  EOL and extract it
bool skip_spaces(std::istream & file, int n);	// returns whether n SPACEs are extracted
bool skip_whitespaces(std::istream & file);		// returns true, skip ANY SPACEs / TABs
bool skip_blank_lines(std::istream & file);		// returns true, skip ANY lines without contents or whitespaces
bool skip_comments(std::istream & file);		// returns true, skip ANY lines begin with a # (use at line starts)

// return whether field exist (field 0 always exist, but could be blank). field# starts from 0.
bool get_field(std::istream & file, int field, std::string& os, const std::string& dlt="\t ", bool skip_eol=false);

// return whether read something (equal to !is_blank_row(file)). field# starts from 0.
bool get_fields(std::istream & file, std::vector<std::string>& ov, const std::string& dlt="\t ", bool skip_eol=false);

// =================== OS Independent file reading functions ====================

// robust to DOS/LINUX/MAC new line characters, and is fast
std::istream& safeGetline    (std::istream& is, std::string& t); // will  clear t
std::istream& safeGetline_add(std::istream& is, std::string& t); // won't clear t

// extract one char, return char by value
int osi_getchar(std::istream & file);

// extracts EOL, modify string
void osi_getline(std::istream & file,std::string & os);

// EndOfLine is not extracted, modify string
void osi_get(std::istream & file,std::string & os);

// extracts EOL, return string by value
std::string osi_getline(std::istream & file);

// EndOfLine is not extracted, return string by value
std::string osi_get(std::istream & file);

// extracts EOL, modify string, using es as EOL character
void osi_getline(std::istream & file,std::string& os,char es);

// EndOfLine is not extracted, modify string, using es as EOL character
void osi_get(std::istream & file,std::string& os,char es);

// copy one line from one file to another, EndOfLine is not extracted
void osi_copy_row(std::istream & ifle, std::ostream & ofle);

// copy one line from one file to another, extracts EOL
void osi_copy_line(std::istream & ifle, std::ostream & ofle);

// extracts EOL, modify string, using endc as EOL character
void osi_getline_allowing_brackets_or_quotes(std::istream & file,std::string& outs,char endc);

// read one line of file and tokenize it in UNIX-style, EOL is extracted, robust to diff EOL / no EOL
// 1) "\c" => \c, except when c = \`$" which results in c only.
// 2)  \c  =>  c, even when c = t n SPACE
// 3) '\c' => \c
// 4) if last char is \ (no preseeding \) then go to next line (Unix do this also to `")
void tokenize(std::istream& file, std::vector<std::string>& result);

// -------------------- serie --------------------

// a serie's format: @#[/option1/option2/..]//  (# is int & >0)
int serie_fill_string(std::string& str,int seq, int val, std::map<int,std::vector<std::string> >& rplmap,char starting_char='@');
int serie_fill_number(std::string& str,int seq, double val,int width=0, char fill='0',char starting_char='@');
int number_of_series(const std::string& str,char starting_char='@');
int max_serie_ID(const std::string& str,char starting_char='@');

// ------------- data input -----------------

class input_exception { public: const char *ShowReason() const { return "unknown answer"; } };

std::string str_YesOrNo(const bool YesOrNo);
void		clear_cin();
std::string input_line(const std::string& s="");
std::string input_paragraph(const std::string& s=""); // Unix: use ^d to input EOF. The last line must has EOL, otherwise it does not go to cin.
bool		input_yes(const std::string& s="");
bool		input_yes(const bool dflt, const std::string& s); // s cannot be default to "" due to the previous func, or input_yes("hello") => input_yes(true)
char		input_char(const std::string& s="");
char		input_char(const char dflt, const std::string& s);// s cannot be default to "" due to the previous func, or input_yes("hello") => input_yes(true)

template <typename T>
inline void input_number(T& dest, const std::string& s="")
{
	std::cout << s << ": ";
	for (;;)
	{
		try { dest=boost::lexical_cast<T>(input_line()); return; }
		catch (boost::bad_lexical_cast &) 
		{
			std::cout<<"Not a required number type, please try again: ";
			continue;
		}
	}
}

template <typename T>
inline void input_number(T& dest, const T dflt, const std::string& s)
{
	std::cout << s << " [" << dflt << "]: ";
	for (;;)
	{
		try
		{
			std::string instr=input_line();
			if (instr.empty()) dest=dflt;
			else dest=boost::lexical_cast<T>(instr);
			return;
		}
		catch (boost::bad_lexical_cast &) 
		{
			std::cout<<"Not a required number type, please try again: ";
			continue;
		}
	}
}

// =================== variable translation ====================

template <typename T1, typename T2>
class TranslateVariable
{
public:
	typedef typename std::map<T1,T2>::const_iterator const_iterator;
	std::map<T1,T2> reference_data;

	TranslateVariable(){}
	TranslateVariable(const std::string& filename) {
		setup(filename);
	}
	TranslateVariable(const std::string& filename, std::istream& alt) {
		boost::iostreams::filtering_istream file;
		if ( openInpFile(file,filename) )	setup(file);
		else								setup(alt);
	}
	void setup	(const std::string& filename) {
		reference_data.clear();
		openInpFile_or_exit(file,filename);
		setup(file);
		closefile(file);
	}
	void setup	(std::istream& file) {
		reference_data.clear();
		for (each_line(file))
		{
			T1 x;
			T2 y;
			file >> x >> y;
			reference_data[x]=y;
		}
	}
	bool exist (const T1& var) {
		const_iterator found = reference_data.find(var);
		return found!=reference_data.end();
	}
	T2		tr	(const T1& var) {
		const_iterator found = reference_data.find(var);
		if (found==reference_data.end())
		{
			std::stringstream ss;
			ss << "cannot translate " << var;
			exit_error(ss.str());
		}
		return found->second;
	}
};

template <typename T>
class TranslateVariableSameType : public TranslateVariable<T,T>
{
public:
	T tr_if_exist	(const T& var) {
		typename std::map<T,T>::const_iterator found = this->reference_data.find(var);
		if (found!=this->reference_data.end()) return found->second;
		return var;
	}
};

class progress_timer {
private:
	struct progress_timer_data;
	progress_timer_data * d;
public:
	progress_timer();
	~progress_timer();
	progress_timer(const progress_timer& othr);
	progress_timer& operator=(const progress_timer& orig);
	
	// usage
	void prefix(const std::string& s);
	void start(unsigned long N);
	unsigned long operator+=( unsigned long increment );
	unsigned long  operator++() { return operator+=( 1 ); }
	void finish();
};

#endif
