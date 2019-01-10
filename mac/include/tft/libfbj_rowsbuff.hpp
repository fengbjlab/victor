// body file: libfbj_rowsbuff.cpp

#ifndef LIBFBJ_ROWS_BUFFER
#define LIBFBJ_ROWS_BUFFER

#include <string>
#include <vector>
#include <deque>

// By default, only the last row is stored, unless set_oldest()/keep_all_input() is called.
// In accessing data:
//   1) row=0 is current row, -1 if last, 1 is next ...
//   2) if there's no data, complains and exit
//   3) if there're data but query an inexistent row or column, throw an exception
// To iterate all rows: use rowID in [0,size).
// Therefore, it's possible to get the next row: keep_all_input, allow_future_rows, then iterate rowID.

class RowsBufferException: public std::exception
{
	virtual const char* what() const throw()
	{
		return "Query an inexistent row or column.";
	}
};

typedef long long RowNum_t;

class RowsBuffer {

public:
	
	typedef std::string string;
	typedef std::vector<string> RowData;
	typedef std::deque<RowData>::iterator				iterator;
	typedef std::deque<RowData>::reverse_iterator		reverse_iterator;
	typedef std::deque<RowData>::const_iterator			const_iterator;
	typedef std::deque<RowData>::const_reverse_iterator	const_reverse_iterator;
	
	// setup
	~RowsBuffer();
	RowsBuffer();
	RowsBuffer(const RowsBuffer& othr);
	RowsBuffer& operator=(const RowsBuffer& orig);
	void		set_oldest(RowNum_t row);				// set the oldest row for query, must be <=0
	void		keep_all_input();						// do not remove the oldest row when adding
	void		allow_future_rows();					// allow get the future rows
	
	// change data
	RowData&	push_back();							// add a row w/o  data
	void		push_back(RowData&);					// add a row with data
	void		pop_back();								// remove the last row, will reduce rows_added
	void		clear();								// remove all rows
	void		clear_ExceptTheLastRow();				// remove all rows but keep only the last one
	
	// get data: row=(-inf,inf), col=[0,inf). row=0 is current, -1 is last, 1 is next.
	RowNum_t&	rowID();								// get/set current row. 0-based
	RowNum_t	rows_added() const;						// total number of rows added -1. Used as a 0-based row number.
	string&		operator() (RowNum_t row, size_t col);	// return data[target_row][col] by reference
	string&		operator[] ( size_t col );				// return data[current_row][col] by reference
	RowData&	target_row(RowNum_t row);				// get target row. row is relative row# -- 0 is current, >0 is future, <0 is previous
	RowData&	current_row();							// get current row. equivalent of target_row(0).
	
	// access data
	iterator				begin()		const;
	iterator				end()		const;
	const_iterator			cbegin()	const;
	const_iterator			cend()		const;
	reverse_iterator		rbegin()	const;
	reverse_iterator		rend()		const;
	const_reverse_iterator	crbegin()	const;
	const_reverse_iterator	crend()		const;
	bool					empty()		const;
	bool					has_room()	const;
	size_t					size()		const;

private:
	struct RowsBufferData;
	RowsBufferData * d;
};

#endif

/* to debug
int main (int argc, char * const argv[]) 
{
	RowsBuffer hr;
	hr.set_oldest(-2);
	vector<string> v1(1,"11"); v1.push_back("12"); v1.push_back("13");
	vector<string> v2(1,"21"); v2.push_back("22"); v1.push_back("23");
	vector<string> v3(1,"31"); v3.push_back("32"); v1.push_back("33");
	vector<string> v4(1,"41"); v4.push_back("42"); v1.push_back("43");
	hr.push_back(v1);
	hr.push_back(v2);
	hr.push_back(v3);
	hr.push_back(v4);
	cout<<hr(-1,2)<<endl;
	return 1;
}//*/
