#ifndef FACTORIAL_CLASS
#define FACTORIAL_CLASS

#include <map>
#include <cmath>
#include <sstream>

class FACTORIAL
{

private:

	typedef unsigned long	FACT_T;
	typedef long double		RESL_T;
	std::map<FACT_T,int> _o; // original, store n!
	bool calcualted;
	RESL_T value;
	
	void _inverse(const std::map<FACT_T,int>& _l, std::map<FACT_T,int>& _r) const {
		// intentionally there's no _r.clear() here, so that _r can be written aggregately
		std::map<FACT_T,int>::const_iterator it;
		for (it=_l.begin(); it!=_l.end(); ++it) _r[it->first] -= it->second;
	}
	
	void init_original	(RESL_T& v) { v=1; }
	void init_log10		(RESL_T& v) { v=0; }
	void init_ln		(RESL_T& v) { v=0; }
	void up_orginal	(RESL_T& v, const FACT_T& j) { v*=j; }
	void dn_orginal	(RESL_T& v, const FACT_T& j) { v/=j; }
	void up_log10	(RESL_T& v, const FACT_T& j) { v+=log10(j); }
	void dn_log10	(RESL_T& v, const FACT_T& j) { v-=log10(j); }
	void up_ln		(RESL_T& v, const FACT_T& j) { v+=log(j); }
	void dn_ln		(RESL_T& v, const FACT_T& j) { v-=log(j); }
	void (FACTORIAL::*init)(RESL_T&);
	void (FACTORIAL::*up)(RESL_T&, const FACT_T&);
	void (FACTORIAL::*dn)(RESL_T&, const FACT_T&);
	
	RESL_T _directly_calculate(std::map<FACT_T,int> _l) 
	{
		RESL_T v;
		(this->*init)(v);
		
		// first take out the max nominator and max denominator together, this not only largely reduce the number of * / operations
		// but also lower the chances that the intermediate value go beyond the extreme (too big or too close to 0)
		for (;;)
		{
			FACT_T max_pos=0;
			FACT_T max_neg=0;
			std::map<FACT_T,int>::reverse_iterator rit;
			for (rit=_l.rbegin(); rit!=_l.rend(); ++rit)
			{
				if		(rit->second>0 && max_pos==0) { max_pos=rit->first; --rit->second; }
				else if (rit->second<0 && max_neg==0) { max_neg=rit->first; ++rit->second; }
				if (max_pos && max_neg) break;
			}
			if (max_pos==0 && max_neg==0) break;
			if (max_pos > max_neg)	for (FACT_T j=max_pos; j>max_neg; --j) (this->*up)(v,j); // v*=j;
			else					for (FACT_T j=max_neg; j>max_pos; --j) (this->*dn)(v,j); // v/=j;
			if (max_pos==0 || max_neg==0) break;
		}
		
		// now do the remaining
		std::map<FACT_T,int>::iterator it;
		for (it=_l.begin(); it!=_l.end(); ++it) {
			if		(it->second>0)	{ for (FACT_T j=2;j<=it->first;++j) for (int i=0;i<it->second;++i) (this->*up)(v,j); } // v*=j; }
			else if (it->second<0)	{ for (FACT_T j=2;j<=it->first;++j) for (int i=0;i>it->second;--i) (this->*dn)(v,j); } // v/=j; }
		}
		
		return v;
	}
	
public:

	FACTORIAL():calcualted(true)	{ up=&FACTORIAL::up_orginal;	dn=&FACTORIAL::dn_orginal;	init=&FACTORIAL::init_original; (this->*init)(value); }
	void cal_original()				{ up=&FACTORIAL::up_orginal;	dn=&FACTORIAL::dn_orginal;	init=&FACTORIAL::init_original; (this->*init)(value); calcualted=false; }
	void cal_log10()				{ up=&FACTORIAL::up_log10;		dn=&FACTORIAL::dn_log10;	init=&FACTORIAL::init_log10;	(this->*init)(value); calcualted=false; }
	void cal_ln()					{ up=&FACTORIAL::up_ln;			dn=&FACTORIAL::dn_ln;		init=&FACTORIAL::init_ln;		(this->*init)(value); calcualted=false; }
	void clear()					{ _o.clear(); calcualted=true; (this->*init)(value); }
	void add(const FACT_T& n)		{ if (n<2) return; ++_o[n]; calcualted=false; }
	void sub(const FACT_T& n)		{ if (n<2) return; --_o[n]; calcualted=false; }
	void remove(const FACT_T& n)	{ if (n<2) return; --_o[n]; calcualted=false; }
	RESL_T result() 
	{
		if (calcualted) return value;
		value=_directly_calculate(_o);
		calcualted=true;
		return value;
	}
	bool operator<= (const FACTORIAL &other)  { 
		std::map<FACT_T,int> _c = this->_o; // combined
		_inverse(other._o , _c);
		if (init==&FACTORIAL::init_original)	return ( _directly_calculate(_c) <= 1 );
		else									return ( _directly_calculate(_c) <= 0 );
	}
	bool operator< (const FACTORIAL &other)  { 
		std::map<FACT_T,int> _c = this->_o; // combined
		_inverse(other._o , _c);
		if (init==&FACTORIAL::init_original)	return ( _directly_calculate(_c) < 1 );
		else									return ( _directly_calculate(_c) < 0 );
	}
	bool operator== (const FACTORIAL &other)  {
		std::map<FACT_T,int> _c = this->_o; // combined
		_inverse(other._o , _c);
		if (init==&FACTORIAL::init_original)	return ( _directly_calculate(_c) == 1 );
		else									return ( _directly_calculate(_c) == 0 );
	}
	std::string print()
	{
		std::stringstream ss;
		
		// nominator
		int id=0;
		for (auto &p:_o)
		{
			if		(p.second>1)	{ if (id++) ss<<"*"; else ss<<"("; ss<<p.first<<"^"<<p.second; }
			else if (p.second==1)	{ if (id++) ss<<"*"; else ss<<"("; ss<<p.first; }
		}
		if (id==0) ss<<"1"; else ss<<")";
		
		// denominator
		id=0;
		for (auto &p:_o)
		{
			if		(p.second<-1)	{ if (id++) ss<<"*"; else ss<<"/("; ss<<p.first<<"^"<<-p.second; }
			else if (p.second==-1)	{ if (id++) ss<<"*"; else ss<<"/("; ss<<p.first; }
		}
		if (id!=0) ss<<")";
		
		return ss.str();
	}
	
};

#endif
// ----------------- removed functions -----------------
// previously use expansion and cancellation to reduce the number of * and / operation
// but the expansion itself is too slow because of the c++ map operation ( in expand(), _r[j] etc)
// so I removed these functions, and use direct calculation from _o. It's much faster now.
// another good feature of direct calculation is that very small value become 0 rather than inf
// as shown in the third case of the example program.

/* private:
RESL_T _calculate_after_expansion(const map<FACT_T,int>& _l) const {
	RESL_T v=1;
	map<FACT_T,int>::const_iterator it;
	for (it=_l.begin(); it!=_l.end(); ++it) {
		if	(it->second>0)	{ for (int i=0;i<it->second;i++) v*=it->first; }
		else				{ for (int i=0;i>it->second;i--) v/=it->first; }
	}
	return v;
}
void _clean(map<FACT_T,int>& _l) {
begin_cleaning:
	map<FACT_T,int>::iterator it;
	for (it=_l.begin(); it!=_l.end(); ++it) 
		if (it->second==0) { _l.erase(it); goto begin_cleaning; }
}
void _expand(map<FACT_T,int> _l, map<FACT_T,int>& _r) { // _l store n! while _r store n(n-1)..2
	// intentionally there's no _r.clear() here, so that _r can be written aggregately
	map<FACT_T,int>::iterator it;
	for (;;)
	{
		FACT_T max_pos=0;
		FACT_T max_neg=0;
		for (it=_l.begin(); it!=_l.end(); ++it) 
		{
			if		(it->second>0)	{ if (it->first>max_pos) max_pos=it->first; }
			else if (it->second<0)	{ if (it->first>max_neg) max_neg=it->first; }
		}
		if (max_pos==0 || max_neg==0) break;
		_l[max_pos]--;
		_l[max_neg]++;
		if (max_pos > max_neg)
			for (FACT_T j=max_pos; j>max_neg; j--) _r[j]++;
		else 
			for (FACT_T j=max_neg; j>max_pos; j--) _r[j]--;
	}
	_clean(_l);
	for (it=_l.begin(); it!=_l.end(); ++it) 
		for (FACT_T j=2;j<=it->first;j++) _r[j]+=it->second;
	_clean(_r);
}

public:
RESL_T result_old() {
	if (calcualted) return value;
	_clean(_o);
	
	// expand to _e
	map<FACT_T,int> _e; // expanded
	_expand(_o,_e);
	
	// calculate result value
	value=_calculate_after_expansion(_e);
	calcualted=true;
	return value;
}
bool operator<= (const FACTORIAL &other)  { 
	map<FACT_T,int> _c = this->_o; // combined
	_inverse(other._o , _c);
	map<FACT_T,int> _e; // expanded
	_expand(_c, _e);
	return ( _calculate_after_expansion(_e) <= 1 );
}
bool operator< (const FACTORIAL &other)  { 
	map<FACT_T,int> _c = this->_o; // combined
	_inverse(other._o , _c);
	map<FACT_T,int> _e; // expanded
	_expand(_c, _e);
	return ( _calculate_after_expansion(_e) < 1 );
}
*/

/*/ factorial example
#include <iostream>
int main()
{	FACTORIAL f;
	f.add(15);
	f.remove(25);
	std::cout<<"15!/25! = "<<f.result()<<endl; // 15!/25! = 8.43051e-14
	f.add(48);
	f.remove(12);
	f.remove(3);
	std::cout<<"15!48!/25!12!3! = "<<f.result()<<endl; // 15!48!/25!12!3! = 3.64145e+38
	f.clear();
	f.add(7556);
	f.add(3480);
	f.remove(11036);
	std::cout << "7556!3480!/(7556+3480)! = "<<f.result()<<endl; // 7556!3480!/(7556+3480)! = inf >>>> should I make it to 0 ????
	return 0;
} //*/
