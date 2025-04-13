#pragma once

#include <iostream>

namespace basis {

template<typename Value>
class blocked_range {
public :
	typedef Value const_iterator;
	typedef std::size_t size_type;
	
	blocked_range(const_iterator b, const_iterator e) : my_begin(b), my_end(e) {}

	const_iterator begin() const {return my_begin;}
	const_iterator end() const {return my_end;}

	const size_type size() const { return size_type(my_end-my_begin); }
	const bool empty() const {return !(my_begin < my_end);}

private : 
	const_iterator my_begin, my_end;
};

};			// basis