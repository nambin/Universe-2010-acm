#pragma once

#include "blocked_range.h"

namespace basis {

template<typename RowValue, typename ColValue=RowValue>
class blocked_range2d {
public :
	typedef blocked_range<RowValue> row_range_type;
	typedef blocked_range<ColValue> col_range_type;

	blocked_range2d( RowValue row_begin, RowValue row_end, ColValue col_begin, ColValue col_end ) : my_rows(row_begin,row_end),	my_cols(col_begin,col_end) {}

	const row_range_type& rows() const {return my_rows;}
	const col_range_type& cols() const {return my_cols;}

	const bool empty() const { return my_rows.empty() || my_cols.empty();	}

private : 
	row_range_type my_rows;
	col_range_type my_cols;
};

};			// basis