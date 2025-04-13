#pragma once

#include "blocked_range.h"

namespace basis {

template<typename PageValue, typename RowValue=PageValue, typename ColValue=RowValue>
class blocked_range3d {
public :
	typedef blocked_range<PageValue> page_range_type;
	typedef blocked_range<RowValue> row_range_type;
	typedef blocked_range<ColValue> col_range_type;

	blocked_range3d( PageValue page_begin, PageValue page_end, RowValue  row_begin,  RowValue row_end, ColValue  col_begin,  ColValue col_end ) 
	: my_pages(page_begin,page_end), my_rows(row_begin,row_end), my_cols(col_begin,col_end) {}

	const page_range_type& pages() const {return my_pages;}
	const row_range_type& rows() const {return my_rows;}
	const col_range_type& cols() const {return my_cols;}

	const bool empty() const { return my_pages.empty() || my_rows.empty() || my_cols.empty();	}

private : 
	page_range_type my_pages;
	row_range_type my_rows;
	col_range_type my_cols;
};

};			// basis