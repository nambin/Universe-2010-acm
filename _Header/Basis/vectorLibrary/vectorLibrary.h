#ifndef _BASIS_VECTOR_LIBRARY_H_
#define _BASIS_VECTOR_LIBRARY_H_

#include <vector>
#include "../FileIO/FileIO.h"

namespace basis {

class vectorLibrary {
public : 
	template<class T> static bool importBinary(std::ifstream&, std::vector<T>&);
	template<class T> static bool exportBinary(std::ofstream&, const std::vector<T>&);

	template<class T> static void zero(std::vector<T> &v);
	template<class T> static void set(std::vector<T> &v, T val);
	template<class T> static bool isContain(const std::vector<T> &a, T e);
	template<class T> static bool add_unique(std::vector<T> &a, T e);
	template<class T> static void insert(std::vector<T> &a, size_t index, T e);
	template<class T> static void erase(std::vector<T> &a, size_t index);
	template<class T> static void erase_swap(std::vector<T> &a, size_t index);
	template<class T> static void erase_unordered(std::vector<T> &a, size_t index);
	template<class T> static void erase_unordered_swap(std::vector<T> &a, size_t index);
	template<class T> static void find_and_erase_unordered(std::vector<T> &a, const T &doomed_element);
	template<class T> static void replace_once(std::vector<T> &a, const T &old_element, const T &new_element);

	template<class T> static void append(std::vector<T>&, const std::vector<T>&);
};

template<class T> 
void vectorLibrary::append(std::vector<T>& dst, const std::vector<T>& src) {
	for(std::size_t i=0;i<src.size();++i)
		dst.push_back(src[i]);
}

template<class T> 
bool vectorLibrary::importBinary(std::ifstream& ifs, std::vector<T>& vec) {
	std::size_t _size;
	::importBinary(ifs, _size);

	vec.clear(); vec.resize(_size);
	for(std::size_t i=0;i<vec.size();++i)
		::importBinary(ifs, vec[i]);
	return true;
}

template<class T> 
bool vectorLibrary::exportBinary(std::ofstream& ofs, const std::vector<T>& vec) {
	::exportBinary(ofs, vec.size());
	for(std::size_t i=0;i<vec.size();++i)
		::exportBinary(ofs, vec[i]);
	return true;
}

template<class T>
void vectorLibrary::zero(std::vector<T> &v) { 
	set(v, T(0.0));
}

template<class T> 
void vectorLibrary::set(std::vector<T> &v, T val) {
	for(size_t i=0;i<v.size();++i) 
		v[i] = val; 
}

template<class T>
bool vectorLibrary::isContain(const std::vector<T> &a, T e) {
	for(size_t i=0;i<a.size();++i) {
		if(a[i] == e) {
			return true;
		}
	}
	return false;
}

template<class T>
bool vectorLibrary::add_unique(std::vector<T> &a, T e) {
	for(size_t i=0;i<a.size();++i) {
		if(a[i] == e) {
			return false;
		}
	}
	a.push_back(e);
	return true;
}

template<class T>
void vectorLibrary::insert(std::vector<T> &a, size_t index, T e) {
	a.push_back(a.back());
	for(size_t i=a.size()-1; i>index; --i)
		a[i] = a[i-1];
	a[index] = e;
}

template<class T>
void vectorLibrary::erase(std::vector<T> &a, size_t index) {
	for(size_t i=index; i<a.size()-1; ++i)
		a[i] = a[i+1];
	a.pop_back();
}

template<class T>
void vectorLibrary::erase_swap(std::vector<T> &a, size_t index) {
	for(size_t i=index; i<a.size()-1; ++i)
		std::swap(a[i], a[i+1]);
	a.pop_back();
}

template<class T>
void vectorLibrary::erase_unordered(std::vector<T> &a, size_t index) {
	a[index] = a.back();
	a.pop_back();
}

template<class T>
void vectorLibrary::erase_unordered_swap(std::vector<T> &a, size_t index) {
	swap(a[index], a.back());
	a.pop_back();
}

template<class T>
void vectorLibrary::find_and_erase_unordered(std::vector<T> &a, const T &doomed_element) {
	for(size_t i=0; i<a.size(); ++i) {
		if(a[i] == doomed_element) {
			erase_unordered(a, i);
			return;
		}
	}
}

template<class T>
void vectorLibrary::replace_once(std::vector<T> &a, const T &old_element, const T &new_element) {
	for(size_t i=0; i<a.size(); ++i) {
		if(a[i] == old_element) {
			a[i] = new_element;
			return;
		}
	}
}

};		// basis

#endif