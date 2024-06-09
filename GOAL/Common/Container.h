////////////////////////////////
/// usage : 1.	extensions to STL containers.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_CONTAINER_H
#define CN_HUST_GOAL_COMMON_CONTAINER_H


#include <algorithm>


namespace goal {

template<typename Container>
typename Container::value_type& newBack(Container& container) {
	container.emplace_back();
	return container.back();
}

template<typename Container>
void unorderedRemove(Container& container, typename Container::iterator pos) {
	*pos = container.back();
	container.pop_back();
}

template<typename Container>
void reset(Container& container, typename Container::size_type size, const typename Container::value_type& value) {
	container.clear();
	container.resize(size, value);
}

template<typename Container>
bool contain(const Container& container, const typename Container::value_type& value) {
	return std::find(container.begin(), container.end(), value) != container.end();
}

template<typename Container>
bool containKey(const Container& container, const typename Container::key_type& value) {
	return container.find(value) != container.end();
}

template<typename Container, typename T>
typename Container::iterator find(Container& container, const T& value) {
	return std::find(container.begin(), container.end(), value);
}

template<typename Container, typename T>
typename Container::const_iterator find(const Container& container, const T& value) {
	return std::find(container.begin(), container.end(), value);
}

template<typename Container, typename Func>
typename Container::iterator findIf(Container& container, Func pred) {
	return std::find_if(container.begin(), container.end(), pred);
}

template<typename Container, typename Func>
typename Container::const_iterator findIf(const Container& container, Func pred) {
	return std::find_if(container.begin(), container.end(), pred);
}

template<typename Container, typename Func>
bool containIf(const Container& container, Func pred) {
	return findIf(container, pred) != container.end();
}

}


#endif // CN_HUST_GOAL_COMMON_CONTAINER_H
