////////////////////////////////
/// usage : 1.	file i/o utilities.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_FILE_H
#define CN_HUST_GOAL_COMMON_FILE_H


#include <iostream>
#include <vector>

#include "GOAL/Typedef.h"


namespace goal {
namespace file {

struct Lines {
	#pragma region wrapper for std::vector<char*>
	using Container = Vec<char*>;
	using SizeType = typename Container::size_type;
	using ValType = typename Container::value_type;
	using IterType = typename Container::iterator;
	using ConstIterType = typename Container::const_iterator;

	SizeType size() const { return lines.size(); }
	bool empty() const { return lines.empty(); }

	void resize(SizeType newSize) { lines.reserve(newSize); }
	void reserve(SizeType newCapacity) { lines.reserve(newCapacity); }

	ValType& operator[](SizeType i) { return lines[i]; }
	const ValType& operator[](SizeType i) const { return lines[i]; }

	IterType begin() { return lines.begin(); }
	ConstIterType begin() const { return lines.begin(); }
	IterType end() { return lines.end(); }
	ConstIterType end() const { return lines.end(); }

	void push_back(const ValType& v) { lines.push_back(v); }
	void push_back(ValType&& v) { lines.push_back(v); }
	void pop_back() { lines.pop_back(); }
	#pragma endregion wrapper for std::vector<char*>

	void copyLines(const Lines& a) {
		lines.resize(a.lines.size());
		for (SizeType i = 0; i < lines.size(); ++i) {
			lines[i] = &charBuf[0] + (a.lines[i] - &a.charBuf[0]);
		}
	}
	void safeMove(Lines& a) {
		const char* oldCharBuf = &a.charBuf[0];
		charBuf = std::move(a.charBuf);
		if (&charBuf[0] == &a.charBuf[0]) {
			lines = std::move(a.lines);
		} else {
			lines.resize(a.lines.size());
			for (SizeType i = 0; i < lines.size(); ++i) {
				lines[i] = &charBuf[0] + (a.lines[i] - oldCharBuf);
			}
		}
	}

	Lines() {}

	Lines(const Lines& a) : charBuf(a.charBuf) { copyLines(a); }
	//Lines(Lines&& a) = default;
	Lines(Lines&& a) { safeMove(a); }

	Lines& operator=(const Lines& a) {
		if (this != &a) {
			charBuf = a.charBuf;
			copyLines(a);
		}
		return *this;
	}
	//Lines& operator=(Lines&& a) = default;
	Lines& operator=(Lines&& a) {
		if (this != &a) { safeMove(a); }
		return *this;
	}

	Str charBuf;
	Vec<char*> lines;
};

Str readAllText(const Str& path);
Lines readAllLines(const Str& path);
Lines readAllLines(const Str& path, char commentChar);
Lines readAllLinesSkipEmpty(const Str& path);
Lines readAllLinesSkipEmpty(const Str& path, char commentChar);


template<typename T>
T next(std::istream& is) {
	T obj;
	is >> obj;
	return obj;
}

inline void rewind(std::ostream& os, int charNum = 1) {
	os.seekp(-charNum, std::ios_base::cur);
}

}
}


#endif // CN_HUST_GOAL_COMMON_FILE_H
