////////////////////////////////
/// usage : 1.	string operation extensions.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_COMMON_STR_H
#define CN_HUST_GOAL_COMMON_STR_H


#include <string>

#include "GOAL/Typedef.h"


namespace goal {

template<typename T>
inline Str toStr(const T& obj) { return std::to_string(obj); }

inline bool startWith(const Str& str, const Str& substr) { return str.compare(0, substr.size(), substr) == 0; }


struct VarStr {
    static constexpr char Delimiter = '_';


    static Str concat(const char *obj) { return Str(obj); }

    template<typename T>
    static Str concat(const T &obj) { return toStr(obj); }

    template<typename T, typename ... Ts>
    static Str concat(const T &obj, Ts ... objs) { return concat(obj) + Delimiter + concat(objs ...); }

    static Str pad(int count, char fillChar = '0') { return Str(count, fillChar); }
};

}


#endif // CN_HUST_GOAL_COMMON_STR_H
