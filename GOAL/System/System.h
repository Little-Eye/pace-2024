////////////////////////////////
/// usage : 1.	operating system utilities.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef CN_HUST_GOAL_SYSTEM_SYSTEM_H
#define CN_HUST_GOAL_SYSTEM_SYSTEM_H


#include <iostream>
#include <cstdlib>
#include <cstdio>

#include "GOAL/Flag.h"
#include "GOAL/Typedef.h"


namespace goal {
namespace os {

static int exec(const Str& cmd) { return system(cmd.c_str()); }

static Str quote(const Str& s) { return ('\"' + s + '\"'); }

class StdIo {
    static int stdIn;
    static int stdOut;
    static int stdErr;

public:
    static void disableStdOut();
    static void restoreStdOut();
    #ifdef _OS_MS_WINDOWS
    static FILE* resetStdOut() { return freopen("CON", "w", stdout); }
    #else
    static FILE* resetStdOut() { return freopen("/dev/tty", "w", stdout); }
    #endif
};

struct MemorySize {
    using Unit = long long;

    static constexpr double Base = 1024;

    friend std::ostream& operator<<(std::ostream& os, const MemorySize& memSize) {
        static auto units = { "B", "KB", "MB", "GB", "TB", "PB" };
        double s = sCast<double>(memSize.size);
        for (auto u = units.begin(); u != units.end(); ++u, s /= Base) {
            if (s < Base) { return os << s << *u; }
        }
        return os;
    }

    Unit size;
};

struct MemoryUsage {
    MemorySize physicalMemory;
    MemorySize virtualMemory;
};

MemoryUsage memoryUsage();
MemoryUsage peakMemoryUsage();


struct ColorStr {
    enum CmdColor {
        Reset = 0,
        // foreground color.
        BlackFG = 30, RedFG, GreenFG, YellowFG, BlueFG, MagentaFG, CyanFG, WhiteFG,
        BrightBlackFG = 90, BrightRedFG, BrightGreenFG, BrightYellowFG, BrightBlueFG, BrightMagentaFG, BrightCyanFG, BrightWhiteFG,
        // background color.
        BlackBG = 40, RedBG, GreenBG, YellowBG, BlueBG, MagentaBG, CyanBG, WhiteBG,
        BrightBlackBG = 100, BrightRedBG, BrightGreenBG, BrightYellowBG, BrightBlueBG, BrightMagentaBG, BrightCyanBG, BrightWhiteBG,
    };


    static std::string get(const std::string& s, CmdColor fcolor, CmdColor bcolor) {
        return colorStr(fcolor) + colorStr(bcolor) + s + "\x1b[0m";
    }
    static std::string get(const std::string& s, CmdColor color) {
        return colorStr(color) + s + "\x1b[0m";
    }

    friend std::ostream& operator<<(std::ostream& os, CmdColor color) {
        return os << "\x1b[" << sCast<int>(color) << "m";
    }


protected:
    static std::string colorStr(CmdColor color) { return "\x1b[" + std::to_string(color) + "m"; }
};

}
}


#endif // CN_HUST_GOAL_SYSTEM_SYSTEM_H
