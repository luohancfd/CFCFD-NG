// Taken from ...
// Andrei Alexandrescu (2003)
// Assertions, C/C++ User's Journal
// May 13, 2003
//
// Included in local collection by
// Rowan J. Gollan
// 30-Apr-2007
//
#ifndef DBC_ASSERT_HH
#define DBC_ASSERT_HH

#include <iostream>
#include <ctime>
#include <cassert>
#include <cctype>
#include <stdexcept>
#include <cstdlib>
using namespace std;

#define BREAK_HERE abort()

template <class E> class AsserterEx;

class Asserter
{
protected:
    const bool holds_;
    static bool ignoreAll_;

    enum Action { giveUp = 'A', ignore = 'I', lineIgnore = 'L', globalIgnore = 'G',
        throwUp = 'T', debug = 'D' };

    static Action AskUser(const char* file, int line, const char* msg)
    {
        cerr << "Assertion failed in file " << file << ", line " << line << '\n'
            << msg;
        return Action(toupper(cin.get()));
    }

    bool DoHandle(Action action, bool& lineFlag) const
    {
        switch (action)
        {
        case giveUp:
            abort();
        case lineIgnore:
            lineFlag = true;
            break;
        case globalIgnore:
            ignoreAll_ = true;
            break;
        case ignore:
            break;
        default:
            return false;
        }
        return true;
    }

public:
    Asserter(bool holds) : holds_(holds)
    {
    }

    // Added R.J.Gollan - to alleviate compiler warnings.
    virtual ~Asserter() {}

    virtual bool Handle(const char* file, int line, bool& lineFlag) const
    {
        if (holds_ || ignoreAll_) return true;
        return DoHandle(AskUser(file, line, 
            "Press 'A' to abort, 'I' to ignore, 'L' to always ignore this "
            "assertion, 'G' to ignore all assertions, anything else to debug. "),
            lineFlag);
    }
        
    template <class Condition>
    static Asserter Make(Condition flag)
    {
        return Asserter(flag != 0);
    }

    template <class E> 
    static AsserterEx<E> Make(bool flag, const char* msg)
    {
        return AsserterEx<E>(flag, msg);
    }
};

// Moved to .cxx file as this will be multiply included.
//bool Asserter::ignoreAll_;

template <class E>
class AsserterEx : public Asserter
{
    const char* const msg_;

public:
    AsserterEx(bool holds, const char* msg) 
        : Asserter(holds), msg_(msg)
    {
    }

    // Added R.J.Gollan - to alleviate compiler warnings.
    virtual ~AsserterEx() {}

    virtual bool Handle(const char* file, int line, bool& ignoreLine) const
    {
        if (holds_ || ignoreAll_) return true;
        const Action action = AskUser(file, line,             
            "Press 'A' to abort, 'I' to ignore, 'L' to always ignore this "
            "assertion, 'G' to ignore all assertions, "
            "'T' to throw an exception, anything else to debug. ");
        if (action == throwUp) throw E(msg_);
        return Asserter::DoHandle(action, ignoreLine);
    }
};

#ifndef DEBUG

#define ASSERT\
    if (true) ; else struct Local {\
        Local(const Asserter& info)\
        {\
            static bool ignore_;\
            if (!ignore_ && !info.Handle(__FILE__, __LINE__, ignore_))\
                BREAK_HERE;\
        }\
    } localAsserter = Asserter::Make

#else

#define ASSERT\
    if (false) ; else struct Local {\
        Local(const Asserter& info)\
        {\
            static bool ignore_;\
            if (!ignore_ && !info.Handle(__FILE__, __LINE__, ignore_))\
                BREAK_HERE;\
        }\
    } localAsserter = Asserter::Make

#endif


#endif // DBC_ASSERT_HH
