# ifndef __CLONEWARS_TYPES_H
# define __CLONEWARS_TYPES_H


/*
 ********************************************************
 Basic integer types used throughout
 ********************************************************
 */

#include <cstdint>


typedef uint_fast8_t uint8;
typedef uint_fast32_t uint32;
typedef int_fast32_t sint32;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;

/*
 For the following usage of the below code...
 ```
 struct A {
   MEMBER(int,x)
   MEMBER(double,y)
 };
 ```
 you can use x_ and y_ inside the class and x() and y() outside,
 and the usage inside is non-const, while usage outside is const.
 (From https://stackoverflow.com/a/14487414)
 */
#define MEMBER(T,x) \
private: T x##_;    \
    public: T const& x () const { return x##_; }




#endif
