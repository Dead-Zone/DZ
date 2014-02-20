//
//  Horspool.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Horspool__
#define __BenchmarkerRandomFromText__Horspool__

#include <iostream>
#include <vector>

class Horspool
{
    std::string p;
    // We need to store the size of the alphabet to build the table.
    int alphabetSize;
    int pSize;
    // Tables to store the values.
    std::vector<int> delta1;
    int matchCount;

public:
    int getCount();
    void buildDelta1(const std::string &p);
    Horspool(const std::string &key);
    void match(const std::string &in);
};

#endif /* defined(__BenchmarkerRandomFromText__Horspool__) */
