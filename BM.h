//
//  BM.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__BM__
#define __BenchmarkerRandomFromText__BM__

#include <iostream>
#include <vector>

class BM
{
    std::string p;
    // We need to store the size of the alphabet to build the table.
    int alphabetSize;
    int pSize;
    // Tables to store the values.
    std::vector<int> delta1;
    std::vector<int> delta2;
    std::vector<int> f;
    int matchCount;
    
public:
    int getCount();
    void buildDelta1(const std::string &p);
    void preProcessDelta2Case1(const std::string &p);
    void preProcessDelta2Case2();
    void buildTables(const std::string &p);
    BM(const std::string &key);
    void match(const std::string &in);
};

#endif /* defined(__BenchmarkerRandomFromText__BM__) */