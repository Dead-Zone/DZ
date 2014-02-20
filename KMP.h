//
//  KMP.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__KMP__
#define __BenchmarkerRandomFromText__KMP__

#include <iostream>
#include <vector>

class KMP
{
	std::string p;
	// We need to store the size of the pattern for table computation.
    int pSize;
    // Table to store the next values.
    std::vector<int> table;
    int matchCount;
    
public:
    void buildTable(const std::string &p);
    int getCount();
    KMP(const std::string &key);
    void match(const std::string &in);
};

#endif /* defined(__BenchmarkerRandomFromText__KMP__) */
