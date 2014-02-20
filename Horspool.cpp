//
//  Horspool.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "Horspool.h"

int Horspool::getCount()
{
    return matchCount;
}

// Initialize the table
void Horspool::buildDelta1(const std::string &p) {
    int k = 0;
    
    // When a character is encountered that does not occur in the pattern, we can safely skip ahead for the whole length of the pattern.
    for (k = 0; k < alphabetSize; k++)
        delta1[k] = pSize;
    
    // Determine how far to skip ahead in the search text.
    for (k = 0; k < pSize - 1; k++)
        delta1[p[k]] = pSize - k - 1;
    
    /*std::cout << "Occurances: " ;
     for (int j = 0; j < delta1.size(); j++) {
     std::cout << delta1[j] << " " ;
     }
     std::cout << std::endl;*/
}

// The constructor takes the keyword and builds the table.
// 128 alphabet size => Unicode.
Horspool::Horspool(const std::string &key) : p(key) {
    alphabetSize = 256;
    delta1.resize(alphabetSize);
    pSize = key.length();
    buildDelta1(key);
}

void Horspool::match(const std::string &in) {
    matchCount = 0;
    
    int k = pSize - 1;
    
    while (k < in.length())
    {
        int j = pSize - 1;
        int i = k;
        while ((j >= 0) && (in[i] == p[j])) {
            j--;
            i--;
        }
        if (j == -1)
        {
            //std::cout << "Yay! Match at " << (k - pSize + 1) << std::endl;
            matchCount++;
            k += delta1[in[k]];
        }
        else
        {
            k += delta1[in[k]];
        }
    }
    
    //std::cout << "Matches found with Horspool: " << matchCount << std::endl;
}