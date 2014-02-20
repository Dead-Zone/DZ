//
//  BM.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "BM.h"

int BM::getCount()
{
    return matchCount;
}

// Build table delta1 for last occurrence of character in the pattern.
// Bad character preprocessing.
void BM::buildDelta1(const std::string &p) {
    for(size_t i = 0; i < alphabetSize; i++){
        delta1[i] = -1;
    }
    
    for(size_t i = 0; i < pSize; i++){
        delta1[p[i]] = i;
    }
    
    /*std::cout << "Occurances: " ;
     for (int j = 0; j < delta1.size(); j++) {
     std::cout << delta1[j] << " " ;
     }
     std::cout << std::endl;*/
}

// The matching suffix occurs somewhere else in the pattern.
// Preprocessing for the good-suffix heuristics.
void BM::preProcessDelta2Case1(const std::string &p)
{
    int i = pSize;
    int j = i + 1;
    f[i] = j;
    while (i > 0)
    {
        while ((j <= pSize) && (p[i - 1] != p[j - 1]))
        {
            if (delta2[j] == 0) {
                delta2[j] = j - i;
            }
            
            j = f[j];
        }
        i--;
        j--;
        f[i] = j;
    }
    
    /*std::cout << "f: " ;
     for (int j = 0; j < f.size(); j++) {
     std::cout << f[j] << " " ;
     }
     std::cout << std::endl;
     std::cout << "delta2: " ;
     for (int j = 0; j < delta2.size(); j++) {
     std::cout << delta2[j] << " " ;
     }
     std::cout << std::endl;*/
}

// Only a part of the matching suffix occurs at the beginning of the pattern => This part is a border of the pattern. The pattern can be shifted as far as its widest matching border allows. For each suffix the widest border of the pattern that is contained in that suffix is determined.
// Preprocessing for the good-suffix heuristics.
void BM::preProcessDelta2Case2()
{
    int i, j;
    // The starting position of the widest border of the pattern at all is stored in f[0].
    j = f[0];
    for (i = 0; i <= pSize; i++)
    {
        // When the suffix of the pattern becomes shorter than f[0], the algorithm continues with the next-wider border of the pattern, i.e. with f[j].
        if (delta2[i] == 0) {
            delta2[i] = j;
        }
        if (i == j) {
            j = f[j];
        }
    }
    /*std::cout << "f: " ;
     for (int j = 0; j < f.size(); j++) {
     std::cout << f[j] << " " ;
     }
     std::cout << std::endl;
     std::cout << "delta2: " ;
     for (int j = 0; j < delta2.size(); j++) {
     std::cout << delta2[j] << " " ;
     }
     std::cout << std::endl;*/
}

// Calls the functions to build tables.
void BM::buildTables(const std::string &p)
{
    buildDelta1(p);
    preProcessDelta2Case1(p);
    preProcessDelta2Case2();
}

// The constructor takes the keyword and builds the tables.
BM::BM(const std::string &key) : p(key)
{
    alphabetSize = 256;
    delta1.resize(alphabetSize);
    pSize = key.length();
    delta2.resize(pSize+1);
    f.resize(pSize+1);
    buildTables(key);
}

void BM::match(const std::string &in)
{
    matchCount = 0;
    int i = 0, j;
    while (i <= in.length() - pSize)
    {
        //std::cout << "Attempting match at " << i << std::endl;
        j = pSize - 1;
        // Compare each character in the text with each character in the pattern, from right to left.
        while ((j >= 0) && (p[j] == in[i + j])) {
            j--;
        }
        // Match => Shift according to how much widest border allows.
        if (j < 0)
        {
            //std::cout << "Yay! Match at " << i << std::endl;
            matchCount++;
            i += delta2[0];
        }
        // Mismatch => Shift by max of delta1 and delta2.
        else {
            i += std::max(delta2[j + 1], j - delta1[in[i + j]]);
        }
    }
    //std::cout << "Matches found with BM: " << matchCount << std::endl;
}