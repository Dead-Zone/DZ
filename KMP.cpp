//
//  KMP.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "KMP.h"

void KMP::buildTable(const std::string &p) {
    // table[0] is always -1 and table[1] is always 0 --> hardcode
    table[0] = -1;
    table[1] = 0;
    int i = 2;
    int j = 0;
    
    // While i less than length of p
    while (i < pSize) {
        // If there is a substring in p
        if (p[i - 1] == p[j]) {
            j = j + 1;
            table[i] = j;
            i = i + 1;
        }
        // There was a substring, but now it has stopped
        else if (j > 0) {
            j = table[j];
        }
        // We were not in a substring, use 0
        else {
            table[i] = 0;
            i = i + 1;
        }
    }
    // This is what the table looks like
    /*std::cout << "Next table: " ;
     for (int j = 0; j < table.size(); j++) {
     std::cout << table[j] << " " ;
     }
     std::cout << std::endl;*/
}

int KMP::getCount()
{
    return matchCount;
}

// The constructor takes the keyword and builds a table as well.
KMP::KMP(const std::string &key) : p(key) {
    pSize = key.length();
    table.resize(key.length());
    buildTable(key);
}

void KMP::match(const std::string &in) {
    int pPos = 0;
    int sPos = 0;
    matchCount = 0;
    // While we don’t exceed the length of the string we are pattern searching
    while (sPos + pPos < in.length()) {
        //std::cout << "Attempting match at " << sPos << std::endl;
        // Match characters in pattern and string
        if (p[pPos] == in[sPos + pPos]) {
            //std::cout << p[pPos] << " " << in[sPos + pPos] << std::endl;
            // Reached end of pattern => match found.
            if (pPos == pSize - 1) {
                matchCount++;
                //std::cout << "Yay! Match at " << sPos << std::endl;
                //new value = j + i - table[i]
                sPos = sPos + pPos - table[pPos];
                // We don’t need to compare some characters of p again,
                // We know that they already match up due to the table.
                // Shift according to table values if there is still space to pattern match
                if (sPos + pSize <= in.length()) {
                    if (table[pPos] > -1) {
                        pPos = table[pPos];
                    }
                    else {
                        pPos = 0;
                    }
                }
                else {
                    //std::cout << "Matches found with KMP: " << matchCount << std::endl;
                    //return;
                }
                
            }
            else {
                pPos++;
            }
        }
        // Characters don't match, shift pattern
        else {
            //new value = j + i - table[i]
            sPos = sPos + pPos - table[pPos];
            // We don’t need to compare some characters of p again,
            // We know that they already match up due to the table.
            if (table[pPos] > -1) {
                pPos = table[pPos];
            }
            else {
                pPos = 0;
            }
        }
    }
    //std::cout << "Matches found with KMP: " << matchCount << std::endl;
}