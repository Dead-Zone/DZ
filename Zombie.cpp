//
//  Zombie.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "Zombie.h"

Zombie::Zombie(const std::string &key) {
    P = key;
    m = P.length();
    // This resembles the Horspool.
    int i;
    for(i=0;i<SIGMA;i++) {
        shl[i]=m-1;
        shr[i]=m;
    }
    for(i=m-1; i>0; --i) {
        shl[P[i]]=i;
    }
    for(i=0; i<m-1; i++) {
        shr[P[i]]=m-1-i;
    }
    
}