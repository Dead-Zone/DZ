//
//  Zombie_rec_nr.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//
#include <cassert>
#include "Zombie_rec_nr.h"

int Zombie_rec_nr::searchrec_nr(int lo, int hi) {
    int i;
    int count;
    int probe = (lo+hi)>>1;
    assert(probe == (lo+hi)/2);
    assert(lo <= probe);
    assert(probe < hi);
    for (i=0; i<m && P[i] == T[probe+i]; i++) {
        // Intentionally empty;
    }
    count = (i==m);
    {
        int kdleft = probe - shl[T[probe]];
        if (lo < kdleft) {
            count += searchrec_nr(lo, kdleft);
        }
    }
    {
        int kdright = probe + shr[T[probe+m-1]];
        if (kdright < hi) {
            count += searchrec_nr(kdright, hi);
        }
    }
    return count;
}

Zombie_rec_nr::Zombie_rec_nr(const std::string &key) : Zombie_recursive(key) {
    // Intentionally empty.
}

int Zombie_rec_nr::search(const std::string &text) {
    n = text.length();
    if (n < m) {
        return 0;
    }
    d = 0;
    T = text.c_str();
    return searchrec_nr(0, n-(m-1));
}