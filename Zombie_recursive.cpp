//
//  Zombie_recursive.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include <cassert>
#include "Zombie_recursive.h"

int Zombie_recursive::searchrec(int lo, int hi)
{
    if (lo >= hi) {
        d = lo;
        return 0;
    } else {
        int i;
        int count = 0;
        int probe = (lo+hi)>>1;
        assert(probe == (lo+hi)/2);
        assert(lo <= probe);
        assert(probe < hi);
        for (i=0; i<m && P[i] == T[probe+i]; i++) {
            // Intentionally empty;
        }
        count = (i==m);
        count += searchrec(lo, probe - shl[T[probe]]);
        count += searchrec(std::max(d, probe + shr[T[probe+m-1]]), hi);
        return count;
    }
}

Zombie_recursive::Zombie_recursive(const std::string &key) : Zombie(key), d(0), n(0), T(0) {
    // Intentionally empty.
}
int Zombie_recursive::search(const std::string &text) {
    d = 0;
    n = text.length();
    T = text.c_str();
    return searchrec(0, n-(m-1));
}