//
//  Zombie_iter_noshare.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "Zombie_iter_noshare.h"

Zombie_iter_noshare::Zombie_iter_noshare(const std::string &key) : Zombie_iterative(key) {
    // Intentionally empty.
}

int Zombie_iter_noshare::search(const std::string &T) {
    struct  {
        int first, second;
    } todo[32];
    int i, count=0;
    int tos=0;
    int n = T.length();
    int lo=0, hi=n-(m-1);
    int kdleft, kdright;
    int probe;
    
    if (n < m) {
        return count;
    }
    
    // Here things get really dirty:
#define EMPTY (tos==0)
#define POP (--tos)
#define TOP (todo[tos])
#define PUSH(x,y) ++tos; TOP.first=(x); TOP.second=(y)
    
    PUSH(0,INT_MAX);
    
    for(;;) {
        probe=(lo+hi)>>1;
        assert(probe == (lo+hi)/2);
        
        assert(lo<hi);
        assert(lo<=probe && probe<hi);
        
        for (i=0; i<m && P[i]==T[probe+i]; i++) {
            // Intentionally empty
        }
        if (i==m) {
            count++;
        }
        // Do shifts
        {
            kdleft = probe - shl[T[probe]];
            kdright = probe + shr[T[probe + m - 1]];
            if (lo < kdleft) {
                // Left is good, so enstack right, regardless if it's good.
                PUSH(kdright, hi);
                hi = kdleft;
            } else {
                // Left is empty.
                assert(lo >= kdleft);
                // ...consider using the right...
                if ((lo=kdright) >= hi) {
                    // Right is also bad...
                    assert(lo >= hi);
                    assert(!EMPTY);
                    while ((TOP.first) >= (TOP.second)) {
                        assert(!EMPTY);
                        POP;
                    }
                    if (TOP.second == INT_MAX) {
                        return count;
                    } else {
                        lo = TOP.first;
                        hi = TOP.second;
                        POP;
                    }
                }
            }
            assert(lo < hi);
        }
    }
#undef EMPTY
#undef POP
#undef TOP
#undef PUSH
}