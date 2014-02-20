//
//  Zombie_iterative.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include <cassert>
#include "Zombie_iterative.h"

Zombie_iterative::Zombie_iterative(const std::string &key) : Zombie(key) {
		// Intentionally empty.
	}

int Zombie_iterative::search(const std::string &T) {
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
					lo = kdright;
					if (kdright >= hi) {
						// Right is also bad...
						//assert(kdright >= high);
						while (!EMPTY && !((lo=(std::max(TOP.first, lo))) < (hi=TOP.second))) {
							assert(!EMPTY);
							POP;
						}
						if (EMPTY) {
							return count;
						} else {
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