//
//  Zombie.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Zombie__
#define __BenchmarkerRandomFromText__Zombie__

#include <iostream>

class Zombie {
protected:
	// Make this protected so that the derived (inheriting) classes can access them.
	static const int SIGMA=256;
	int shl[SIGMA], shr[SIGMA];
	std::string P;
	int m;
    
public:
    Zombie(const std::string &key);
    virtual int search(const std::string &text) = 0;
};

#endif /* defined(__BenchmarkerRandomFromText__Zombie__) */
