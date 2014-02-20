//
//  Zombie_recursive.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Zombie_recursive__
#define __BenchmarkerRandomFromText__Zombie_recursive__

#include <iostream>
#include "Zombie.h"

class Zombie_recursive : public Zombie
{
protected:
	int d, n;
	const char *T;
    
public:
    int searchrec(int lo, int hi);
	Zombie_recursive(const std::string &key);
    virtual int search(const std::string &text);
};

#endif /* defined(__BenchmarkerRandomFromText__Zombie_recursive__) */
