//
//  Zombie_rec_nr.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Zombie_rec_nr__
#define __BenchmarkerRandomFromText__Zombie_rec_nr__

#include <iostream>
#include "Zombie_recursive.h"

class Zombie_rec_nr : protected Zombie_recursive
{
public:
	int searchrec_nr(int lo, int hi);
    Zombie_rec_nr(const std::string &key);
    virtual int search(const std::string &text);
};
#endif /* defined(__BenchmarkerRandomFromText__Zombie_rec_nr__) */
