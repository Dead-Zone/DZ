//
//  Zombie_iter_noshare.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Zombie_iter_noshare__
#define __BenchmarkerRandomFromText__Zombie_iter_noshare__

#include <iostream>
#include <cassert>
#include "Zombie_iterative.h"

class Zombie_iter_noshare : protected Zombie_iterative
{
public:
    Zombie_iter_noshare(const std::string &key);
    virtual int search(const std::string &T);
    
};
#endif /* defined(__BenchmarkerRandomFromText__Zombie_iter_noshare__) */
