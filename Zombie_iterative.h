//
//  Zombie_iterative.h
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#ifndef __BenchmarkerRandomFromText__Zombie_iterative__
#define __BenchmarkerRandomFromText__Zombie_iterative__

#include <iostream>
#include "Zombie.h"

class Zombie_iterative : public Zombie
{
public:
    Zombie_iterative(const std::string &key);
    virtual int search(const std::string &T);
};
#endif /* defined(__BenchmarkerRandomFromText__Zombie_iterative__) */
