//
//  RandomGenerator.cpp
//  BenchmarkerRandomFromText
//
//  Created by Melanie on 2014/02/20.
//
//

#include "RandomGenerator.h"

RandomGenerator::RandomGenerator()
{}

int RandomGenerator::generateRandomIndex(int maxLength)
{
    return rand() % maxLength;
}