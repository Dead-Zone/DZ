#include <iostream>
#include<fstream>
#include <string>
#include <vector>
#include <stack>
#include <algorithm>
#include <cassert>
#include <stdint.h> /* for uint64_t */
#include <time.h>  /* for struct timespec */
#include <mach/mach_time.h>
#include <utility>

#include "RandomGenerator.h"

#include "BM.h"
#include "Horspool.h"
#include "KMP.h"

// DZ includes
#include "Zombie_iter_noshare.h"
#include "Zombie_recursive.h"
#include "Zombie_rec_nr.h"

int main (int argc, char * const argv[]) {
    RandomGenerator rg;
    srand(time(NULL));
    
    std::cout << "" ;
    system("pwd");
    
    std::fstream myFile;
    myFile.open("bible.txt");
    std::string line;
    std::string input("");
    if (myFile.is_open())
    {
        while ( myFile.good() )
        {
            getline (myFile,line);
            input += line;
        }
        myFile.close();
    }
    else
    {
        std::cout << "Unable to open file" << std::endl;
        return 0;
    }
    
    for (int cnt = 0; cnt < 120; cnt++)
    {
        int patternSize = 4;
        while (patternSize <= 16384)
        {
            
            // Generate random index in file to get pattern from
            int num = rg.generateRandomIndex(input.length() - patternSize);
            std::string temp = input.substr(num, patternSize);
            std::string keyword(temp);
            
            //std::string input("abracadabrabra abracad abracadabrabracadabrabracadabra brabra abracadabra");
            //std::string keyword("abra");
            // There are 13 instances of the pattern "abra".
            
            for (int rep = 0; rep < 30; rep++)
            {
                
                //std::cout << "Input string = \"" << input << "\"" << std::endl;
                std::cout << cnt + 1 << "\t" ;
                std::cout << rep + 1 << "\t" ;
                
                std::cout << keyword << "\t" ;
                
                std::cout << patternSize << "\t" ;
                
                mach_timebase_info_data_t info;
                mach_timebase_info(&info);
                
                /* BM */
                
                BM bm(keyword);
                uint64_t start = mach_absolute_time();
                //std::cout << start << std::endl;
                bm.match(input);
                uint64_t stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                uint64_t duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t" ;
                
                /* Horspool */
                
                Horspool h(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                h.match(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t" ;
                
                /* KMP */
                
                KMP kmp(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                kmp.match(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t" ;
                
                // Zombies
                
                Zombie_recursive zrec(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                int zrec_count = zrec.search(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t";
                
                Zombie_rec_nr zrec_nr(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                int zrec_nr_count = zrec_nr.search(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t";
                
                Zombie_iterative ziter(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                int ziter_count = ziter.search(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t";
                
                Zombie_iter_noshare ziter_noshare(keyword);
                // Do some matching.
                start = mach_absolute_time();
                //std::cout << start << std::endl;
                int ziter_noshare_count = ziter_noshare.search(input);
                stop = mach_absolute_time();
                //std::cout << stop << std::endl;
                
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << std::endl;
                
                
                int counts[2100];
                int counter = 0;
                counts[counter++] = bm.getCount();
                counts[counter++] = h.getCount();
                counts[counter++] = kmp.getCount();
                
                counts[counter++] = zrec_count;
                counts[counter++] = zrec_nr_count;
                counts[counter++] = ziter_count;
                counts[counter++] = ziter_noshare_count;
            }
            
            patternSize *= 2;
        }
    }
    
	return 0;
}
