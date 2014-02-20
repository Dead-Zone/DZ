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

#include "BM.h"
#include "Horspool.h"
#include "KMP.h"

// DZ includes
#include "Zombie_iter_noshare.h"
#include "Zombie_recursive.h"
#include "Zombie_rec_nr.h"

/************
 RANDOM
 ***********/

class RandomGenerator
{
public:
    RandomGenerator()
    {}
    
    int generateRandomIndex(int maxLength)
    {
        return rand() % maxLength;
    }
    
};

/*
 *  main.cpp
 *  Deadzoner
 *	Give a short example of how to do matching using a safe shifter.
 *
 *  Created by Bruce W. Watson on 2012/02/17.
 *  Copyright 2012 FASTAR. All rights reserved.
 *
 *	TODO: Get this into C, and use the interface in the SMART environment.
 *	TODO: Remove hard-coding, and either get something from the command-line, or standard input, etc.
 *	TODO: Take a good look at whether to change over to C-style strings?
 *	TODO: We probably want to introduce "match orders" (aka "mo") into this.
 *		DONE!!
 *	TODO: Let's think about how/if/when to split this into several files.
 *	TODO: Tune the hell out of it.
 *	TOOD: Urgently add some other shifters, perhaps even some dumb ones.
 *	TODO: Very urgently pass information back up the recursion tree.
 *		DONE!!
 *	TODO: See other TODO's below.
 *	NB: Zones (live and dead) here are in [) format.
 */

// Need this for pair<>.

// Here are a few match orders. Note that they have the same interface (NB!!!).
// They always take the keyword p as an argument, but don't always need it.
// Super simple left-to-right match order.
struct MO_fwd {
	// Define this to make it look like a function. This is a common C++ idiom.
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// Don't use argument p, since it's not needed.
		return i;
	}
};

// Somewhat simple right-to-left match order.
struct MO_rev {
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// Here we really do need p, since we're indexing from right to left, so we need p.length.
		return p.length() - 1 - i;
	}
};

// Here's a more complex outside-in-left-to-right.
struct MO_OILR {
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// For even numbers we index from the left, for odd from the right. Here we count on integer arithmetic in C/C++ rounding down.
		// NB: we could've really done this with that cool ?: idiom, but that's only for nerds who don't understand that compilers do the same thing, so let's stick with readability.
		if (i % 2 == 0) {
			return i / 2;
		} else {
			return p.length() - i / 2 - 1;
		}
	}
};

// Here's a more complex outside-in-right-to-left.
struct MO_OIRL {
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// This looks a lot like MO_OILR because it's the mirror image.
		if (i % 2 == 0) {
			return p.length() - i / 2 - 1;
		} else {
			return i / 2;
		}
	}
};

// Here's a really freaky one where we index in lock-step in the first half (l-to-r) and same in the second half.
struct MO_lockstep_LR {
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// Use evenness to decide whether we're doing the first half or the second.
		if (i % 2 == 0) {
			return i / 2;
		} else {
			return (p.length() + i) / 2;
		}
	}
};

// Here's an inside out left-to-right.
struct MO_IOLR {
	inline int operator()(int i, const std::string &p) const {
		//assert(i >= 0 && i < p.length());
		// Use evenness to decide whether we go left or right.
		if (i % 2 == 0) {
			return (p.length() / 2) - ((i + 1) / 2);
		} else {
			return (p.length() / 2) + ((i + 1) / 2);
		}
	}
};

// Define a shifter class, which is constructed from the keyword, and gives a safe left- and right-shift distance after a match attempt.
// TODO: there will be many kinds of shifters, some of which take more information (e.g. mismatching character) to make a shift. See SPARE Parts, or any of the taxonomies for some ideas.
template<typename MOrder>
class Shifter_kmp_bm {
private:
	// Keep the keyword length for speed and so we can use it in assertions later during shifts.
	// NB: This must come before declarations of the other fields, as the constructor's field initialization depends on this order.
	int plen;
	// TODO: these really should be vectors or proper arrays, but leave it dirty for now.
	// These must be arrays indexable by [0,|p|] to give the shift distance. NB: this is one _longer_ than |p|, as the subscript is the index of the first mismatch, which is |p| if p matched.
	int *shl, *shr;
	// NB: we don't keep a copy of the keyword, since we're just responsible for shifting...not matching.
	// Make sure that Shifter_kmp_bm's are not assignable by putting the assignment operator and copy-constructor here...also a C++ idiom.
	Shifter_kmp_bm(const Shifter_kmp_bm &);
	const Shifter_kmp_bm &operator=(const Shifter_kmp_bm &);
public:
	// Construct a shifter (a pair of shift tables, actually).
	Shifter_kmp_bm(const std::string &p, const MOrder &mo) : plen(p.length()), shl(new int[plen+1]), shr(new int[plen+1]) {
		// Build the shift tables...the memory is already allocated in the shl, shr initializer-list above.
		for (int i(0); i <= plen; ++i) {
			// What do we know after we've gotten a mismatch at p[i]? We need to do a shift left by shl[i] and right by shr[i].
			// We know that p[0..i) matched, and that p[i] mismatched.
            
			// NB: shl[i] >= shl[i-1] because more information (ie. more progress in a match) is always better and gives longer shifts. Same applies to shr!
			// This is used below in linear searches to define the starting point (of the linear search) in a smarter way.
			
			// So, how far can we safely shift to the left? At least one symbol...but we can shift more:
			// smallest k>0 such that
			//	p[0..i) (the part already matched) aligns with the left-moving p, at p[k..i+k) (while taking care of the bounds)...here, I didn't take care of writing out mo()
			// Now we can write this down using the match-order as well:
			//	written as (\forall n : 0 <= n < i : mo(n,p)+k < plen => p[mo(n,p)] == p[mo(n,p)+k])
			//	(The funny stuff on the left of the implication (in ascii-art as =>...which is not a funny comparison symbol!!!) is the bounds-check...cf. "while taking care of bounds", above.)
			// Do this with a linear-search to find the smallest k>0 such that... Recall we start with shl[i-1]
			for (int k(i == 0 ? 1 : shl[i-1]); ; ++k) {
				// Declare n here so we can use it after the for-loop.
				int n;
				// Looking at the above \forall, and recalling the definition of implication, we know that A => B \equiv !(A && !B) \equiv !A || B (thanks to the definition of =>).
				// Copying from the above \forall, we have (mo(n,p)+k < plen => p[mo(n,p)] == p[mo(n,p)+k]), which is equivalent to (mo(n,p)+k >= plen || p[mo(n,p)] == p[mo(n,p)+k])
				for (n = 0; n < i && (mo(n,p)+k >= plen || p[mo(n,p)] == p[mo(n,p)+k]); ++n) {
					// Intentionally empty
				}
				// If n == i, we have a successful alignment.
				if (n == i) {
					shl[i] = k;
					break;
				}
			}
			// Since we exited the loop, we have a successful alignment at shift k.
			
			// Now let's do the right-shift:
			// How far can we safely shift to the right? Also at least one symbol...but, again, we can shift more:
			// smallest k>0 such that
			//	p[0..i) (the part already matched) aligns with the right-moving p, at p[0-k..i-k) (which is of course p[0..i-k) because of bounds-checking).
			//	written as (\forall n : 0 <= n < i : mo(n,p)-k >= 0 => p[mo(n,p)] == p[mo(n,p)-k])
			// The rationale for this is above. Do this also with a linear-search.
			for (int k(i == 0 ? 1 : shr[i-1]); ; ++k) {
				int n;
				// Copying from the above \forall, we have (mo(n,p)-k >= 0 => p[mo(n,p)] = p[mo(n,p)-k]), which is equivalent to mo(n,p)-k < 0 || p[mo(n,p)] == p[mo(n,p)-k]
				for (n = 0; n < i && (mo(n,p)-k < 0 || p[mo(n,p)] == p[mo(n,p)-k]); ++n) {
					// Intentionally empty
				}
				// If n == i, we have a successful alignment.
				if (n == i) {
					shr[i] = k;
					break;
				}
			}
			// Since we exited the loop, we have a successful alignment at shift k.
            /*#ifndef NDEBUG
             std::cerr << "i = " << i << "\tshl[i] = " << shl[i] << "\tshr[i] = " << shr[i] << std::endl;
             #endif*/
		}
        /*#ifndef NDEBUG
         std::cerr << std::endl;
         #endif*/
	}
	// Destructor to get rid of the tables.
	~Shifter_kmp_bm() {
		// Free up the heap-allocated memory, using the funny form of delete [] because they were allocated using the array form of "new".
		delete [] shl;
		delete [] shr;
	}
	// The two shift functions, based on the keyword. Note that we give them some additional information, such as where in the input string the match attempt was.
	int shift_left(int i, const std::string &in, int j) const {
		//assert(i >= 0 && i <= plen && j >= 0 && j < in.length());
		return shl[i];
	}
	int shift_right(int i, const std::string &in, int j) const {
		//assert(i >= 0 && i <= plen && j >= 0 && j < in.length());
		return shr[i];
	}
};

// Here's a super dumb shifter which doesn't know anything, and just does the minimal shift of one symbol always.
class Shifter_naive {
public:
	// Construct a shifter, which is actually empty because it's just dumb. It's forced to take a MOrder (a match order) as constructor argument because the other shifters do, and they share a common interface.
	template<typename MOrder>
	Shifter_naive(const std::string &p, const MOrder &mo) {
		// Intentionally empty.
	}
	// The two shift functions, always just shifting by one.
	int shift_left(int i, const std::string &in, int j) const {
		//assert(j >= 0 && j < in.length());
		return 1;
	}
	int shift_right(int i, const std::string &in, int j) const {
		//assert(j >= 0 && j < in.length());
		return 1;
	}
};

// Here's a shifter depending on only the character aligned with p[0] and p[p.length()-1].
// TODO: Normally, we would build and maintain a pair of tables here, one to be indexed by the input char aligned with p[0], resp. p[p.length()-1].
//	To be honest, I was lazy, and wanted to get this going soon!
class Shifter_end_chars {
private:
	// Keep a local copy of p for future shifts.
	std::string pcopy;
public:
	template<typename MOrder>
	Shifter_end_chars(const std::string &p, const MOrder &mo) : pcopy(p) {
		// Intentionally empty.
	}
	int shift_left(int i, const std::string &in, int j) const {
		//assert(j >= 0 && j < in.length());
		// Take a look at in[j] (which is aligned with p[0]) and see where it can also fit in p.
		// Do this with a linear-search using k, declared here because we need it afterwards.
		int k;
		for (k = 1; k < pcopy.length() && in[j] != pcopy[k]; ++k) {
			// Intentionally empty.
		}
		// We've either found a next available match of in[j] in p, or we've run off the end of p, in which case we shift the full p.length().
		return k;
	}
	int shift_right(int i, const std::string &in, int j) const {
		//assert(j >= 0 && j < in.length());
		// This obviously resembles shift_left.
		int k;
		for (k = 1; k < pcopy.length() && in[j + pcopy.length() - 1] != pcopy[pcopy.length() - 1 - k]; ++k) {
			// Intentionally empty.
		}
		return k;
	}
};

// Here's a shifter depending on only the character aligned with p[0] and p[p.length()-1].
class Shifter_end_chars_cached {
private:
	// Keep a pair of shift tables (one left and one right). These will be indexed by the chars aligned with the first (resp. last) chars of the keyword.
	// This is also hideous...I should have used a vector<int> here, but let's just get this going soon.
	int *left, *right;
	// This is hideous style...I really should just get this from limits.h:
	static const int lastchar = 0xff;
	const int plen;
public:
	template<typename MOrder>
	Shifter_end_chars_cached(const std::string &p, const MOrder &mo) : left(new int[lastchar+1]), right(new int[lastchar+1]), plen(p.length()) {
		// Tables allocated, so now populate them, indexed by the character.
		// First set the shifts to p.length() as the default.
		for (int k(0); k <= lastchar; k++) {
			left[k] = right[k] = plen;
		}
		// Now go right to left through p and set the left-shifter. Don't go down to 0, as the left-most character in p doesn't contribute to the shift.
		// See my thesis or other toolkit (SPARE or Eindhoven Pattern Kit) docs for more.
		for (int j(plen-1); j > 0; --j) {
			left[/*reinterpret_cast<unsigned char> */(p[j])] = j;
		}
		// Now do the same thing left to right to get the right shifts...see explanation above.
		for (int j(0); j < plen-1; ++j) {
			right[/*reinterpret_cast<unsigned char> */(p[j])] = plen - 1 - j;
		}
	}
	virtual ~Shifter_end_chars_cached() {
		delete [] left;
		delete [] right;
	}
	// We've allocated some memory using new, so need to delete it too.
	int shift_left(int i, const std::string &in, int j) const {
        //		assert(j >= 0 && j < in.length());
		return left[/*reinterpret_cast<unsigned char> */(in[j])];
	}
	int shift_right(int i, const std::string &in, int j) const {
        //		assert(j >= 0 && j < in.length());
		return right[/*reinterpret_cast<unsigned char> */(in[j + plen - 1])];
	}
};

// This class combines the best of two shifters using max.
template<typename Sh1, typename Sh2>
class Shifter_combiner {
private:
	// Keep local copies of the two shifters.
	Sh1 sh1;
	Sh2 sh2;
public:
	template<typename MOrder>
	Shifter_combiner(const std::string &p, const MOrder &mo) : sh1(p, mo), sh2(p, mo) {
		// Intentionally empty.
	}
	// The two shift functions, always just maximizing the shifts
	int shift_left(int i, const std::string &in, int j) const {
		assert(j >= 0 && j < in.length());
		return std::max(sh1.shift_left(i, in, j), sh2.shift_left(i, in, j));
	}
	int shift_right(int i, const std::string &in, int j) const {
		assert(j >= 0 && j < in.length());
		return std::max(sh1.shift_right(i, in, j), sh2.shift_right(i, in, j));
	}
};

// Here's another function object (see the "match orders" above for what that means in C++) to choose the probe point in [low, high).
struct Probe_chooser_mid {
	inline int operator()(int low, int high) const {
		//assert(low < high);
		// NB: it's important that integer arithmetic rounds down (as it does in C/C++), in case we've got a live-zone of one-index (rounding up would be evil).
		return (low + high) /2;
	}
};

// Here's another probe chooser which goes a third of the way from low to high.
struct Probe_chooser_first_third {
	inline int operator()(int low, int high) const {
		//assert(low < high);
		// NB: it's important that integer arithmetic rounds down (as it does in C/C++), in case we've got a live-zone of one-index (rounding up would be evil).
		return low + (high - low) / 3;
	}
};

// Here's a dumb probe chooser, which always chooses the low...thereby giving traditional left-to-right algorithms (over the whole input).
struct Probe_chooser_LR {
	inline int operator()(int low, int high) const {
		//assert(low < high);
		return low;
	}
};

// And another dumb probe chooser, which always chooses the high-1...thereby giving right-to-left algorithms (over the whole input).
struct Probe_chooser_RL {
	inline int operator()(int low, int high) const {
		//assert(low < high);
		return high - 1;
	}
};

// For a given keyword, we have a "pattern matcher" object, which embodies all the shift knowledge for that keyword.
// It's template parameterized by the match order (for performance), and also by the Shifter. These two types have defaults as well, in case you don't know what to choose.
template<typename MOrder = MO_rev, typename Shifter = Shifter_kmp_bm<MOrder>, typename Probe_chooser = Probe_chooser_mid>
class Pattern_matcher {
private:
	// Ranges are a recurring theme, so us a typedef to give it a short name.
	typedef std::pair<int, int> index_range_t;
	
	// We need to keep a copy of the keyword being matched.
	// TODO: do we want to have a separate "match attempter" object?
	std::string p;
	// Keep a local copy of the match order object.
	// NB: !!! It's critical that this is declared before the Shifter, since mo is used in the initialization of the Shifter.
	MOrder mo;
	
	// We need a shifter object, which is built around the same match order, otherwise this won't work at all!!
	Shifter sh;
	
	// Let's keep an attempt counter here, and make it mutable since it's a trivial data-member and not related to deep functionality.
	mutable int attempt_count;
	
	// Helper function to do the recursive matching. Give it the same name, thus overloading it (the compiler chooses based on the arguments).
	// Argument explanation:
	//	in is the input string
	//	[low, high) is the live-zone to explore. NB: the upperbound is NOT inclusive. TODO: perhaps bundle these together in an index_range_t.
	//	out is the output stream to print debugging info to
	//	debug_indent is how far to indent the debugging messages.
	// Return value explanation:
	//	Return [l, h) of a zone now known to be dead. NB: this will contain the [low, high) passed in.
	// Note that this function is const, meaning it doesn't change the Pattern_matcher object's fields.
	index_range_t match(const std::string &in, int low, int high, std::ostream &out, const std::string &debug_indent = "") const {
		// Prepare the return value of known-dead.
		index_range_t known_dead;
        
		// If we've got an empty live-range, then just return.
		if (high <= low) {
            /*#ifndef NDEBUG
             //			out << debug_indent << "\tEmpty range, so doing nothing and returning right away. This doesn't count as an attempt!" << std::endl << std::endl;
             #endif*/
			// Return whatever was passed in, as we know for sure that's a dead-zone.
			known_dead.first = low;
			known_dead.second = high;
			return known_dead;
		}
        
        /*#ifndef NDEBUG
         out << std::endl << debug_indent << "Invoked with a live-zone of [" << low << ',' << high << ')' << std::endl;
         #endif*/
		
		// We know we have at least one place to check now.
		//assert(low < high);
		++attempt_count;
        
		// Our aim now is to check all of [low, high), reporting matches as we go.
		// Let's start with some point in this live-zone, e.g. the middle, but let's determine it using our Probe_chooser function object.
		int probe(Probe_chooser()(low, high));
		
        /*#ifndef NDEBUG
         // Do some pretty printing here, giving a picture of the input string with the keyword aligned at the right place.
         out << debug_indent << "Attempting a match at " << probe << std::endl;
         out << debug_indent << "Picture" << std::endl;
         out << debug_indent << '\t' << in << std::endl << debug_indent << '\t';
         for (int z(0); z < probe; ++z) {
         out << ' ';
         }
         out << p << std::endl;
         #endif*/
		
		// We'll need an index into p for matching. Declare it here since we'll need it after the for-loop here.
		int i;
		// Here, we make use of the "match-order" function "mo", which is just a permutation on [0, p.length()).
		// Keep going until we have a mismatch:
		for (i = 0; i < p.length() && p[mo(i,p)] == in[probe+mo(i,p)]; ++i) {
			// Intentionally empty.
		}
        /*#ifndef NDEBUG
         out << debug_indent << "Match got as far as i = " << i << std::endl;
         #endif*/
		// Report any match, which is the case if we got all the way to i == p.length().
		if (i == p.length()) {
			//out << debug_indent << "Yay - match starting at index " << probe << std::endl;
		}
        
		// Based on that match attempt (and where it went wrong), what do we now know, in terms of live/dead, etc.?
		//	We have tried at position probe, so that's now dead...ie. [probe, probe] is a dead-zone.
		//	We can shift left (from probe) by shl(i) to get the next possible attempt position to the left. That means that (probe - shl(i), probe] is a dead-zone.
		//	We can shift right (from probe as well) by shr(i) to get the next possible attempt position to the right, meaning that [probe, probe + shr(i)) is also dead.
		//	Taken together, we know that (probe - shl(i), probe + shr(i)) is dead...or in our [) form: [probe - shl(i) + 1, probe + shr(i)) is dead.
		// NB: I've written shl(i) here (same for shr), but really mean sh.shift_left(i, in, probe).
		// Let's update our known-dead:
		known_dead.first = probe - sh.shift_left(i, in, probe) + 1;
		known_dead.second = probe + sh.shift_right(i, in, probe);
		// What remains for us to check? Well, we started with [low, high) and now we have the known_dead, which might split [low, high)
		//	Left of probe, we should go and check in [low, known_dead.first)
		//	Right of probe, we should go and check in [known_dead.second, high)
		//		More on this later, because we may gain some info after having gone left.
		
        /*#ifndef NDEBUG
         // Let's output those shifts for debugging purposes.
         out << debug_indent << "Will now shift left/right by " << sh.shift_left(i, in, probe) << '/' << sh.shift_right(i, in, probe) << std::endl;
         out << debug_indent << "Known dead-zone is [" << known_dead.first << ',' << known_dead.second << ')' << std::endl;
         #endif
         #ifndef NDEBUG
         out << debug_indent << "Left will be [" << low << ',' << known_dead.first << ") and we think (for now) right will be [" << known_dead.second << ',' << high << ')' << std::endl;
         #endif
         */
		// As we go left first with a recursive invocation. That invocation will give us back its own "definitely_dead":
		index_range_t left_known_dead;
		left_known_dead = match(in, low, known_dead.first, out, debug_indent/* + "\t\t"*/);
		// What do we now know? There's a chance that this invocation passed back a dead-zone which is bigger than the zone we passed into it. As a result, merge zones:
		known_dead.first = std::min(known_dead.first, std::min(low, left_known_dead.first));
		known_dead.second = std::max(known_dead.second, left_known_dead.second);
		// NB: we may have shifted known_dead.second, which is what we're about to pass into our second recursive invocation (see above) to check live zone [known_dead.second, high), which is really cool.
		
		// Again, let's collect what the recursive invocation tells us, since it may be useful information.
		index_range_t right_known_dead;
		right_known_dead = match(in, known_dead.second, high, out, debug_indent/* + "\t\t"*/);
		// What do we now know? As before, there's a chance that this invocation passed back a dead-zone which is bigger than the zone we passed into it. As a result, merge zones:
		known_dead.first = std::min(known_dead.first, right_known_dead.first);
		known_dead.second = std::max(known_dead.second, std::max(high, right_known_dead.second));
		
        /*#ifndef NDEBUG
         out << debug_indent << "Guaranteed dead-zone = [" << known_dead.first << ',' << known_dead.second << ')' << std::endl;
         #endif*/
		return known_dead;
	}
public:
	// The constructor takes the keyword and builds a shifter as well.
	Pattern_matcher(const std::string &key) : p(key), sh(key, mo), attempt_count(0) {
		// Intentionally empty...all fields initialized in the initializer list.
	}
	// Do some matching.
	// TODO: see if we really want to register matches via an ostream. In reality, we'd probably want to do some array or a counter, etc.
	void match(const std::string &in, std::ostream &out) const {
		// Set our attempt count to 0, since this is a "per match" attempt counter.
		attempt_count = 0;
		// The initial live zone isn't [0,in.length()), but a bit less on the upper side, as the last p.length()-1 symbols cannot possibly start a match.
		// Recall that zones are in the [) format.
		match(in, 0, in.length()-(p.length()-1), out);
        //#ifndef NDEBUG
		//out << "Total attempts = " << attempt_count << std::endl << std::endl;
        //#endif
	}
};

// This (template) class has the same interface as Pattern_matcher (more or less), but is implemented with iteration instead of recursion.
// For details of the template parameters, see the above.
template<typename MOrder = MO_rev, typename Shifter = Shifter_kmp_bm<MOrder>, typename Probe_chooser = Probe_chooser_mid>
class Pattern_matcher_iterative {
private:
	// Ranges are a recurring theme, so us a typedef to give it a short name.
	typedef std::pair<int, int> index_range_t;
	
	// We need to keep a copy of the keyword being matched.
	// TODO: do we want to have a separate "match attempter" object?
	std::string p;
	// Keep a local copy of the match order object.
	// NB: !!! It's critical that this is declared before the Shifter, since mo is used in the initialization of the Shifter (in the constructor).
	MOrder mo;
	
	// We need a shifter object, which is built around the same match order, otherwise this won't work at all!!
	Shifter sh;
	
public:
	// The constructor takes the keyword and builds a shifter as well.
	Pattern_matcher_iterative(const std::string &key) : p(key), sh(key, mo) {
		// Intentionally empty...all fields initialized in the initializer list.
	}
	
	// Do some matching and report how many matches were found.
	// TODO: see if we really want to register matches via an ostream. In reality, we'd probably want to do some array or a counter, etc.
	int match(const std::string &in, std::ostream &out) const {
		// Set our attempt count to 0, since this is a "per match" attempt counter.
		int match_count(0);
		
		// Also keep track of the part which is known to be dead...[0,d)...initially [0,0)
		int d(0);
		
		// Use a local stack to keep track of the ranges that still need to be examined:
		std::stack<index_range_t> todo;
		// Load it with the initial range, leaving out the stuff at the end, since it can't contribute to a match:
		todo.push(index_range_t(0, in.length()-(p.length()-1)));
		
		while (!todo.empty()) {
			// Get the next one to do...
			index_range_t curr(todo.top());
			// ...give them nice names...
			int low(curr.first), high(curr.second);
			// ...and remove it from the stack.
			todo.pop();
			// Adjust this element according to the already-known-dead [0,d)
			low = std::max(low, d);
			
			// If we've got an empty live-range, then just return.
			if (high <= low) {
#ifndef NDEBUG
				//			out << debug_indent << "\tEmpty range, so doing nothing and returning right away. This doesn't count as an attempt!" << std::endl << std::endl;
#endif
				d = low;
				// Intentionally empty.
			} else {
#ifndef NDEBUG
				out << std::endl << "Invoked with a live-zone of [" << low << ',' << high << ')' << std::endl;
#endif
				// We know we have at least one place to check now.
				assert(low < high);
				
				// Our aim now is to check all of [low, high), reporting matches as we go.
				// Let's start with some point in this live-zone, e.g. the middle, but let's determine it using our Probe_chooser function object.
				int probe(Probe_chooser()(low, high));
#ifndef NDEBUG
				// Do some pretty printing here, giving a picture of the input string with the keyword aligned at the right place.
				out << "Attempting a match at " << probe << std::endl;
				out << "Picture" << std::endl;
				out << '\t' << in << std::endl << '\t';
				for (int z(0); z < probe; ++z) {
					out << ' ';
				}
				out << p << std::endl;
#endif
				
				// We'll need an index into p for matching. Declare it here since we'll need it after the for-loop here.
				int i;
				// Here, we make use of the "match-order" function "mo", which is just a permutation on [0, p.length()).
				// Keep going until we have a mismatch:
				for (i = 0; i < p.length() && p[mo(i,p)] == in[probe+mo(i,p)]; ++i) {
					// Intentionally empty.
				}
#ifndef NDEBUG
				out << "Match got as far as i = " << i << std::endl;
#endif
				// Report any match, which is the case if we got all the way to i == p.length().
				if (i == p.length()) {
					match_count++;
				}
				
				index_range_t known_dead;
				// Based on that match attempt (and where it went wrong), what do we now know, in terms of live/dead, etc.?
				//	We have tried at position probe, so that's now dead...ie. [probe, probe] is a dead-zone.
				//	We can shift left (from probe) by shl(i) to get the next possible attempt position to the left. That means that (probe - shl(i), probe] is a dead-zone.
				//	We can shift right (from probe as well) by shr(i) to get the next possible attempt position to the right, meaning that [probe, probe + shr(i)) is also dead.
				//	Taken together, we know that (probe - shl(i), probe + shr(i)) is dead...or in our [) form: [probe - shl(i) + 1, probe + shr(i)) is dead.
				// NB: I've written shl(i) here (same for shr), but really mean sh.shift_left(i, in, probe).
				// Let's update our known-dead:
				known_dead.first = probe - sh.shift_left(i, in, probe) + 1;
				known_dead.second = probe + sh.shift_right(i, in, probe);
				// What remains for us to check? Well, we started with [low, high) and now we have the known_dead, which might split [low, high)
				//	Left of probe, we should go and check in [low, known_dead.first)
				//	Right of probe, we should go and check in [known_dead.second, high)
				//		More on this later, because we may gain some info after having gone left.
				
#ifndef NDEBUG
				// Let's output those shifts for debugging purposes.
				out << "Will now shift left/right by " << sh.shift_left(i, in, probe) << '/' << sh.shift_right(i, in, probe) << std::endl;
				out << "Known dead-zone is [" << known_dead.first << ',' << known_dead.second << ')' << std::endl;
#endif
#ifndef NDEBUG
				out << "Left will be [" << low << ',' << known_dead.first << ") and we think (for now) right will be [" << known_dead.second << ',' << high << ')' << std::endl;
#endif
				// Push the two new ranges that we need to do...make sure to do them in this order, so that the left range is dealt-with first.
				todo.push(index_range_t(known_dead.second, high));
				todo.push(index_range_t(low, known_dead.first));
			}
		}
		
#ifndef NDEBUG
		out << "Total matches = " << match_count << std::endl << std::endl;
#endif
		return match_count;
	}
};

// This (template) class has the same interface as Pattern_matcher_iterative, but does several things manually.
// For details of the template parameters, see the above.
template<typename MOrder = MO_rev, typename Shifter = Shifter_kmp_bm<MOrder>, typename Probe_chooser = Probe_chooser_mid>
class Pattern_matcher_iterative_raw {
private:
	// Ranges are a recurring theme, so us a typedef to give it a short name.
	struct index_range_t {
		int first, second;
	};
	
	// We need to keep a copy of the keyword being matched.
	// TODO: do we want to have a separate "match attempter" object?
	std::string p;
	// Keep a local copy of the match order object.
	// NB: !!! It's critical that this is declared before the Shifter, since mo is used in the initialization of the Shifter (in the constructor).
	MOrder mo;
	
	// We need a shifter object, which is built around the same match order, otherwise this won't work at all!!
	Shifter sh;
	
public:
	// The constructor takes the keyword and builds a shifter as well.
	Pattern_matcher_iterative_raw(const std::string &key) : p(key), sh(key, mo) {
		// Intentionally empty...all fields initialized in the initializer list.
	}
	
	// Do some matching and report how many matches were found.
	// TODO: see if we really want to register matches via an ostream. In reality, we'd probably want to do some array or a counter, etc.
	int match(const std::string &in, std::ostream &out) const {
		// Set our attempt count to 0, since this is a "per match" attempt counter.
		register int match_count(0);
		
		// Also keep track of the part which is known to be dead...[0,d)...initially [0,0)
		register int d(0);
		
		// Use a local stack to keep track of the ranges that still need to be examined:
		// This stack is done manually, so really bad style :-?
		static index_range_t todo[1024];
		register int tos(0);
#define STACK_EMPTY (tos == 0)
#define STACK_POP (--tos)
#define STACK_TOP (todo[tos])
#define STACK_PUSH(x,y) ++tos; STACK_TOP.first=(x); STACK_TOP.second=(y)
		
		// Load it with the initial range, leaving out the stuff at the end, since it can't contribute to a match:
		STACK_PUSH(0, in.length()-(p.length()-1));
		
		while (!STACK_EMPTY) {
			// Get the next one to do...
			// ...give them nice names...
			register int low(STACK_TOP.first), high(STACK_TOP.second);
			// ...and remove it from the stack.
			STACK_POP;
			// Adjust this element according to the already-known-dead [0,d)
			low = std::max(low, d);
			
			// If we've got an empty live-range, then just return.
			if (high <= low) {
#ifndef NDEBUG
				//			out << debug_indent << "\tEmpty range, so doing nothing and returning right away. This doesn't count as an attempt!" << std::endl << std::endl;
#endif
				d = low;
				// Intentionally empty.
			} else {
#ifndef NDEBUG
				out << std::endl << "Invoked with a live-zone of [" << low << ',' << high << ')' << std::endl;
#endif
				// We know we have at least one place to check now.
				assert(low < high);
				
				// Our aim now is to check all of [low, high), reporting matches as we go.
				// Let's start with some point in this live-zone, e.g. the middle, but let's determine it using our Probe_chooser function object.
				register int probe(Probe_chooser()(low, high));
#ifndef NDEBUG
				// Do some pretty printing here, giving a picture of the input string with the keyword aligned at the right place.
				out << "Attempting a match at " << probe << std::endl;
				out << "Picture" << std::endl;
				out << '\t' << in << std::endl << '\t';
				for (int z(0); z < probe; ++z) {
					out << ' ';
				}
				out << p << std::endl;
#endif
				
				// We'll need an index into p for matching. Declare it here since we'll need it after the for-loop here.
				int i;
				// Here, we make use of the "match-order" function "mo", which is just a permutation on [0, p.length()).
				// Keep going until we have a mismatch:
				for (i = 0; i < p.length() && p[mo(i,p)] == in[probe+mo(i,p)]; ++i) {
					// Intentionally empty.
				}
#ifndef NDEBUG
				out << "Match got as far as i = " << i << std::endl;
#endif
				// Report any match, which is the case if we got all the way to i == p.length().
				if (i == p.length()) {
					match_count++;
				}
				
				register index_range_t known_dead;
				// Based on that match attempt (and where it went wrong), what do we now know, in terms of live/dead, etc.?
				//	We have tried at position probe, so that's now dead...ie. [probe, probe] is a dead-zone.
				//	We can shift left (from probe) by shl(i) to get the next possible attempt position to the left. That means that (probe - shl(i), probe] is a dead-zone.
				//	We can shift right (from probe as well) by shr(i) to get the next possible attempt position to the right, meaning that [probe, probe + shr(i)) is also dead.
				//	Taken together, we know that (probe - shl(i), probe + shr(i)) is dead...or in our [) form: [probe - shl(i) + 1, probe + shr(i)) is dead.
				// NB: I've written shl(i) here (same for shr), but really mean sh.shift_left(i, in, probe).
				// Let's update our known-dead:
				known_dead.first = probe - sh.shift_left(i, in, probe) + 1;
				known_dead.second = probe + sh.shift_right(i, in, probe);
				// What remains for us to check? Well, we started with [low, high) and now we have the known_dead, which might split [low, high)
				//	Left of probe, we should go and check in [low, known_dead.first)
				//	Right of probe, we should go and check in [known_dead.second, high)
				//		More on this later, because we may gain some info after having gone left.
				
#ifndef NDEBUG
				// Let's output those shifts for debugging purposes.
				out << "Will now shift left/right by " << sh.shift_left(i, in, probe) << '/' << sh.shift_right(i, in, probe) << std::endl;
				out << "Known dead-zone is [" << known_dead.first << ',' << known_dead.second << ')' << std::endl;
#endif
#ifndef NDEBUG
				out << "Left will be [" << low << ',' << known_dead.first << ") and we think (for now) right will be [" << known_dead.second << ',' << high << ')' << std::endl;
#endif
				// Push the two new ranges that we need to do...make sure to do them in this order, so that the left range is dealt-with first.
				STACK_PUSH(known_dead.second, high);
				STACK_PUSH(low, known_dead.first);
			}
		}
#undef STACK_PUSH
#undef STACK_TOP
#undef STACK_POP
#undef STACK_EMPTY
		
#ifndef NDEBUG
		out << "Total matches = " << match_count << std::endl << std::endl;
#endif
		return match_count;
	}
};

// This (template) class has the same interface as Pattern_matcher_iterative, but does several things manually.
// For details of the template parameters, see the above.
template<typename MOrder = MO_rev, typename Shifter = Shifter_kmp_bm<MOrder>, typename Probe_chooser = Probe_chooser_mid>
class Pattern_matcher_iterative_raw_nr {
private:
	// Ranges are a recurring theme, so us a typedef to give it a short name.
	struct index_range_t {
		int first, second;
	};
	
	// We need to keep a copy of the keyword being matched.
	// TODO: do we want to have a separate "match attempter" object?
	std::string p;
	// Keep a local copy of the match order object.
	// NB: !!! It's critical that this is declared before the Shifter, since mo is used in the initialization of the Shifter (in the constructor).
	MOrder mo;
	
	// We need a shifter object, which is built around the same match order, otherwise this won't work at all!!
	Shifter sh;
	
public:
	// The constructor takes the keyword and builds a shifter as well.
	Pattern_matcher_iterative_raw_nr(const std::string &key) : p(key), sh(key, mo) {
		// Intentionally empty...all fields initialized in the initializer list.
	}
	
	// Do some matching and report how many matches were found.
	// TODO: see if we really want to register matches via an ostream. In reality, we'd probably want to do some array or a counter, etc.
	int match(const std::string &in, std::ostream &out) const {
		// Set our attempt count to 0, since this is a "per match" attempt counter.
		register int match_count(0);
		
		
		// Use a local stack to keep track of the ranges that still need to be examined:
		// This stack is done manually, so really bad style :-?
		static index_range_t todo[1024];
		register int tos(0);
#define STACK_EMPTY (tos == 0)
#define STACK_POP (--tos)
#define STACK_TOP (todo[tos])
#define STACK_PUSH(x,y) ++tos; STACK_TOP.first=(x); STACK_TOP.second=(y)
		
		// Load it with the initial range, leaving out the stuff at the end, since it can't contribute to a match:
		if (in.length() < p.length()) {
			return match_count;
		}
		STACK_PUSH(0, in.length()-(p.length()-1));
		
		while (!STACK_EMPTY) {
			// Get the next one to do...
			// ...give them nice names...
			register int low(STACK_TOP.first), high(STACK_TOP.second);
			// ...and remove it from the stack.
			STACK_POP;
			
#ifndef NDEBUG
			out << std::endl << "Invoked with a live-zone of [" << low << ',' << high << ')' << std::endl;
#endif
			// We know we have at least one place to check now.
			assert(low < high);
			
			// Our aim now is to check all of [low, high), reporting matches as we go.
			// Let's start with some point in this live-zone, e.g. the middle, but let's determine it using our Probe_chooser function object.
			register int probe(Probe_chooser()(low, high));
#ifndef NDEBUG
			// Do some pretty printing here, giving a picture of the input string with the keyword aligned at the right place.
			out << "Attempting a match at " << probe << std::endl;
			out << "Picture" << std::endl;
			out << '\t' << in << std::endl << '\t';
			for (int z(0); z < probe; ++z) {
				out << ' ';
			}
			out << p << std::endl;
#endif
			
			// We'll need an index into p for matching. Declare it here since we'll need it after the for-loop here.
			int i;
			// Here, we make use of the "match-order" function "mo", which is just a permutation on [0, p.length()).
			// Keep going until we have a mismatch:
			for (i = 0; i < p.length() && p[mo(i,p)] == in[probe+mo(i,p)]; ++i) {
				// Intentionally empty.
			}
#ifndef NDEBUG
			out << "Match got as far as i = " << i << std::endl;
#endif
			// Report any match, which is the case if we got all the way to i == p.length().
			if (i == p.length()) {
				match_count++;
			}
			
			register index_range_t known_dead;
			// Based on that match attempt (and where it went wrong), what do we now know, in terms of live/dead, etc.?
			//	We have tried at position probe, so that's now dead...ie. [probe, probe] is a dead-zone.
			//	We can shift left (from probe) by shl(i) to get the next possible attempt position to the left. That means that (probe - shl(i), probe] is a dead-zone.
			//	We can shift right (from probe as well) by shr(i) to get the next possible attempt position to the right, meaning that [probe, probe + shr(i)) is also dead.
			//	Taken together, we know that (probe - shl(i), probe + shr(i)) is dead...or in our [) form: [probe - shl(i) + 1, probe + shr(i)) is dead.
			// NB: I've written shl(i) here (same for shr), but really mean sh.shift_left(i, in, probe).
			// Let's update our known-dead:
			known_dead.first = probe - sh.shift_left(i, in, probe) + 1;
			known_dead.second = probe + sh.shift_right(i, in, probe);
			// What remains for us to check? Well, we started with [low, high) and now we have the known_dead, which might split [low, high)
			//	Left of probe, we should go and check in [low, known_dead.first)
			//	Right of probe, we should go and check in [known_dead.second, high)
			//		More on this later, because we may gain some info after having gone left.
			
#ifndef NDEBUG
			// Let's output those shifts for debugging purposes.
			out << "Will now shift left/right by " << sh.shift_left(i, in, probe) << '/' << sh.shift_right(i, in, probe) << std::endl;
			out << "Known dead-zone is [" << known_dead.first << ',' << known_dead.second << ')' << std::endl;
#endif
#ifndef NDEBUG
			out << "Left will be [" << low << ',' << known_dead.first << ") and we think (for now) right will be [" << known_dead.second << ',' << high << ')' << std::endl;
#endif
			// Push the two new ranges that we need to do...make sure to do them in this order, so that the left range is dealt-with first.
			if (known_dead.second < high) {
				STACK_PUSH(known_dead.second, high);
			}
			if (low < known_dead.first) {
				STACK_PUSH(low, known_dead.first);
			}
		}
#undef STACK_PUSH
#undef STACK_TOP
#undef STACK_POP
#undef STACK_EMPTY
		
#ifndef NDEBUG
		out << "Total matches = " << match_count << std::endl << std::endl;
#endif
		return match_count;
	}
};

/******************
 MAIN
 ******************/


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
                
                /* BAD!! Runs out of stack:
                 // Bruce Horspool
                 //std::cout << "Horspool matching:" << std::endl;
                 Pattern_matcher<MO_rev, Shifter_end_chars, Probe_chooser_LR> m_horspool(keyword);
                 start = mach_absolute_time();
                 m_horspool.match(input, std::cout);
                 top = mach_absolute_time();
                 duration = stop - start;
                 
                 // Convert to nanoseconds
                 duration *= info.numer;
                 duration /= info.denom;
                 
                 std::cout << duration << "," ;*/
                
                // Deadzone
                // Do forward matching, the best shifter and choosing mid-point probing (default)?
                //std::cout << "(Forward, end-point character, mid-point probe) matching:" << std::endl;
                Pattern_matcher_iterative_raw<MO_fwd, Shifter_kmp_bm<MO_fwd> > m_bw4(keyword);
                start = mach_absolute_time();
                m_bw4.match(input, std::cout);
                stop = mach_absolute_time();
                duration = stop - start;
                
                // Convert to nanoseconds
                duration *= info.numer;
                duration /= info.denom;
                
                std::cout << duration << "\t";
                
                
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
