# mspbwt

Current
=============
 - Can handle 0 in hapMatcher or Z
 - Need to write this into the C++ bit of building the index
 - At some point, need to write matching into C++ too
 - but first test a bit!
 - Then maybe some mild profiling etc
 - Then broadly ready for use!



## Notes on some variables / operations

Suppose we have large set of haplotypes with K rows and nSNPs columns, stored as binary values (called rhi_t)
We can convert that to a binary representation with K rows and nGrids columns, with nGrids being ceiling(nSNPs / 32) (called rhb_t)
And we can re-map that so that each column has values that start at 1, and increment by 1, such that the number of remapped i's is >= i + 1 (e.g. there are more 1's than 2's, and more 2's than 3's, etc) (called hapMatcher)
Now for reasons not covered here, values not remapped are set to a 0
Now we use this hapMatcher as input to build the mspbwt indices, along with a helper list called all_symbols with a description of what symbols are found in each grid, and how many of them
In building the indices, all 0's are re-mapped to the highest value observed in a given grid (e.g. if in the first grid we see (normal) values of 1,2,3,4, and there are some 0's, those are stored in hapMatcher as 0's, but used in building the indices as 5's)
For matching, for some new Z in binary representation, this can be mapped to the new hapMatcher symbol base, in a straightforward way, except for two cases
In the case that Z can't be mapped to the non-zero values, but there are 0 values, it is set to 0
In the case that Z can't be mapped to the non-zero values, and there are no 0 values, a closest match is found
As such, the re-mapped Z can't have symbols not found in hapMatcher at any grid point, and the algorithm should work without issue

## Old scratch

 



Old todo:
I have the first part working, in that the multi-symbol version is the same as the two symbol version, when run on two symbols
It also works on a larger SNP set, when there are a few symbols, which is great

Now, for the more efficient version, that is now available as an encoding (and is working)
i.e. there is a storage_info matrix, then an out_matrix, then out vector
Want to store the ms-PBWT indices using this and try it out!

The next steps are to make it more efficient, figure out a better storage mechanism
Then if that makes sense, port it to C++ for some more tests!

High level goals for TODO, with partial information
 - DONE - regular encoding
   - written in c++ and R, works, tests
 - super encoding, on its own
   - DONE written in R and tested, works
   - TODO, add in C++ for decoding using offsets and test, for both individual functions, and whole thing. make sure pass by value being used effectively here!
 - mspbwt
   - integrate newer super encoding into tests (in R)
   - write C++ version of building indices, performing super encoding, and then the matching
   - test them!




todo:
I can now make a and d efficiently (how? same approach?)
I now have the alternative version for "us" now which is usge_all
Only unclear bit in building is whether the usg matrix gets too big / unweildly? Not sure. Probably worth profiling etc at some point
Next is to look at the R code a bit more and simplify
Next is to turn into C++
Works pretty well it looks like!


TODO
================
 - make encoding of d an option
 - build the full thing then take a look
 - To do, get even larger, full sized version
 - a I am generally screwed, I think, unless I can truly figure out building on the fly etc
 - but it's not so bad?


 - Ensure functions to access u, etc, can be easily accessed and swapped out
 - Does it make sense to do a and d in the same way, or somehow take their difference? Should I have done that for u? (or problem - having to go too far?)
So I need to encode u, a, d 
Note that a, d contain a lot of identical-ness?

In test-unit-mspbwt.R, in mspbwt.R, change it from returning usge_all to returning the super encoding (this can be done at the end? possibly with some C++ thrown in)
Then try on the test stuff at least, with maybe some more C++, that it works
Then see if a and d can also be made in the same way?

Then make sure the new super format works

Figure out exactly how I store a and d, is it the same way?
