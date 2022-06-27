# mspbwt

══ Failed ══════════════════════════════════════════════════════════════════════
── 1. Failure (test-unit-encode.R:352:9): can encode and decode a complete encod
usg[K + 1, ] not equivalent to a[, 2].
Lengths differ: 12 is not 11



Current
=============
 - OK so now have a bug, probably a big one
 - I think the updating is wrong for ms-pbwt
 - Better understand current (normal) PBWT, and think about what should happen given arbitrary (not just 0-1) break / swap

 - Then
 - Now have the whole building in C++
 - So could test on an e.g. 2 Mbp region, see if it can work
 - Then maybe some mild profiling
 - Then broadly ready for use!




 

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
