# mspbwt

todo:
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


About where I am
================
In test-unit-mspbwt.R, in mspbwt.R, change it from returning usge_all to returning the super encoding (this can be done at the end? possibly with some C++ thrown in)
Then try on the test stuff at least, with maybe some more C++, that it works
Then see if a and d can also be made in the same way?

Then make sure the new super format works

Figure out exactly how I store a and d, is it the same way?
