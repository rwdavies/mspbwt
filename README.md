# mspbwt

todo:
I have the first part working, in that the multi-symbol version is the same as the two symbol version, when run on two symbols
It also works on a larger SNP set, when there are a few symbols, which is great

The next steps are to make it more efficient, figure out a better storage mechanism
Then if that makes sense, port it to C++ for some more tests!

First, fix bug, then keep going

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
I can now make a and d efficiently
I now have the alternative version for "us" now which is usge_all
Only unclear bit in building is whether the usg matrix gets too big / unweildly? Not sure. Probably worth profiling etc at some point
Next is to look at the R code a bit more and simplify
Next is to turn into C++
Works pretty well it looks like!
