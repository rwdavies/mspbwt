# mspbwt

todo:
I have the first part working, in that the multi-symbol version is the same as the two symbol version, when run on two symbols
It also works on a larger SNP set, when there are a few symbols, which is great

The next steps are to make it more efficient, figure out a better storage mechanism
Then if that makes sense, port it to C++ for some more tests!

todo:
Start writing components in C++, starting with the encoding
 - DONE - complicated encoding
 - TODO - simple encoding
 - TODO - decoding
 - TODO - putting those both together
 - TODO - puttingt the whole building the indices into C++ (confirm I can avoid type conversion stuff with lists?)


todo:
I can now make a and d efficiently
I now have the alternative version for "us" now which is usge_all
Only unclear bit in building is whether the usg matrix gets too big / unweildly? Not sure. Probably worth profiling etc at some point
Next is to look at the R code a bit more and simplify
Next is to turn into C++
Works pretty well it looks like!
