# mspbwt

todo:
I have the first part working, in that the multi-symbol version is the same as the two symbol version, when run on two symbols
It also works on a larger SNP set, when there are a few symbols, which is great

The next steps are to make it more efficient, figure out a better storage mechanism
Then if that makes sense, port it to C++ for some more tests!

todo:
Start writing components in C++, starting with the encoding
 - DONE - complicated encoding
 - DONE - simple encoding
 - DONE - decoding
 - TODO - can't take lists out of Rcpp without copying! so convert the list of lists into 3 matrices
   - DONE - write R code to take list of lists and make single matrices
   - TODO - check this works
   - DONE - in R, make and check minimal decoding work with offsets
   - TODO - in C++, make and check minimal decoding work with offsets   
   - TODO - do the same thing for maximal decodings
   - what is the most efficient way to build the single matrices on the fly? can I use a combo of Rcpp and R to more efficiently build them, passing by reference etc. do I need to, is it fast enough in R already
 - TODO - putting those both together
 - TODO - puttingt the whole building the indices into C++ (confirm I can avoid type conversion stuff with lists?)


todo:
I can now make a and d efficiently
I now have the alternative version for "us" now which is usge_all
Only unclear bit in building is whether the usg matrix gets too big / unweildly? Not sure. Probably worth profiling etc at some point
Next is to look at the R code a bit more and simplify
Next is to turn into C++
Works pretty well it looks like!
