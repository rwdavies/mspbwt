# mspbwt

todo:
I have the first part working, in that the multi-symbol version is the same as the two symbol version, when run on two symbols
It also works on a larger SNP set, when there are a few symbols, which is great

The next steps are to make it more efficient, figure out a better storage mechanism
Then if that makes sense, port it to C++ for some more tests!


todo:
I can now make a and d efficiently
I need to figure out an alterative version for "us" now
First step is to abstract out us a bit
Then maybe test some ideas
I think the idea of store as
1) if rare store all instances of where it occurs
2) if common store every th value (e.g. every 1000th), then between, store runs of ups and downs
I think this should be highly efficient
Also ideally this would work natively on rhb_t because that has ALL symbols - can I do this? Need to index symbols somehow? Shouldn't be too hard I imagine?