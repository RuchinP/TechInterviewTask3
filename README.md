# TechInterviewTask3
# Purpose: This algorithm splits chromosome 21 into 1kb non-overlapping bins and prints out the coordinates of all
# pairs of bins less than 1Mb apart that both have >50% GC and GC content within 2% of each other into an output file
# before finally printing the total number of bins that meet these criteria. The coordinates are 0-based and the 1Mb
# distance is measured between the starts of the bin intervals. The program has also been optimized for efficiency.

# Method: Everything is  hardcoded in for efficiency (even if the improvement is very little). The fasta file is read
# one character at a time (to prevent complete reduplication) in 1000 character sections. If the character is G,
# C or S (either G or C), a gc counter increments and if the gc counter is over 500 at the end of the 1kb section,
# the starting index/1000 is stored in a dictionary of deques in the deque that is keyed to the sections gc value. By
# recording it this way, there is no need to do up to 999 comparisons for each bin to find each pair and instead you
# only need to access 39 dictionary values each time. The index and gc are then stored in a pair of deques so that
# when the stored become >1Mb away from the current one, they can be removed from the dictionary without searches
# otherwise comparisons would have to be done to avoid printing invalid pairs. Before insertion into the dictionary
# and the list, the dictionary is searched for keys within 20 of the gc value and the returned deques of bins are
# iterated through and the coordinates are printed to the out file. The bin indexes are put into a set to keep track of
# how many different bins are in the output file. Deques are heavily used instead of lists since there is a lot of
# removal of data from the beginning of lists and random access is unneeded
