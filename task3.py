from collections import deque
import gzip

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

in_file = gzip.open('chr21.fa.gz', 'rt')  # Open input file
in_file.readline()  # Discard title line
out_file = open('paired_bins.tsv', 'w')  # Create output file
out_file.write(f'bin1_start\tbin1_end\tbin2_start\tbin2_end\n')  # Print header
high_GC_dict = {}  # Create dictionary to hold relevant bins
high_GC_list = (deque(), deque())  # Keep track of bins to discard
bins_in_pairs_set = set()  # Create a set to hold all bins in pairs
for i in range(46709):  # There are 46709 full 1kb bins
    gc = 0  # Counter for GC bases
    for j in range(1000):  # Iterate over each 1kb bin base-by-base
        base = in_file.read(1)
        if base == 'G' or base == 'C' or base == 'S':
            gc += 1  # Increment if the base is strong
    if gc > 500:  # If GC content is >50%, the bin qualifies for further consideration
        for k in range(gc - 19, gc + 20):  # Get lists of all bins with |ΔGC%| < 2%
            for m in high_GC_dict.get(k, deque()):  # Print coordinates of all bins that form a valid pair with this bin
                out_file.write(f'{m * 1000}\t{(m + 1) * 1000}\t{i * 1000}\t{(i + 1) * 1000}\n')
                bins_in_pairs_set.add(m)  # Add both bins to the bin set
                bins_in_pairs_set.add(i)
        high_GC_dict.setdefault(gc, deque()).append(i)  # Add this bin to the deque keyed to the right gc
        high_GC_list[0].append(i)  # Add this bin and the gc pair to the deques of all qualified bins
        high_GC_list[1].append(gc)
    if len(high_GC_list[0]) > 0 and high_GC_list[0][0] == (i - 999):  # Remove bin >1Mb away from all deques
        high_GC_dict[high_GC_list[1][0]].popleft()
        high_GC_list[0].popleft()
        high_GC_list[1].popleft()
gc = 0  # Repeat for final incomplete bin
for j in range(983):
    base = in_file.read(1)
    if base == 'G' or base == 'C' or base == 'S':
        gc += 1
gc *= 1000/983  # Adjust gc count to account for lower length
if gc > 500:
    for k in range(gc - 19, gc + 20):
        for m in high_GC_dict.get(k, deque()):
            out_file.write(f'{m * 1000}\t{(m + 1) * 1000}\t{46709000}\t{46709983}\n')
            bins_in_pairs_set.add(m)
            bins_in_pairs_set.add(46709)
in_file.close()
out_file.close()
print(len(bins_in_pairs_set))

# Average Time Complexity: O(n) - Unlikely to ever near worst case of O(n^2)
# Average Space Complexity: O(n) - Due to need to track of all bins that are printed else it would be O(1)
