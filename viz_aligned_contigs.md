# Visualizing aligned contigs and alignments

Yuying Rong
20250407

After aligning split contigs all-vs-all, I wanted to view how contigs aligned to each other and where the aligned regions are within the contigs. 

Therefore, I carved out a random contig from the (filtered) alignment file. 

```sh
# bash

# first, filter the alignment file
# by removing alignments <100 kb
# this reduces the alignment file size from 59 G to 7.6 M,
# faster for the demonstration
cat $path/outs/Altus.split.self.paf | awk -F "  " '$11 >= 100000' > $path/outs/Altus.split.self_100000.paf && echo "$(date +"%T"): awk: done"

# then, extract all alignments to the chosen query contig
# each row in the .paf is an alignment; column 1 contains the query contigs
# let's choose contig "Altus_ch02_h0_290:1-5199970" to visualize
cat Altus.split.self_100000.paf | grep "^Altus_ch02_h0_290:1-5199970" | cut -f 1-12 > Altus_ch02_h0_290_1-5199970.paf
```

Now I have a file containing all alignments to contig "Altus_ch02_h0_290:1-5199970" as the query contig.


### Preparations
From this file, I want to get the start & end coordinates of alignments and of contigs. 

Then I will calculate the coordinates in relation to the chosen query contig. 

```py
# python

# this script reads in the alignment file and
# generates an array capturing the coordinates of
# all contigs and their alignments to the query contig

def array_contigs_alignments(rows):
    
    # arrays of contigs and alignments as pairs of start and end coords
    contigs = []
    alignments = []
    
    # array of flattened intervals; add the start coord of target contig
    array = [1,]

    # read file by row and add coord info
    for row in rows:
        row = row.split()
        direction = str(row[4])
        qstart = int(row[2])
        qend = int(row[3])
        qlen = int(row[1])
        tstart = int(row[7])
        tend = int(row[8])
        tlen = int(row[6])
        
        alignment = (qstart, qend)# this is not a true alignment, but just to add the query in for viz
        alignments.append(alignment)

        array.append(alignment[0])# also define start & end interval borders for the query contig
        array.append(alignment[1])
        
        # add sense direction target contig and alignment intervals
        if direction == '+':
            contig = ((qstart - tstart), (qstart - tstart + tlen))
            contigs.append(contig)
            array.append(contig[0])
            array.append(contig[1])
            
        # add antisense direction target contig and alignment intervals
        # note: first nucleotide in negative coord of the query contig (extending left side) starts at 0, not -1
        elif direction == '-':
            contig = ((qstart - (tlen - tend)), (tend + qstart))
            contigs.append(contig)
            array.append(contig[0])
            array.append(contig[1])
            
            
    # also append the end coord of query contig to interval array
    array.append(qlen)
    
    # also append the interval of query contig to list of contig tuples
    contigs.append((1, qlen))
    
    # also append the interval of alignment on the query contig
    # don't want to visualize the alignments in the query contig, so (0, 0)
    alignments.append((0, 0))
    
    array = sorted(list(set(array)))
    
    return array, contigs, alignments


f = open('/Users/shell/MSc/WUR/paf/Altus_ch02_h0_290_1-5199970_100000.paf', 'r')

array, contigs, alignments = array_contigs_alignments(f)

f.close()


print(f'array: {array}')# contains all interval start & end coordinates, flattened
print(f'array: {contigs}')# contains all contig start & end coordinates
print(f'alignments: {alignments}')# contains all alignment start & end coordinates

# sanity check
print(f'len(contigs) == len(alignments)? {len(contigs) == len(alignments)}')
# should be True
# there should be as many entries of contigs as entries of alignments
```

Output:
```
array: [-3103213, -2897023, -734672, -49999, -30776, -7521, -772, -75, -52, 0, 1, 14, 71, 30776, 52510, 83286, 183723, 357784, 372329, 448457, 475140, 559083, 848296, 922447, 939245, 941282, 1273230, 1273976, 1286612, 1291148, 1298017, 1316061, 1316259, 1372597, 1394221, 1437517, 1610459, 1949395, 2007211, 2098379, 2175218, 2256331, 2376769, 2376774, 2586143, 2601730, 2612659, 2613321, 2618417, 2655953, 2655981, 2677255, 2707176, 2737098, 2842063, 3017290, 3092126, 3108544, 3116017, 3130485, 3130503, 3148128, 3148194, 3167506, 3306177, 3306197, 3405824, 3410590, 3416000, 3437978, 3437995, 3516497, 3683232, 3694672, 3701892, 3702552, 3728854, 3848427, 4067005, 4111552, 4124135, 4128716, 4224181, 4229652, 4266911, 4266928, 4436489, 4489519, 4625249, 4747820, 5149971, 5169194, 5192449, 5199198, 5199895, 5199918, 5199963, 5199970, 5200041, 5230746, 5238787, 5280574, 6193606, 6565514, 8816886, 13820768]
array: [(0, 5199970), (4266911, 5280574), (4266911, 5238787), (1316061, 2677255), (2612659, 3092126), (0, 848296), (2655953, 3167506), (939245, 1298017), (939245, 1291148), (3148128, 3405824), (3437978, 3694672), (3437978, 3683232), (0, 475140), (3130485, 3306197), (-30776, 5169194), (30776, 5230746), (2376769, 2618417), (-49999, 5149971), (-7521, 5192449), (3728854, 6565514), (3017290, 6193606), (1949395, 4436489), (2737098, 13820768), (-734672, 2098379), (-2897023, 2256331), (-52, 5199918), (-75, 5199895), (-772, 5199198), (-3103213, 8816886), (71, 5200041), (1, 5199970)]
alignments: [(14, 5199963), (4266928, 5199963), (4266928, 5199963), (1316259, 2586143), (2613321, 3108544), (14, 922447), (2655981, 3116017), (941282, 1286612), (941282, 1286612), (3148194, 3410590), (3437995, 3701892), (3437995, 3702552), (14, 372329), (3130503, 3306177), (52510, 357784), (83286, 183723), (2376774, 2601730), (4111552, 4747820), (2007211, 2175218), (4124135, 4224181), (4489519, 4625249), (3416000, 3516497), (3416000, 3516497), (1273976, 1372597), (448457, 559083), (1437517, 1610459), (3848427, 4067005), (2707176, 2842063), (4128716, 4229652), (1273230, 1394221), (0, 0)]len(contigs) == len(alignments)? True
len(contigs) == len(alignments)? True
```


### Visualizing contigs and alignments

Now, I want to visualize the contigs and alignments (in a crude way).

I will draw a contig as a stretch of dots and an alignment as a stretch of squares.

Each character (dot or square) is a "pixel".

Since I cannot plot the entire length of tens of thousands of bps of alignments or contigs, I will shrink the pixel lengths proportionally by bin_size. 
```py
# python

def view_contigs_alignments(contigs, alignments, bin_size):
    # first coordinate
    min_start = (min(contigs))[0]
    print(f'first coordinate: {min_start}')
    print(f'number of intervals (including target contig): {len(contigs)}')
    for contig, alignment in zip(contigs, alignments):
        print(f'contig: {contig}; alignment: {alignment}')###
        print((' ' * ((contig[0] - min_start)//bin_size) +
               '•' * ((alignment[0] - contig[0])//bin_size) +
               '■' * ((alignment[1] - alignment[0])//bin_size) +
               #' ' * ((alignment[1] - alignment[0])//bin_size) +# this shows purged alignments as empty stretches
               '•' * ((contig[1] - alignment[1])//bin_size)))# the last row is the example contig
    return 0

# I manually set a bin_size of 200 kb
view_contigs_alignments(contigs, alignments, 200000)
```

Output:
```
first coordinate: -3103213
number of intervals (including target contig): 31
contig: (0, 5199970); alignment: (14, 5199963)
               ■■■■■■■■■■■■■■■■■■■■■■■■■
contig: (4266911, 5280574); alignment: (4266928, 5199963)
                                    ■■■■
contig: (4266911, 5238787); alignment: (4266928, 5199963)
                                    ■■■■
contig: (1316061, 2677255); alignment: (1316259, 2586143)
                      ■■■■■■
contig: (2612659, 3092126); alignment: (2613321, 3108544)
                            ■■
contig: (0, 848296); alignment: (14, 922447)
               ■■■■
contig: (2655953, 3167506); alignment: (2655981, 3116017)
                            ■■
contig: (939245, 1298017); alignment: (941282, 1286612)
                    ■
contig: (939245, 1291148); alignment: (941282, 1286612)
                    ■
contig: (3148128, 3405824); alignment: (3148194, 3410590)
                               ■
contig: (3437978, 3694672); alignment: (3437995, 3701892)
                                ■
contig: (3437978, 3683232); alignment: (3437995, 3702552)
                                ■
contig: (0, 475140); alignment: (14, 372329)
               ■
contig: (3130485, 3306197); alignment: (3130503, 3306177)
                               
contig: (-30776, 5169194); alignment: (52510, 357784)
               ■••••••••••••••••••••••••
contig: (30776, 5230746); alignment: (83286, 183723)
               •••••••••••••••••••••••••
contig: (2376769, 2618417); alignment: (2376774, 2601730)
                           ■
contig: (-49999, 5149971); alignment: (4111552, 4747820)
               ••••••••••••••••••••■■■••
contig: (-7521, 5192449); alignment: (2007211, 2175218)
               •••••••••••••••••••••••••
contig: (3728854, 6565514); alignment: (4124135, 4224181)
                                  ••••••••••••
contig: (3017290, 6193606); alignment: (4489519, 4625249)
                              ••••••••••••••
contig: (1949395, 4436489); alignment: (3416000, 3516497)
                         •••••••••••
contig: (2737098, 13820768); alignment: (3416000, 3516497)
                             ••••••••••••••••••••••••••••••••••••••••••••••••••••••
contig: (-734672, 2098379); alignment: (1273976, 1372597)
           •••••••••••••
contig: (-2897023, 2256331); alignment: (448457, 559083)
 ••••••••••••••••••••••••
contig: (-52, 5199918); alignment: (1437517, 1610459)
               ••••••••••••••••••••••••
contig: (-75, 5199895); alignment: (3848427, 4067005)
               •••••••••••••••••••■•••••
contig: (-772, 5199198); alignment: (2707176, 2842063)
               ••••••••••••••••••••••••
contig: (-3103213, 8816886); alignment: (4128716, 4229652)
••••••••••••••••••••••••••••••••••••••••••••••••••••••••••
contig: (71, 5200041); alignment: (1273230, 1394221)
               •••••••••••••••••••••••••
contig: (1, 5199970); alignment: (0, 0)
               •••••••••••••••••••••••••
```

### Depth of coverage by contigs

Now I am also interested in knowing how many contigs overlap with each other per-interval:

This is the per-interval (or per-base) coverage by contigs.

```py
# python

# each interval is a contig

def interval_coverage(array, intervals):
    # array of coverage of each interval
    coverage = [0] * (len(array) - 1)
    
    
    # use interval info to update coverage
    for interval in intervals:
        for i in range(len(array)):
            
            # start adding 1 once new_interval[0] == interval[0]
            if interval[0] == array[i]:
                
                # while new_interval[1] is larger than interval[1], keep adding 1
                while interval[1] > array[i]:
                    try:
                        #print(('array[i], interval[0], interval[1]', array[i], interval[0], interval[1]))
                        coverage[i] += 1
                        i += 1
                    
                    except IndexError:
                        # when i += 1 shoots out of the last item in coverage
                        break
    
    return coverage


coverage = interval_coverage(array, contigs)
print(f'coverage: {coverage}')

# coverage is calculated for the interval between each coord in array, therefore:
print(f'len(array)-1 == len(coverage)? {len(array)-1 == len(coverage)}')
# should be True
```


Output:
```py
coverage: [3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
len(array)-1 == len(coverage)? True
```

Finally, to plot the numerical per-base depth of coverage by contigs (shrinked proportionally by bin_size):
```py
# python

def coverage_histogram(array, coverage, bin_size):
    max_coverage = 0
    for i in range(len(array)-1):
        max_coverage = max(coverage[i], max_coverage)
        
        # end = "" ensures printing without newline between loops
        print(f"{'■' * coverage[i]}\n" * ((array[i+1] - array[i] + 1)//bin_size), end = "")
        #print(f"{'■' * coverage[i]}\n" * ((array[i+1] - array[i] + 1)//bin_size))

    print(f"max_coverage: {max_coverage}")
    
    return 0

# bin_size = 200 kb
coverage_histogram(array, coverage, 200000)
```


Output:
```
■
■■
■■
■■
■■
■■
■■
■■
■■
■■
■■
■■■
■■■
■■■
■■■■■■■■■■■■■■
■■■■■■■■■■■■■■■
■■■■■■■■■■■■■■
■■■■■■■■■■■■■■
■■■■■■■■■■■■■■■
■■■■■■■■■■■■■■■■
■■■■■■■■■■■■■■■■
■■■■
■■■■
■■■■
■■■■
■■■
■■
■■
■■
■■
■■
■■
■■
■■
■■
■■
■■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
■
max_coverage: 17
```

Here I presented my scripts to quickly generate a simple visualization of the contigs and alignments.

Based on the visualizations for a few more query contigs, I got an idea on how contigs align to each other in the complex all-vs-all situation. 
