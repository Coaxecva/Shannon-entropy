import sys, math
from itertools import groupby, islice
from collections import Counter

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

def sliding_window(seq, n=2):
    """Returns a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   
       this function courtesy http://stackoverflow.com/questions/6822725/rolling-or-sliding-window-iterator-in-python
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def Shannon_entropy(seq, k):
    c = Counter()    
    for kmer in sliding_window(seq, k):        
        new_count = c[kmer] + 1
        c[kmer] = new_count
    s = -0
    for i in c:
        print i, c[i]
        p = float(c[i]) / (len(seq) - k + 1)
        s += p * math.log(p, 2)
    return (-s)

if __name__ == '__main__':
    fasta, k = sys.argv[1], sys.argv[2]

    for v in fasta_iter(fasta):
        print Shannon_entropy(v[1], int(k))

