
----
I am writing this prior to this project for record keeping and to track my thoughts.

Using the AASequence tool on my github along with the biopython library I need to take input from a `.fasta` file, in this case `DRR128265.fasta`

----
First, we need to parse the file. In order to do this we use the SeqIO function in Biopython.

```python
from Bio import SeqIO 


for seq_record in SeqIO.parse("DRR128265.fasta", "fasta"):
    print(seq_record)
```
Unfortunately, for many full sized `.fasta` files there are tens of millions of lines to run through and holding that much information in memory just simply isn't realistic. The most reasonable approach would be to split one large file into several, smaller files for handling. For our file `DRR128265.fasta` we have a staggering ~92 million lines broken up into these sequence chunks. 
```
>DRR128265.1 1 length=301
CTTTCATCGTTTACTCTGGCAAAATGAATTGTAACCCTCCATTCTCCTCCTTTAACTTAT
TTCTACTCAACTTTTTCTTCCCTTCTTCCTCCTTTTTCACATTTTTATTCTCCCCTCTCA
TCCTCTCCACTTCCTTTTTCTTATCCTACTTCTTTCTAACTCTCATCTTCTCCTTCTTTG
TTTTCTCAATTCTCTCCCCCCTCTTCAACACCTTCTTTCTTCCGTTTACCCCTTGTTTTT
CTTAAGTGCACTCTTCTTCTCTCCTTTCTCTTTTTCTTTTCTCTTTTTCCTTCTTCTTCT
T
```
While I could force users to use SRA Toolkit or other programs to split the `.fasta` files into those limits while downloading, it seemed nicer to just do it for them and let them toss in a full sizes file no matter the size. This will also force a standard naming convention for the rest of the program to follow. In order to do this, we will use the method described here https://biopython.org/wiki/Split_large_file

```python
def batch_iterator(iterator, batch_size):...


record_iter = SeqIO.parse(open("DRR128265.fasta"), "fasta")
for i, batch in enumerate(batch_iterator(record_iter, 25000)):
    filename = "group_%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
```
Copy pasting the code from the site didn't work for me, and I found after some trial and error that you can't do `iterator.next()` in the `batch_iterator` function. Instead, `next(iterator)` worked for me.

I also found that for my 4GB file, doing a batch size of 1000 created like...7000 files so I changed it to 25,000 