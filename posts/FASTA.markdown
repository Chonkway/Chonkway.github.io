
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
I decided 512MB file sizes would be a reasonable split size, and while I could force users to use SRA Toolkit or other programs to split the `.fasta` files into those limits, it seemed nicer to just do it for them and let them toss in a full sizes file no matter the size. In order to do this, we will use the method described here https://biopython.org/wiki/Split_large_file

```python
def batch_iterator(iterator, batch_size):...


record_iter = list(SeqIO.parse(open("DRR128265.fasta"), "fasta"))
for i, batch in enumerate(batch_iterator(record_iter, 10000)):
    filename = "group_%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
```
I found that if you don't do `list(SeqIO.parse)` you run into an issue of python not being able to apply `.next()`. If you copy paste from the site it won't work, though maybe that was more of an issue of me tossing stuff together.