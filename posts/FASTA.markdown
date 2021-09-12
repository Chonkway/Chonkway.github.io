
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

Now that we have a more managable ~300 files, we need to isolate the part we care about: the sequence.

Since we have all these smaller `.fasta` files, we can more reasonably commit them to memory. My plan is to use SeqIO to parse each file, strip the sequence and append it to a list in a similar fashion to what we just did. I'll set a batch_size again to prevent the list from getting too large, write the list to a text file, and repeat until we have multiple text files. We will finally pass these text files through my program, sum each instances Nitrogen, Phosphorous and Carbon count and create a final log file of the results.

```python
ext = ('.fasta', '.fna', '.fnn', '.faa', '.frn' , '.fa')
```
FASAT files do not strictly come as `.fasta`, so I will use the `os` module to account for the common extension aliases.


~~The idea is that we scan the directory for any common FASTA file extensions. The headers of these files all start with that `>` character:
> '>DRR128265.1 1 length=301

~~So we look for any line that contains this character, ignore it and then create a newfile without it. This will give us files without any information we don't want to pass through our tool. We add a check to ensure our original fasta file is untouched by the header stripping.~~



Originally, the code looked like this:

```python
for files in os.listdir():
    if files.endswith(ext) and files != filename:
        with open(files, "r+") as oldfile, open(files, "r+") as cleanfile:
            lines = oldfile.readlines()
            for line in lines:
                if '>' in line: 
                    pass
                else:
                    cleanfile.write(line)
                    print("Succesfully removed FASTA header from" + "" + str(oldfile))

    else:
        pass
```

And this DOES work, but over the course of a huge file like an entire organism it would take way, way too long.

Instead, after a while of digging I found that the `SeqIO.parse()` returns a `SeqRecord` object and I can just take that object and call the sequence itself. Since I figure SeqIO is much more optimized than my loop, it should work much faster. This could certainly be used to remove the need for the second chunk of code but as of now I can't be bothered.

```python
for files in os.listdir():
    if files.endswith(ext):  #Ensures original FASTA is safe?
        sequences = []
        ffile = open(files)
        for seq_record in SeqIO.parse(ffile, "fasta"):
            for n in seq_record:
                s = seq_record.seq
                sequences.append(SeqRecord(s, id="", description=""))
            SeqIO.write(sequences, "testfile.fasta", "fasta")
        print("Stripped header from " + str(files))
    else:
        pass
```
As of now, the main problems are as follows:

* I realized that there is nothing preventing the program from reading output files if they're fasta files, meaning it will indefinitely loop over the new files with their header stripped. This will result in duplicate strings and, more pressingly, an infinitely running program.

* This still feels slow and I have a feeling that there are libraries that will do everything I want rendering my written "library" and this jumbled loop useless. But it's a learning experience
