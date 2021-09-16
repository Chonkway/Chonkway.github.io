
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
FASTA files do not strictly come as `.fasta`, so I will use the `os` module to account for the common extension aliases.


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

* This still feels slow and I have a feeling that there are libraries that will do everything I want rendering my written "library" and this jumbled loop useless. But it's a learning experience. 


------------
***UPDATE*** 

I was right. I basically rendered all the above work useless lmao.

This does it all, and it doesn't even need my AASequence code either.
 
```python
print("--------------")
filename = input("Enter your filename (including extension). Ensure it is in the root directory.")
finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
q = input("Is your sequence an mRNA sequence that needs translated? (y/n)") #Accounts for different sequence types (mRNA, DNA)
for entry in SeqIO.parse(filename, "fasta"):
    seqstring = Seq(entry.seq) #Grabs the sequence for each entry and turns it into a Seq object
    if q.lower() == "y": 
        seqstring = Seq(entry.seq).translate() #Translates mRNA to relevant protein sequence
    else:
        pass #Placeholder for other conversion types
    analyzed_seq = ProteinAnalysis(str(seqstring))
    aminoacids = analyzed_seq.count_amino_acids()
    for key in finalseqcount:
        if key in analyzed_seq.amino_acids_content:
            finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
```

This part sets a dictionary of all 20 amino acids with a value of 0. We use `SeqIO` to parse the file and create a `Seq Object` of the sequence string. We add an option to use the `.translate()` function in Biopython just incase the fasta file the user is inputting is mRNA. We convert the sequence into a string so we can use `.count_amino_acids()` and `amino_acids_content()`.

FInally, we update the amino acid value in the dictionary to get a final count of all the different acids present in the sequence.


```python
with open('AminoAcids.json') as json_file:
    data = json.load(json_file)

CCount = 0
NCount = 0
PCount = 0

for i in finalseqcount.keys():
    for key in data:
        if i == key in data:
            CCount = CCount + (finalseqcount[i]*data[key][0])
            NCount = NCount + (finalseqcount[i]*data[key][2])
            PCount = PCount + finalseqcount[i]
```
This part gives us the Carbon, Nitrogen and Phosphorous count of the sequence. By using the `json` module, we open the `AminoAcids.json` file that contains information on the structure of each amino acid. We check to make sure the key in both dictonaries are the same and then apply the relevant operations to update our final counts.


## Update (9-16-21)

After some playing around, reading and time I have a program I'm happy enough with (as far as the bones go. I'm unsure of the accuracy and the speed needs reworking some day)


```python
#Reimpliement batch iteration -

def batch_iterator(iterator, batch_size):
    """
    Returns multiple lists of length batch_size.
    """

    entry = True #Ensures the loop runs once
    while entry:
        batch = []
        while len(batch) < batch_size: 
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break#End of file
            batch.append(entry)
        if batch:
            yield batch

print("--------------")
filename = input("Enter your filename(including extension). Ensure it is in the root directory.") #Sets target file
query = input("If the file is significantly large, it is recommended that you split it. Would you like to split the file?(y/n)")
seqtype = input("Does your file need transcribed?(y/n)") #Used for a check below to use SeqIO's .translate() module

if query.lower() == "y": #Checks if you want/need to split the file
    record_iter = SeqIO.parse(open(filename), "fasta") #Uses SeqIO to parse the file as the iterator (subject to change to a different fasta parser)
    for i, batch in enumerate(batch_iterator(record_iter, 60000)):
        file = "group_%i.fasta" % (i + 1)
        with open(file, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, file))
```

I started by re-implementing the option to use the batch_iterator. This is only an option for those who want to do larger files on it, though I don't recomment it at this point in time. The batch size of `60000` was for a very, very large file and it took several hours.

There is some flavor text for readability, and I have it handle whether the user needs their file translated or not all at the start.

```python
#Scan directory for common fasta file extensions using OS module -
ext = ('.fasta', '.fna', '.fnn', '.faa', '.frn' , '.fa')

PCount = 0 # Phosphorous count will be sum of sequence length pre-translation, creates variable.
finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

#Loop over all files, count their Amino Acids and create a total - 
print("Now parsing files. Depending on the number of files and their size, this may take a while...")
if query == 'y':
    for files in os.listdir(): #Scan for all files (applicable only if files are split)
        if files.endswith(ext):
            if files != filename:
                fileparse = SeqIO.parse(files, "fasta") #Begins the file parse. Trying to keep this as the only parse step to avoid costly time
else:
    fileparse = SeqIO.parse(filename, "fasta") #If there's only one file used, no need to scan directory
```

This portion is mostly the same, though now we only scan the directory if the user decided to say `yes` to splitting the file. Otherwise, it loads the entire file into memory if they decided it was sufficiently small.

```python
for entry in fileparse: #Depending on if the file is large/a large series of a lot of smaller files, this can take a long time
    PCount = PCount + len(str(entry.seq))

    if seqtype.lower() == "y": # Check for mRNA
        mRNA_translate = Seq(str(entry.seq)).translate()
        analyzed_seq = ProteinAnalysis(str(mRNA_translate)) #Allows .count_amino_acids() to apply to the files

    else:
        analyzed_seq = ProteinAnalysis(str(entry.seq))

    aacount = analyzed_seq.count_amino_acids()

    for key in finalseqcount:
        if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
            finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]

    else:
        pass

print("Finished.")
print("\n")

print("Your final amino acid counts are: \n")
print(finalseqcount)
log = input("Would you like to write these to a log file?(y/n)")

with open('AminoAcids.json') as json_file:
    data = json.load(json_file)

NCount = 0
for i in finalseqcount.keys():
    for key in data:
        if i == key in data:
            NCount = NCount + (finalseqcount[i]*data[key][2])
            
if log.lower() == "y":   
    with open("results_'{}'.txt".format(filename), "x") as results: #Logs results because I'm sick of waiting for console
        results.write("Amino acid count:\n")
        results.write(str(finalseqcount))
        results.write("\n")
        results.write("Phosphorous count: \n")
        results.write(str(PCount))
        results.write("\n")
        results.write("Nitrogen count: \n")
        results.write(str(NCount))
```

Finally, the script generates a dictionary of the acid counts and then allows the user to write a `results.txt` logfile for their run if they want. It takes the name of `results_{filename}` for ease of keeping track. 