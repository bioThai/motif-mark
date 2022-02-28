# Motif Mark OOCA: Thai Nguyen

List of classes to include (explained in detail below):
- Gene
- Motif
- MotifMarker


## Gene class

Each Gene object holds data about:
  - its name or a unique gene ID value to store
  - the start and stop position of that gene (so that intron/exon locations can be labelled correctly)
  - its nucleotide sequence? (is this even needed? or would this take up too much memory to be practical?)
  - and intron/exon locations (splice site 5' and 3' locations, so store each set of splice site data in a tuple?)
      - hold this data in a dictionary 
          - key: strings with value of "introns" or "exons", 
          - value: list of tuples, with each tuple containing start and end positions of each intron or exon)

      - but what if the locations are out of order, how would you sort???


## Motif class

Each Motif object holds data about:
  - its name?
  - which binding protein it is associated with
  - its sequence 
  - which gene it belongs to
  - and its beginning/end locations on a gene (what if it's found on multiple genes, then what?)


## MotifMarker class (drives the whole motif marker program)

Contains argparse method, FASTA and text file reading method, and methods to determine how to find motifs within each gene, and draws their positions on each gene. Instantiates Gene and Motif objects, saves them in data structures, and accesses them to draw Pycairo images.

Algorithm:
- First, read in list of motifs from a txt file, 
  - for each line/motif in file, instantiates a Motif object with data about its sequence, name, etc
  - stores each motif object in a big list of motifs?

- read each inputed FASTA file
- is each FASTA file all the seqs for a single gene? If so, need a method to instantiate a Gene object with data about its name, sequence, etc.
  - for each sequence in gene, 
      - read through each nucleotide (incrementing a position counter for each NT)
          - have an overall general variable to track NT position currently being read
          - if NT is uppercase or lowercase, then set an uppercase flag to true or false, respectively. Also save this intron/exon start position in a variable.
          - then keep reading next NTs to figure out at which position the intron/exon ends. 
              - Position where NT changes from uppercase to lowercase (or vice versa) will be saved as start position of new intron/exon
                   - end position of the just-read intron/exon will be (this position) minus 1
  
  - iterate through list of all Motif objects. For each motif object, search for it in each gene in gene list.
    - use regex to get all start and end positions of each motif
    - need a dictionary (call it gene_motifs_dict):
      -key: gene objs
      -value: list of motif objects and locations of each motif object along each gene
          - nested list, with each inner list holding motif objects in element 0 and tuples of start/stop positions of each motif in subsequent elements. [[motif1, (start1, stop1), (start2, stop2), etc], [motif2, (start1, stop1), (start2, stop2), etc]]

  - call an generate_img() function to generate a pycairo image with all genes and all locations of motifs along each gene (this info was saved in gene_motifs_dict)
    - basically loop through each item in dict and plot it

  - to be continued...
