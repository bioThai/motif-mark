#!/usr/bin/python

from sys import flags
import cairo    #pycairo library to draw images
import argparse
import re       # regular expressions


# ---------- Class definitions ----------
class Motif:
    '''Class for Motif objects for each motif read in from a text file.'''
    # specify Motif object attributes (specifying slots reduces memory usage)
    __slots__ = ['seq']

    # ----- Constructor -----
    def __init__(self, sequence: str):
        '''Constructor method to initialize new Motif objects.'''

        self.seq = sequence

    # ----- Other methods -----
    def get_regex_motif(self) -> str:
        '''Take in a raw motif sequence as read in from a text file, and converts the characters in the sequence into their IUPAC equivalents. Returns the converted sequence as a string.'''
        #local variables
        iupac_to_regex_dict: dict = {"W":"[ATU]", "S":"[CG]", "M":"[AC]", "K":"[GTU]", "R":"[AG]", "Y":"[CTU]", "B":"[CGTU]", "D":"[AGTU]", "H":"[ACTU]", "V":"[ACG]", "N":"[ACGTU]"}
        regex_str: str = ""

        #convert raw motif sequence to uppercase for consistency, 
        # and iterate through each character to translate it into its IUPAC equivalent
        for char in self.seq.upper():
            #if char is not in the regex_dict, then assume it's a regular nucleotide (ACGTU) that doesn't need IUPAC conversion
            if char not in iupac_to_regex_dict:
                regex_str += char
            #else, convert char into IUPAC equivalent
            else:
                regex_str += iupac_to_regex_dict[char]

        return regex_str




class Gene:
    '''Class for Gene objects read in from each record in a FASTA file.'''
    # specify Gene object attributes (specifying slots reduces memory usage)
    __slots__ = ['name', 'chrom_name', 'start_pos','stop_pos','intron_exon_dict']

    # ----- Constructor -----
    def __init__(self, gene_name: str, chromosome_name: str, start_position: int, stop_position: int, intron_exon_dictionary: dict):
        '''Constructor method to initialize new Gene objects.'''
        self.name = gene_name
        self.chrom_name = chromosome_name
        self.start_pos = start_position
        self.stop_pos = stop_position
        self.intron_exon_dict = intron_exon_dictionary
            # {"intron": list(), "exon": list()}
            # key: "intron" or "exon"
            # value: list of tuples containing start and stop positions of each intron or exon encountered in a gene sequence
            # note that the positions stored in the tuples are lower bound inclusive and upper bound exclusive
    
    # # ----- Other methods -----
    # def find_motifs(self, regex_motif_list: list):
    #     '''Takes in a list of motif sequences that have been converted to their regular expression equivalent, and searches for each regex motif in the Gene object's sequence.'''




# # function to draw a sample pycairo image for OOCA
# def draw_sample_img():
#     # create an svg surface object (like a canvas) on which to draw things
#     surface1 = cairo.SVGSurface("ooca_image.svg", 100, 100)

#     # create a context object and instantiate it with the surface that it will use
#     context1 = cairo.Context(surface1)

#     # draw a line (not at the origin position) on the image surface
#     context1.set_line_width(3)
#     context1.set_source_rgb(0, 0, 0.7) # set RGB colors of line
#     context1.move_to(50, 50)
#     context1.line_to(50, 80)
#     context1.stroke()

#     # draw a rectangle (not at the origin position) on the image surface
#     context1.set_line_width(1)
#     context1.set_source_rgba(0.4, 0, 0.5, 0.5) # set RGB colors and alpha level of rectangle
#     context1.rectangle(25, 35, 50, 30) #create rectangle with top left position at (25, 35), width of 50, height of 30
#     context1.fill() # fill the rectangle in with the same colors as border of rectangle
#     context1.stroke() # needed if you are drawing the outline of the rectangle but have no fill set

#     #write svg surface to png (need to do this before finishing the surface)
#     surface1.write_to_png("thai_ooca_image.png")

#     # when done drawing shapes, finish the surface so it can't be edited anymore
#     surface1.finish()


# ---------- Global function definitions ----------
def get_args():
    '''Defines/sets possible command line arguments for script'''
    parser = argparse.ArgumentParser("A program to parse text file(s) of motif sequences and input FASTA file(s) and outputs graphical images of each motif's positions along each FASTA sequence.")
    parser.add_argument("-f", "--file", nargs="+", help="Specifies input FASTA file(s). FASTA files must be unzipped and must be properly formatted.", type=str, required=True)
    parser.add_argument("-m", "--motif", nargs="+", help="Specifies input text file(s) containing motif sequences. File(s) must contain only one motif on each line in the file.", type=str, required=True)
    #--help argument is included by default, so is not specified here
    return parser.parse_args()


# def read_file(filename: str, filetype: str):
#     '''Takes in a file name and file type and parses the file accordingly.'''
#     pass


def oneline_fasta(filename: str) -> str:
    '''Reads through a FASTA file and concatenates each multi-line sequence into a single sequence line for each record. Outputs the new FASTA file in the current working directory. Returns the name of the new version of the FASTA file so it can be referenced later in the script without having to manually input it.'''
    #local variables
    output_fname: str = "oneline_" + filename.split("/")[-1]
    curr_header: str = ""
    curr_seq: str = ""      #holds current sequence line being read
    record_num: int = 0
    header_lines_list: list = []
    seq_lines_list: list = []

    with open(filename, "r") as input_fh, open(output_fname, "w") as output_fh:
        for line in input_fh:
            line = line.strip()
            
            if line.startswith(">"):    #if header line is being read
                record_num += 1
                curr_seq = ""           #for each header line read, reset the curr_seq back to "" for the next seq read
                curr_header = line
                header_lines_list.append(curr_header)
            else:
                curr_seq += line

            if curr_seq != "" and len(seq_lines_list) == record_num - 1:
                seq_lines_list.append(curr_seq)
            elif curr_seq != "" and len(seq_lines_list) == record_num:
                seq_lines_list[record_num - 1] = curr_seq    

        #write FASTA headers and concatenated sequences to output file
        for i in range(len(header_lines_list)):
            output_fh.write(header_lines_list[i] + "\n" + seq_lines_list[i] + "\n")

    return output_fname


def create_image(gene_motifs_dict: dict, output_image_name: str):
    '''Takes in an output image filename and a dictionary with Gene objects as keys and values as nested list, with each inner list holding motif objects in element 0 and tuples of start/stop positions of each motif in subsequent elements (i.e., [[motif1, (start1, stop1), (start2, stop2), etc], [motif2, (start1, stop1), (start2, stop2), etc]]). Creates a Pycairo image of all the genes, including their introns and exons, along with locations of each motif along each gene.Saves the Pycairo image into a file with the given output image filename.'''

    
    


def main():
    '''Drives order of execution for program.'''    
    #local variables
    args = get_args()
    fasta_file_list: list = args.file
    motif_file_list: list = args.motif

    motif_obj_list: list = []
    output_img_fname: str = ""
    line_tokens: list = []
    gene_name: str = ""
    gene_chrom_name: str = ""
    gene_start_pos: int = 0
    gene_stop_pos: int = 0
    gene_intron_exon_dict: dict = {}
        # key: "intron" or "exon"
        # value: list of tuples containing start and stop positions of each intron or exon encountered in a gene sequence
            # note that the positions stored in the tuples are lower bound inclusive and upper bound exclusive
    
    gene_obj: Gene = Gene("", "", 0, 0, dict())
    gene_motifs_dict: dict = {}
        # key: gene objects
        # value: nested list, with each inner list holding motif objects in element 0 and tuples of start/stop positions of each motif in subsequent elements. [[motif1, (start1, stop1), (start2, stop2), etc], [motif2, (start1, stop1), (start2, stop2), etc]]
    regex_motif: str = ""
    motif_positions_list: list = []
        # list holding motif objects in element 0 and tuples of start/stop positions of each motif in subsequent elements.
    

   
    # read motifs from file into Motif objects, and put the Motif objs in a list
    #turn this into a get_motifs_from_file(filename) function
    for input_motif_file in motif_file_list:
        with open(input_motif_file, "r") as motif_fh:
            for line in motif_fh:
                line = line.strip()
                motif_obj_list.append(Motif(line))

        for motif_obj in motif_obj_list:
            print(motif_obj.seq)
            
                

    for fasta_filename in fasta_file_list:
        #for every file, assign new empty data structures to variables holding data structures to avoid issues with unwanted/residual data being saved to them
        gene_motifs_dict = {}
        motif_positions_list = []

        #convert input fasta into file with one line per sequence instead of multiple lines per sequence, then open the new FASTA file for parsing
        with open(oneline_fasta(fasta_filename), "r") as fasta_fh:
            for line in fasta_fh:
                line = line.strip()

                if line.startswith(">"): #if header line is being read
                    line_tokens = line.split()
                    gene_name = line_tokens[0].strip(">")
                    gene_chrom_name = line_tokens[1].split(":")[0]
                    gene_start_pos = int(line_tokens[1].split(":")[1].split("-")[0])
                    gene_stop_pos = int(line_tokens[1].split(":")[1].split("-")[1])
                    gene_intron_exon_dict = {"intron": list(), "exon": list()}
                    gene_obj = Gene(gene_name, gene_chrom_name, gene_start_pos, gene_stop_pos, gene_intron_exon_dict)
                    
                    # print(gene_obj.name, gene_obj.chrom_name, gene_obj.start_pos, gene_obj.stop_pos, gene_obj.intron_exon_dict)
                    
                else:
                    #find positions of regex matches for all lowercase NTs (introns)
                    for match in re.finditer(r'[a-z]+', line):
                        match_position_tuple = (match.start(), match.end())
                        gene_obj.intron_exon_dict["intron"].append(match_position_tuple)
                                            
                    #find positions of regex matches for all uppercase NTs (exons)
                    for match in re.finditer(r'[A-Z]+', line):
                        match_position_tuple = (match.start(), match.end())
                        gene_obj.intron_exon_dict["exon"].append(match_position_tuple)

                    # print(gene_obj.name, gene_obj.chrom_name, gene_obj.start_pos, gene_obj.stop_pos, gene_obj.intron_exon_dict)

                    #find positions of motifs in the gene sequence being read
                    gene_motifs_dict[gene_obj] = []
                    for motif in motif_obj_list:
                        regex_motif = motif.get_regex_motif()
                        motif_positions_list = [motif.seq]

                        # use capture group lookahead (?=...) with regex f-string to search for overlapping motifs
                        for match in re.finditer(f'(?={regex_motif})', line, flags=re.IGNORECASE):
                            match_position_tuple = (match.start(), match.start()+len(motif.seq))
                            motif_positions_list.append(match_position_tuple)

                        gene_motifs_dict[gene_obj].append(motif_positions_list)

            # for gene_obj in gene_motifs_dict:
            #     print(gene_obj.name)
            #     print(gene_motifs_dict[gene_obj])
                            

        # generate pycairo image of motifs and their locations along all genes in file
        output_img_fname = fasta_filename.split(".fa")[0] + ".png"
        create_image(gene_motifs_dict, output_img_fname)



if __name__ == "__main__":
    main()
