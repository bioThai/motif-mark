#!/usr/bin/python

import cairo    #pycairo library to draw images
import argparse
import re       #regular expressions
import seaborn  #need for color palette


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
        iupac_to_regex_dict: dict = {"W":"[ATU]", "S":"[CG]", "M":"[AC]", "K":"[GTU]", "R":"[AG]", "Y":"[CTU]", "B":"[CGTU]", "D":"[AGTU]", "H":"[ACTU]", "V":"[ACG]", "N":"[ACGTU]", "U":"[TU]", "T":"[TU]"}
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
    __slots__ = ['name', 'chrom_name', 'start_pos','stop_pos','seq', 'intron_exon_dict']

    # ----- Constructor -----
    def __init__(self, gene_name: str, chromosome_name: str, start_position: int, stop_position: int, sequence: str, intron_exon_dictionary: dict):
        '''Constructor method to initialize new Gene objects.'''
        self.name = gene_name
        self.chrom_name = chromosome_name
        self.start_pos = start_position
        self.stop_pos = stop_position
        self.seq = sequence
        self.intron_exon_dict = intron_exon_dictionary
            # {"intron": list(), "exon": list()}
            # key: "intron" or "exon"
            # value: list of tuples containing start and stop positions of each intron or exon encountered in a gene sequence
            # note that the positions stored in the tuples are lower bound inclusive and upper bound exclusive
    
    # # ----- Other methods -----
    # def find_motifs(self, regex_motif_list: list):
    #     '''Takes in a list of motif sequences that have been converted to their regular expression equivalent, and searches for each regex motif in the Gene object's sequence.'''



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


def create_image(gene_motifs_dict: dict, motif_object_list: list, output_image_name: str):
    '''Takes in an output image filename and a dictionary with Gene objects as keys and values as nested list, with each inner list holding motif objects in element 0 and tuples of start/stop positions of each motif in subsequent elements (i.e., [[motif1, (start1, stop1), (start2, stop2), etc], [motif2, (start1, stop1), (start2, stop2), etc]]). Creates a Pycairo image of all the genes, including their introns and exons, along with locations of each motif along each gene.Saves the Pycairo image into a file with the given output image filename.'''

    #local variables
    GENE_SURFACE_WIDTH: float = 800         #width of surface object for each gene drawing in points 
    GENE_SURFACE_HEIGHT: float = 150        #height of surface object for each gene drawing in points 
    LEGEND_SURFACE_WIDTH: float = 200       #width of legend surface object in points
    TITLE_SURFACE_WIDTH: float = GENE_SURFACE_WIDTH
    TITLE_SURFACE_HEIGHT: float = 50
    PARENT_SURFACE_WIDTH: float = GENE_SURFACE_WIDTH +  LEGEND_SURFACE_WIDTH     #width (in points) of parent surface object that will hold all other drawing sub-surfaces
    PARENT_SURFACE_HEIGHT: float = 0        #initialize to 0 then calculate below

    longest_gene_seq: str = ""
    longest_motif_seq: str = ""
    padding_dict: dict = {"top": 20, "right": 20, "bottom": 20, "left": 20} #number of points to pad the top, right, bottom, and left sides of a gene image element
    colors_dict: dict = {"intron": (0.5, 0.5, 0.5), "exon": (0.35, 0.35, 0.35)} 
        #holds color scheme for images. Initialized with intron and exon colors. Motif colors will be added later.
        #key: motif seq or the literal strings "exon" or "intron"
        #value: tuple holding RGB color scheme for motifs and introns/exons (R,G,B)


    #before creating image surface, get length of longest gene seq and longest motif seq
    for gene_obj in gene_motifs_dict:
        if len(gene_obj.seq) > len(longest_gene_seq):
            longest_gene_seq = gene_obj.seq
    
    for motif_obj in motif_object_list:
        if len(motif_obj.seq) > len(longest_motif_seq):
            longest_motif_seq = motif_obj.seq

    #create parent surface
    PARENT_SURFACE_HEIGHT = TITLE_SURFACE_HEIGHT + (len(gene_motifs_dict) * GENE_SURFACE_HEIGHT)
    parent_surface = cairo.ImageSurface(cairo.Format.RGB24, PARENT_SURFACE_WIDTH, PARENT_SURFACE_HEIGHT)

    #add title surface to parent surface
    title_surface = parent_surface.create_for_rectangle(0, 0, TITLE_SURFACE_WIDTH, TITLE_SURFACE_HEIGHT)
    title_context = cairo.Context(title_surface)
    title_context.select_font_face("monospace", cairo.FontSlant.NORMAL, cairo.FontWeight.BOLD)
    title_context.set_source_rgba(0.35, 0.35, 0.35, 1) #set label color to light gray
    title_context.set_font_size(20)
    title_context.move_to(padding_dict["left"], padding_dict["top"])
    title_context.show_text(output_image_name)

    #add legend surface to parent surface
    legend_surface = parent_surface.create_for_rectangle(TITLE_SURFACE_WIDTH, 0, LEGEND_SURFACE_WIDTH, PARENT_SURFACE_HEIGHT)
    legend_context = cairo.Context(legend_surface)
    legend_context.select_font_face("monospace", cairo.FontSlant.NORMAL, cairo.FontWeight.BOLD)
    legend_context.set_source_rgba(0.8, 0.8, 0.8, 1) #set label color to light gray
    legend_context.set_font_size(20)
    legend_context.move_to(padding_dict["left"], TITLE_SURFACE_HEIGHT + padding_dict["top"])
    legend_context.show_text("LEGEND:")

    #display exon and intron color-codes in legend
    legend_context.set_font_size(20)
    legend_context.set_source_rgba(colors_dict["intron"][0], colors_dict["intron"][1], colors_dict["intron"][2], 0.7)
    legend_context.move_to(padding_dict["left"], TITLE_SURFACE_HEIGHT + (padding_dict["top"]*2))
    legend_context.set_line_width(3)
    legend_context.line_to(padding_dict["left"] + 16, TITLE_SURFACE_HEIGHT + (padding_dict["top"]*2))
    legend_context.stroke()
    legend_context.move_to(padding_dict["left"] + 20, 4 + TITLE_SURFACE_HEIGHT + (padding_dict["top"]*2))
    legend_context.show_text("intron")

    legend_context.set_font_size(20)
    legend_context.set_source_rgba(colors_dict["exon"][0], colors_dict["exon"][1], colors_dict["exon"][2], 0.7)
    legend_context.move_to(padding_dict["left"], TITLE_SURFACE_HEIGHT + (padding_dict["top"]*3))
    legend_context.set_line_width(16)
    legend_context.line_to(padding_dict["left"] + 16, TITLE_SURFACE_HEIGHT + (padding_dict["top"]*3))
    legend_context.stroke()
    legend_context.move_to(padding_dict["left"] + 20, 4 + TITLE_SURFACE_HEIGHT + (padding_dict["top"]*3))
    legend_context.show_text("exon")


    #add motif sequences to legend
    motif_color_palette = seaborn.color_palette(None, len(motif_object_list))
    for motif_obj in motif_object_list:
        motif_index = motif_object_list.index(motif_obj)

        #add each motif's RGB color scheme to colors dict
        colors_dict[motif_obj.seq] = motif_color_palette[motif_index]

        #display color-coded motifs in legend
        legend_context.set_font_size(20)
        legend_context.set_source_rgba(colors_dict[motif_obj.seq][0], colors_dict[motif_obj.seq][1], colors_dict[motif_obj.seq][2], 0.7)
        legend_context.move_to(padding_dict["left"], TITLE_SURFACE_HEIGHT + (padding_dict["top"]*(motif_index + 4)))
        legend_context.set_line_width(16)
        legend_context.line_to(padding_dict["left"] + 16, TITLE_SURFACE_HEIGHT + (padding_dict["top"]*(motif_index + 4)))
        legend_context.stroke()

        legend_context.move_to(padding_dict["left"] + 20, 4 + TITLE_SURFACE_HEIGHT + (padding_dict["top"]*(motif_index + 4)))
        legend_context.show_text(motif_obj.seq)


    #add each gene surface to parent surface
    gene_scale_factor =  GENE_SURFACE_WIDTH / len(longest_gene_seq)
    gene_obj_list = list(gene_motifs_dict.keys())
    for gene_obj in gene_obj_list:

        gene_index = gene_obj_list.index(gene_obj)
        gene_surface = parent_surface.create_for_rectangle(padding_dict["left"], TITLE_SURFACE_HEIGHT + (gene_index*GENE_SURFACE_HEIGHT), GENE_SURFACE_WIDTH, GENE_SURFACE_HEIGHT)
        gene_context = cairo.Context(gene_surface)

        #draw gene label
        gene_context.set_source_rgba(0.8, 0.8, 0.8, 1) #set label color to light gray
        gene_context.select_font_face("monospace", cairo.FontSlant.NORMAL, cairo.FontWeight.BOLD)
        gene_context.set_font_size(20)
        gene_context.move_to(0, padding_dict["top"])
        gene_context.show_text(gene_obj.name + " " + gene_obj.chrom_name + ":" + str(gene_obj.start_pos) + "-" + str(gene_obj.stop_pos))

        #draw introns
        gene_context.set_source_rgba(colors_dict["intron"][0], colors_dict["intron"][1], colors_dict["intron"][2], 0.8) #set color to light gray
        for intron_region in gene_obj.intron_exon_dict["intron"]:
            gene_context.set_line_width(5)
            gene_context.move_to(gene_scale_factor * intron_region[0], GENE_SURFACE_HEIGHT / 2)
            gene_context.line_to(gene_scale_factor * intron_region[1], GENE_SURFACE_HEIGHT / 2)
            gene_context.stroke()

        #draw exons
        gene_context.set_source_rgba(colors_dict["exon"][0], colors_dict["exon"][1], colors_dict["exon"][2], 0.8) #set label color to light gray
        for exon_region in gene_obj.intron_exon_dict["exon"]:
            gene_context.set_line_width(65)
            gene_context.move_to(gene_scale_factor * exon_region[0], GENE_SURFACE_HEIGHT / 2)
            gene_context.line_to(gene_scale_factor * exon_region[1], GENE_SURFACE_HEIGHT / 2)
            gene_context.stroke()

        #draw motifs
        for motif_positions_list in gene_motifs_dict[gene_obj]:
            motif_seq = motif_positions_list[0]
            gene_context.set_source_rgba(colors_dict[motif_seq][0], colors_dict[motif_seq][1], colors_dict[motif_seq][2], 0.7) #set motif color to correct color
            gene_context.set_line_width(48)

            for position_tuple in motif_positions_list[1:]:
                gene_context.move_to(gene_scale_factor * position_tuple[0], GENE_SURFACE_HEIGHT / 2)
                gene_context.line_to(gene_scale_factor * position_tuple[1], GENE_SURFACE_HEIGHT / 2)
                gene_context.stroke()



    parent_surface.write_to_png(output_image_name)



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
    gene_sequence: str = ""
    gene_intron_exon_dict: dict = {}
        # key: "intron" or "exon"
        # value: list of tuples containing start and stop positions of each intron or exon encountered in a gene sequence
            # note that the positions stored in the tuples are lower bound inclusive and upper bound exclusive
    
    gene_obj: Gene = Gene("", "", 0, 0, "", dict())
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

        # for motif_obj in motif_obj_list:
        #     print(motif_obj.seq)
            
                

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
                    gene_obj = Gene(gene_name, gene_chrom_name, gene_start_pos, gene_stop_pos, "", gene_intron_exon_dict)
                    
                    # print(gene_obj.name, gene_obj.chrom_name, gene_obj.start_pos, gene_obj.stop_pos, gene_obj.seq, gene_obj.intron_exon_dict)
                    
                else:
                    gene_obj.seq = line
                    #find positions of regex matches for all lowercase NTs (introns)
                    for match in re.finditer(r'[a-z]+', line):
                        match_position_tuple = (match.start(), match.end())
                        gene_obj.intron_exon_dict["intron"].append(match_position_tuple)
                                            
                    #find positions of regex matches for all uppercase NTs (exons)
                    for match in re.finditer(r'[A-Z]+', line):
                        match_position_tuple = (match.start(), match.end())
                        gene_obj.intron_exon_dict["exon"].append(match_position_tuple)

                    #print(gene_obj.name, gene_obj.chrom_name, gene_obj.start_pos, gene_obj.stop_pos, gene_obj.seq, gene_obj.intron_exon_dict)

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
        create_image(gene_motifs_dict, motif_obj_list, output_img_fname)


if __name__ == "__main__":
    main()
