# Import Required Packages
from Bio import SeqIO
from Bio import Entrez
import csv
from csv import DictReader

# Function content which splits a sequence into triplets
def content(seq, frame):
    for i in range(frame - 1, len(seq) - 2, 3):
        yield seq[i:i + 3]


def codon(sasp_list, codon_box_list, hgnc,  refseq_id, kozak):
    # stores 'ENSG ID', 'RefSeq ID', 'CDS Length'
    new_row = [hgnc, refseq_id, len(sasp_list) * 3]
    seq_row = [refseq_id]

    # counts frequency of each codon in sequence and adds to new row for the output
    for codon in codon_box_list:
        new_row.append(sasp_list.count(codon))
    for codon in sasp_list:
        seq_row.append(box_dict.get(codon))

    new_row.extend(kozak)

    writer.writerow(new_row)  # adds new row to output csv file for RefSeqID read
    n_writer.writerow(seq_row)

email = 'meghanakarumuri@gmail.com'
fname = "Gestalt_dgtpn.csv"
Entrez.email = email
added_ref_seq = []

box_dict = {
    "TTT": "A", "TTC": "B", "TTA": "C", "TTG": "D", "TCT": "E", "TCC": "F", "TCA": "G", "TCG": "H",
    "TAT": "I", "TAC": "J", "TAA": "K", "TAG": "L", "TGT": "M", "TGC": "N", "TGA": "O", "TGG": "P",
    "CTT": "Q", "CTC": "R", "CTA": "S", "CTG": "T", "CCT": "U", "CCC": "V", "CCA": "W", "CCG": "X",
    "CAT": "Y", "CAC": "Z", "CAA": "a", "CAG": "b", "CGT": "c", "CGC": "d", "CGA": "e", "CGG": "f",
    "ATT": "g", "ATC": "h", "ATA": "i", "ATG": "j", "ACT": "k", "ACC": "l", "ACA": "m", "ACG": "n",
    "AAT": "o", "AAC": "p", "AAA": "q", "AAG": "r", "AGT": "s", "AGC": "t", "AGA": "u", "AGG": "v",
    "GTT": "w", "GTC": "x", "GTA": "y", "GTG": "z", "GCT": "0", "GCC": "1", "GCA": "2", "GCG": "3",
    "GAT": "4", "GAC": "5", "GAA": "6", "GAG": "7", "GGT": "8", "GGC": "9", "GGA": "@", "GGG": "$"
}

# Creates output csv file for frequencies/kozak -- CHANGE NAME OF OUTPUT FILE TO DESIRED NAME
with open('freq_' + str(fname), 'w', newline='') as f:
    writer = csv.writer(f)
    row = ['hgnc_symbol', 'RefSeq ID', 'CDS Length']  # initializes first three columns

    # uses looping to add in name of each codon as a column
    for i in box_dict.keys():
        row.append(i)
    for i in range(1, 11):
        row.append('Kozak ' + str(i))
    # writes titles of columns into output csv file
    writer.writerow(row)
    # Creates output csv file for 64 char encoding -- CHANGE NAME OF OUTPUT FILE TO DESIRED NAME
    with open('seq_' + str(fname), 'w', newline='') as nf:
        n_writer = csv.writer(nf)
        n_row = ['RefSeq ID']

        n_writer.writerow(n_row)

        with open(fname, 'r') as read_obj:
            csv_dict_reader = DictReader(read_obj)

            # Runs code to fetch sequence for each RefSeqID
            for row in csv_dict_reader:
                refseq = row['refseq_mrna']
                hgnc = row['hgnc_symbol']

                # Switch the line below for a gtf file (handle3 = name of genbank file) to directly open a gbk file
                handle3 = Entrez.efetch(db="nucleotide", id=refseq, rettype="gb", retmode="text")
                record3 = SeqIO.read(handle3, "gb")

                # Finds out where the sequence starts and ends
                cds_found = False  # Keeps track if CDS with start and end of the sequence is found
                for feature in record3.features:
                    if feature.type == "CDS":
                        cds_found = True
                        start_location = feature.location.start  # Stores the start location of the strand
                        end_location = feature.location.end  # Stores the end location of the strand
                        break

                if not cds_found:
                    # print("CDS not found for " + gene + "\n")
                    continue  # moves on to next RefSeqID
                else:
                    seq = str(record3.seq[start_location:end_location])  # obtain sequence from record
                    kozak = str(record3.seq[(start_location - 6): (start_location + 4)])
                    # print("Sequence Obtained!\n")
                    # print("Sequence: " + seq)
                    print("Kozak Sequence: " + kozak)
                SASP = list(content(seq, 1))  # SASP has a list of separated codons in the sequence
                codon(SASP, box_dict.keys(), hgnc, refseq, list(kozak))
print("Finished Codon Processing! See <seq_" + str(fname) + " freq_" + str(fname) + "> for outputs.")

