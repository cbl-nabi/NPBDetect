import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def convert_to_fasta(in_genbank, out_fasta):
    SeqIO.convert(in_genbank, "genbank", out_fasta, "fasta")
        
def extract_proteins(gbk, out_fasta):
    handle = open(out_fasta,"w")
    prots = gbk.features
    cont = 0
    for P in prots:
        if(P.type == 'CDS'):
            if(P.qualifiers.get('protein_id') != None ):
                pname = P.qualifiers['protein_id'][0]
            else:
                pname = 'P' + str(cont)
                cont += 1
            mseq = SeqRecord(Seq(P.qualifiers['translation'][0]), id = pname, name= '',
                description="")
            SeqIO.write(mseq, handle, "fasta")