from Bio import SeqIO
import numpy
import itertools
import re
import pandas as pd
import os

class nucmer:
    def __init__(self, fasta_loc, base_direc, activity_classes):
        self.fasta_loc= fasta_loc
        self.data_direc = base_direc + 'profiles/nucmer/'
        self.activity_classes = activity_classes
        self.word_size = [6, 8]
    
    def load_into_dict(self, activity, word_size):
        hexdata = dict() 
        score = pd.read_csv(self.data_direc + activity + '_' + str(word_size) + '.csv')
        for i in range(len(score)):
            t = score.iloc[i]
            hexdata.update({t['hexammer']:[t[activity],t['not_' + activity]]})
        hexprof = dict()
        for i in itertools.product(['A','C','G','T'],repeat=word_size):
            n = ''.join(i)
            hexprof.update({n:0})
        return hexdata,hexprof
    
    def compute_hexammer(self, activity, word_size):
        all_seqs = list(SeqIO.parse( self.fasta_loc, 'fasta'))[0]
        nseq = str(all_seqs.seq).upper()
        nseq = re.sub(r'[^ACGT]', '', nseq)
        hexdata,hexprof = self.load_into_dict(activity, word_size)
        step_size = 1
        for x in range(0,len(nseq) - word_size, step_size):
            sub_seq = nseq[x:x+word_size]
            hexprof[sub_seq] += 1
        sum_lr = 0
        tot_windows = 0
        for x in hexprof:
            if(hexprof[x] != 0):
                vals = hexdata[x]
                if(vals[0] == 0 and vals[1] == 0):
                    continue
                if(vals[0] > 0 and vals[1] == 0):
                    sum_lr +=1 
                if(vals[0] == 0 and vals[1] > 0):
                    sum_lr -=1
                if(vals[0] > 0 and vals[1] > 0):
                    sum_lr  +=  numpy.log(vals[0] / vals[1])
                tot_windows += 1
        if(sum_lr != 0):
            nfeat = sum_lr/tot_windows
        else:
            nfeat = 0
        return nfeat
      
    def get_profiles(self):
        nfeat = numpy.zeros((1, 16))
        idx = 0
        for wd in self.word_size:
            for ac in self.activity_classes:
                nfeat[0, idx] = self.compute_hexammer(ac, wd)
                idx += 1
        return nfeat

class protmer:
    def __init__(self, prot_loc, base_direc, activity_classes):
      self.prot_loc= prot_loc
      self.data_direc = base_direc + 'profiles/protmer/'
      self.activity_classes = activity_classes
      self.word_size = 2

    def load_into_dict(self, activity):
        hexdata = dict() 
        score = pd.read_csv(self.data_direc + activity + '_' + str(self.word_size) + '.csv', na_filter = False)
        for i in range(len(score)):
            t = score.iloc[i]
            hexdata.update({t['hexammer']:[t[activity],t['not_' + activity]]})
        hexprof = dict()
        for i in itertools.product(['A','R','N','D','B','C','E','Q','Z','G','H','I','L','K','M','F','P','S','T','W','Y','V'],repeat=self.word_size):
            n = ''.join(i)
            hexprof.update({n:0})
        return hexdata,hexprof
    
    def compute_hexammer(self, activity):
        nfeat = numpy.zeros((1,50))
        
        all_seqs = list(SeqIO.parse( self.prot_loc, 'fasta'))
        for nidx, X in enumerate(all_seqs):
            if(nidx > 49):
                break
            hexdata,hexprof = self.load_into_dict(activity)
            step_size = 1
            seq = str(X).upper()
            nseq = re.sub(r'[^ARNDBCEQZGHILKMFPSTWYV]', '', seq)
            for x in range(0,len(nseq) - self.word_size, step_size):
                sub_seq = nseq[x : x + self.word_size]
                hexprof[sub_seq] += 1
            sum_lr = 0
            tot_windows = 0
            for x in hexprof:
                if(hexprof[x] != 0):
                    vals = hexdata[x]
                    if(vals[0] == 0 and vals[1] == 0):
                        continue
                    if(vals[0] > 0 and vals[1] == 0):
                        sum_lr +=1 
                    if(vals[0] == 0 and vals[1] > 0):
                        sum_lr -=1
                    if(vals[0] > 0 and vals[1] > 0):
                        sum_lr  +=  vals[0] / vals[1]
                    tot_windows += 1
            if(sum_lr != 0):
                nfeat[0, nidx] += sum_lr/tot_windows
            else:
                nfeat[0, nidx] += 0
        return nfeat
    
    def get_profiles(self):        
        nfeat = numpy.zeros((1, 400))
        idx = 0
        for ac in self.activity_classes:  
            nfeat[0, idx:idx+50] = self.compute_hexammer(ac)
            idx += 50
        return nfeat