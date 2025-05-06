import os
from Bio import SeqIO
import numpy

class pfam_domains:
    def make_matrix(self, domain_list, domain_counts):
        ins = len(domain_list)
        feats = len(domain_counts)
        data = numpy.zeros((ins,feats))
        flist = list(domain_counts.keys())
        for ipx, eachd in enumerate(domain_list):
            for i in eachd:
                if (i in flist):
                    idx = flist.index(i)
                    val = eachd[i]
                    data[ipx, idx] = val
        return data
    
    
    def get_PFAM_domains(self, gbk, base_dir, score_thresh = 40):
        counts = {}
        with open(base_dir + 'profiles/PFAM/40_4.txt', 'r') as fp:
            for I in fp.readlines():
                counts.update({I[:-1]:0})
        feat = gbk.features
        dom_dic = {}
        domain_list = []
        for F in feat:
            if(F.type == 'PFAM_domain'):
                score = float(F.qualifiers["score"][0])
                if score < score_thresh:
                    continue
                domain = F.qualifiers["description"][0]
                if(domain not in dom_dic):
                    dom_dic.update({domain:1})
                else:
                    dom_dic[domain] += 1
                if(domain not in counts):
                    continue
                else:
                    counts[domain] += 1
        domain_list.append(dom_dic)
        pfam_matrix = self.make_matrix(domain_list, counts)
        return pfam_matrix