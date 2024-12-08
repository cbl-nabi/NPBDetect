#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 30 14:13:30 2024

@author: dsing243
"""

import sys
import os
import argparse, textwrap

import numpy
import pandas as pd
import tempfile

from Bio import SeqIO
import torch

from bin.gbk_to_fa import convert_to_fasta, extract_proteins
from bin.kmerprofiles import nucmer, protmer
from bin.PFAM_feats import pfam_domains
from bin.model import get_model

__VERSION__ = '1.0.0'



def make_predictions(gbk, prediction_type, output_dir, verbose):
    if( os.path.isfile(gbk) == False):
        print('Incorrect path! check the path again')
        sys.exit(2)
    
    if( gbk.split('.')[-1] != 'gbk'):
        print('Unknown format detected!!! File is not GBK')
        sys.exit(2)
    
    gbk = os.path.abspath(gbk)
    gbk_name = gbk.split('/')[-1].split('.')[0]
    base_dir = os.getcwd()
    if(base_dir[-1] != '/'):
        base_dir += '/'
    temp_dir = tempfile.TemporaryDirectory()
    activity_classes = ['antibacterial', 'antifungal', 'cytotoxic_antitumor', 'inhibitor', 'surfactant', 'antiprotozoal', 'antiviral', 'siderophore']
    model = get_model(base_dir, verbose)
    if(model == None):
        sys.exit(2)
    rec = SeqIO.read(open(gbk, 'r'),"genbank")
    if(verbose > 0):
        print('Extracting features from GBK file')
    
    if(verbose > 1):
        print('Extracting fasta sequence')
    nucl_file = temp_dir.name + '/' + gbk_name + '.fa'
    if(verbose > 1):
        print('Writing Nucleotide file:', nucl_file)
    convert_to_fasta(gbk, nucl_file)
    
    
    prot_file = temp_dir.name + '/p' + gbk_name + '.fa'
    if(verbose > 1):
        print('Writing Protein file:', nucl_file)
    extract_proteins(rec, prot_file)
    
    if(verbose > 1):
        print('Extracting PFAM domains')
    PFAMD = pfam_domains()
    pfam_matrix = PFAMD.get_PFAM_domains(rec, base_dir)
    
    if(verbose > 1):
        print('Extracting nucleotide features')
    NUCFEAT = nucmer(nucl_file, base_dir, activity_classes)
    nuc_matrix = NUCFEAT.get_profiles()

    if(verbose > 1):
        print('Extracting protein features')
    PROTFEAT = protmer(prot_file, base_dir, activity_classes)
    prot_matrix = PROTFEAT.get_profiles()

    input_matrix = numpy.concatenate((pfam_matrix, nuc_matrix, prot_matrix), axis = 1)
    logits = model(torch.Tensor(input_matrix).to(dtype=torch.float32))
    sigm = 1/(1 + numpy.exp(-logits.detach().numpy()))
    predictions = pd.DataFrame(sigm).T
    
    predictions.columns = ['Probabilities']
    predictions.index = activity_classes
    predictions['Predictions'] = numpy.int_(sigm > 0.5).T
    if(verbose > 0):
        print('Predicting bioactivities')
        
    if(prediction_type == 'ORG'):
        if(verbose > 1):
            print('Original')
        if(output_dir != None):
            if(verbose > 1):
                print('Writing outputs in csv file')
            predictions.to_csv(output_dir + '/' + gbk_name + '.csv')
        else:
            print(predictions)
    elif(prediction_type == 'HC'):
        if(verbose > 1):
            print('High confidence clases')
        high_conf = predictions.loc[['antibacterial', 'antifungal', 'cytotoxic_antitumor', 'siderophore']]
        if(output_dir != None):
            if(verbose > 1):
                print('Writing outputs in csv file')
            high_conf.to_csv(output_dir + '/' + gbk_name + '.csv')
        else:
            print(high_conf)
    temp_dir.cleanup()


def print_help( argv ):
    parser  =  argparse.ArgumentParser(
                        prog = 'NPBDetect',
                        add_help = False,
                        description = 'A Neural network model to detect bioactivities',
                        epilog = 'By Dalwinder Singh and Hemant Goyat')
    subparsers = parser.add_subparsers(dest = "commands")
    parser._subparsers.title  =  "commands"
    
    progress_parser = subparsers.add_parser('predict',
                                     add_help = False,
                                     formatter_class = argparse.RawTextHelpFormatter,
                                     help = '''Predict bioactivities''',
                                     epilog='Example: NPBDetect predict --gbk <local_path_to_gbk> --pred HC --out_dir <local_path_to_dir>')
    group_usage = progress_parser.add_argument_group('Usage [options]',
                                                     None)
    
    group_usage.add_argument('--gbk',
                             dest = 'gbk',
                             action = 'store',
                             type = str, 
                             help = '''Path to GBK file for prediction''')
    
    group_usage.add_argument('--pred',
                             dest = 'ptype',
                             action = 'store',
                             choices = ['ORG', 'HC'], 
                             default = 'HC',
                             type = str,
                             help = '''Output predictions for original or high \nconfidence classes only. (Default: HC) \nORG: Predictions for 8 classes.\nHC: Predictions for top 4 classes only.
                                       ''')
    
    group_usage.add_argument('--out_dir',
                             dest = 'output_dir',
                             action = 'store',
                             type = str, 
                             default = None, 
                             help = '''Output directory to write predictions in csv format\n(Default: print to console)''')
    
    group_usage.add_argument('-v','--verbose',
                         dest = 'v',
                         action = 'store',
                         type = int, 
                         default = 0,
                         help = textwrap.dedent('''Controls the verbosity: \nhigher means more messages (Default: 0)\n=0 Silent mode: print predictions\n>0 Main steps: Model Loading, feature extraction and\n   prediction\n>1 More info. about feature extraction and prediction'''))
                                    
    group_usage.add_argument('-h','--help',
                        action='help',
                        # action='store_true',
                        dest='show_help_predict',
                        help = "Print help")
    
    
    parser.add_argument('-h','--help',
                        action = "help", 
                        help = "Basic help")
    
    parser.add_argument('-v','--version', 
                        action = 'version', 
                        version = __VERSION__, 
                        help = 'Print current version')
    
    if( len(argv) == 1):
        parser.parse_args(['-h'])
        sys.exit(1)
    
    # args  =  parser.parse_args('predict -h '.split())
    args  =  parser.parse_args(argv[1:])
    if(args.gbk != None):
        make_predictions(args.gbk, args.ptype, args.output_dir, args.v)
    
if __name__ == "__main__":
    # print_help('')
    try:
        print_help(sys.argv)
    except:
        sys.exit(2)