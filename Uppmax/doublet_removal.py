# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 10:25:37 2022

@author: matil
"""
import os
import scipy.io
import doubletdetection
from scipy import sparse
import numpy as np
import argparse
import gzip
import csv

#how to run: python doublet_removal -i (input directory where matrix file is located) -o (wanted output directory)

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--indir", help="Input direcotry where matrix file is found")
parser.add_argument("-o", "--outdir", help="Output directory to place doublet free matrix in")

args = parser.parse_args()

matrix_dir = args.indir
# read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx.gz"))

# read barcode file
barcodes_path = os.path.join(matrix_dir, "barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

#predict doublets
clf = doubletdetection.BoostClassifier(n_iters=50)
doublets_mat = clf.fit(np.transpose(mat)).predict(p_thresh=1e-7, voter_thresh=0.8)

#remove doublets
matrix_no_dbs =np.transpose(mat.toarray())[ ~doublets_mat.astype(bool)]

#write outfile
scipy.io.mmwrite(os.path.join(args.outdir,"matrix_no_doublets.mtx"),np.transpose(sparse.csr_matrix(matrix_no_dbs).tocoo()))

#write barcodes with doublets removed
barcodes_out_path = os.path.join(args.outdir, "barcodes_no_doublets.tsv")
with open(barcodes_out_path, 'w') as csvfile:
    # creating a csv writer object
    csvwriter = csv.writer(csvfile, delimiter ='\n')
    # writing the data rows
    csvwriter.writerows([np.array(barcodes)[ ~doublets_mat.astype(bool)]])

#plot convergence
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=1e-7, voter_thresh=0.8)
