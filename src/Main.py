#!/usr/bin/env python3
""" Own functions """
from read_clean import *

""" Modules """
from Bio import SeqIO
from Bio.Seq import Seq
import pdb
from dendropy.calculate import treecompare
import dendropy
from dendropy import Tree
import tempfile
from subprocess import Popen, PIPE, run
import os
import re
from shutil import copyfile
import time
import csv
import glob

def fastPhylo(normal_file, tempdir):
        """ Reduced noise tree"""
        with open(tempdir+'/red_fast', 'w') as file:
            run(['fastprot','reduced_file'], stdin=PIPE, stdout=file, stderr=PIPE, encoding = 'utf8', cwd=tempdir)

        with open(tempdir+'/red_tree', 'w') as file:
            run(['fnj','red_fast', '-O', 'newick'], stdin=PIPE, stdout=file, stderr=PIPE, encoding = 'utf8', cwd=tempdir)

        """ Normal tree"""
        with open(tempdir+'/normal_fast', 'w') as file:
            run(['fastprot', normal_file], stdin=PIPE, stdout=file, stderr=PIPE, encoding = 'utf8')

        with open(tempdir+'/normal_tree', 'w') as file:
            run(['fnj','normal_fast', '-O', 'newick'], stdin=PIPE, stdout=file, stderr=PIPE, encoding = 'utf8', cwd=tempdir)

        return tree_compare (tempdir)

def tree_compare(tempdir):
        tns = dendropy.TaxonNamespace()
        tree1 = Tree.get_from_path(
                tempdir+"/ref.tree",
                "newick",
                taxon_namespace=tns)
        tree2 = Tree.get_from_path(
                tempdir+"/normal_tree",
                "newick",
                taxon_namespace=tns)
        tree3 = Tree.get_from_path(
                tempdir+"/red_tree",
                "newick",
                taxon_namespace=tns)
        tree1.encode_bipartitions()
        tree2.encode_bipartitions()
        tree3.encode_bipartitions()
        distance_normal = treecompare.symmetric_difference(tree1, tree2)
        distance_reduced = treecompare.symmetric_difference(tree1, tree3)
        return distance_normal, distance_reduced

def main():
    os.chdir('..')
    dirpath = os.path.dirname(os.path.realpath(__file__))
    directory = dirpath + '/data'
    # Choose the most recent data folder
    directory = max(glob.glob(os.path.join(directory, '*/')), key=os.path.getmtime)
    all_testdata =[]

    for data in os.listdir(directory):
        data = directory +'/' + data
        os.chdir(data)
        # Create a temporary directory where we will work
        with tempfile.TemporaryDirectory(dir = data) as tempdir:
            # Find and copy the reference tree to temporary directory
            for filename in os.listdir(data):
                if os.path.isfile(filename):
                    if filename[0] != '.':
                        if re.search(r'tree', filename):
                            copyfile(filename,tempdir+'/ref.tree')
                            statistics =[(filename,filename+'_reduced')]
                            break

            for filename in os.listdir(data):
                if os.path.isfile(filename):
                    if filename[0] != '.' and not re.search(r'tree', filename):
                            # Create the reduced file
                            if read_file(filename, tempdir):
                                # Run the subprocesses
                                statistics.append(fastPhylo(filename, tempdir))
                            else:
                                pass
        all_testdata.append(statistics)
    localtime = time.localtime(time.time())
    # Save results in results folder
    os.chdir(dirpath + '/results')
    today = str(localtime.tm_year)+'_'+str(localtime.tm_mon)+'_'+str(localtime.tm_mday)
    if not os.path.exists(today):
        os.mkdir(today)
    os.chdir(os.getcwd()+'/'+today)
    with open('results.csv','w') as f:
        for row in all_testdata:
            writer = csv.writer(f)
            writer.writerows(row)

if __name__ == '__main__':
    main()
