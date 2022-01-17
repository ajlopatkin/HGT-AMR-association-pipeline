#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 2 15:29:23 2021

@author: clairejohnson
"""

from reportlab.lib import colors
import argparse
import os
import pandas as pd
import collections
from reportlab.lib import colors
from Bio import SeqIO, SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics.GenomeDiagram import FeatureSet as GDFeatureSet, Diagram as GDDiagram, Track as GDTrack
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.graphics.shapes import Rect
from reportlab.pdfgen.canvas import Canvas
import random

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data_dir', dest = 'd', type=str, default='/Users/clairejohnson/Documents/pipeline_test')
parser.add_argument('-f', '--fasta_dir', dest = 'f', type=str, default='/Users/clairejohnson/Documents/plasmid_project/assemblies')
parser.add_argument('-v', '--visuals', dest = 'v', type=int, default=1)
parser.add_argument('-t', '--threshold', dest='t', type=int, default=2)
args = parser.parse_args()
parent_dir = args.d + "/raw_output"
input_dir = args.f
threshold=args.t
resist_map = {}
args_resist = {}

def abricate_parser(in_file, db):
    df = pd.read_csv(in_file, '\t')
    df["#FILE"] = df["#FILE"].str.split('/').str[-1]
    df = df[["#FILE","SEQUENCE", "GENE", "DATABASE", "START", "END", "STRAND", "%IDENTITY", "%COVERAGE", "RESISTANCE"]]
    df.columns = ["file_id","contig_id", "name", "src", "start", "end", "strand", "identity", "coverage", "resistance"]
    
    df['e_val'] = ''
    df['mobility'] = ''
    df['inc_rep'] = ''
    if db != 'plasmidfinder':
        df['type'] = 'ARG'
    else:
        df['inc_rep'] = df['name'].str.split('_').str[0]
        df['type'] = 'plasmid'
    
    df = df.reindex(['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'resistance', 'e_val', 'identity', 'coverage'], axis = 1)
    return df

def mefinder_parser(in_dir):
    df = pd.DataFrame(columns = ['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'resistance', 'e_val', 'identity', 'coverage'])
    for file in os.listdir(in_dir):
        if '.csv' in file:
            with open(in_dir + "/" + file) as mef_results:
                for line in mef_results.readlines():
                    if line.split(',')[0][0] != '#' and line.split(',')[0][0] != 'm':
                        line = line.split(',')
                        strand = '+'
                        if line[13] > line[14]:
                            strand = '-'
                        if len(line[2]) > 0:
                            line[2] = ";" + line[2]
                        new_row = {'file_id': file.split('.')[0], 'contig_id':line[12], 'start': str(line[13]), 'end': str(line[14]), 'strand': strand, 'type': line[4], 'name': line[1] + line[2], 'identity': line[8], 'e_val': line[7], 'coverage': line[9], 'inc_rep': '-', 'mobility': '-', 'resistance': '-', 'src': 'mefinder'}
                        df = df.append(new_row, ignore_index=True)
    return df
    

def mobtyper_parser(in_dir):
    df = pd.DataFrame(columns = ['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'resistance', 'e_val', 'identity', 'coverage'])
    for folder in os.listdir(in_dir):
        if "DS" not in folder:
            if 'mobtyper_results.txt' in os.listdir(in_dir + folder):
                mobr_df = pd.read_csv(in_dir + folder + "/contig_report.txt", "\t")
                mobt_df = pd.read_csv(in_dir + folder + "/mobtyper_results.txt", "\t")
                mobt_df['cluster_id'] = mobt_df['sample_id'].str.split(':', expand=True)[1]
                mobt_df = mobt_df[['cluster_id', 'predicted_mobility']]
                mob_df = pd.merge(mobr_df, mobt_df, left_on= 'primary_cluster_id', right_on = 'cluster_id', how='inner')
                mob_df['src'] = 'MOBrecon'
                mob_df['file_id'] = mob_df['sample_id'].str.split(":").str[0]

                mob_df = mob_df[['file_id','molecule_type','contig_id', 'primary_cluster_id','rep_type(s)','predicted_mobility_y', 'src']]

                mob_df.columns = ['file_id', 'type', 'contig_id', 'name','inc_group','mobility', 'src']

                mob_df['coverage'] = ''
                mob_df['identity'] = ''
                mob_df['e_val'] = ''
                mob_df['start'] = ''
                mob_df['end'] = ''
                mob_df['strand'] = ''
                mob_df['resistance'] = ''
                mob_df = mob_df.reindex(['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'resistance', 'e_val', 'identity', 'coverage'], axis = 1)
                df = pd.concat([df, mob_df])
    return df

def make_color_map(args_df):
    random.seed(2)
    
    #Make a dict for args and resistance
    global args_resist
    args_resist = dict(zip(args_df['name'], args_df['resistance']))
    
    #One color for each resistance profile
    resist_profiles = list(set(args_df['resistance'].to_list()))
    color_categories = [",".join(sorted(x.split(','))) for x in resist_profiles]
    col_list = []
    step = int(200/len(color_categories))
    r = 130
    while len(col_list) < len(color_categories):
        g = random.randint(80, 100)
        b = random.randint(30, 90)
    
        col = colors.toColor('hsl(' + str(r) + ',' + str(g) + '%,' + str(b) + '%)')
        if col not in col_list:
            col.alpha = 0.7
            col_list.append(col)
            r += step
            
    global resist_map
    resist_map = dict(zip(color_categories, col_list))
    
    
def make_full_output(full_df, visuals, plas_list):
    full_df['resistance'] = full_df['resistance'].str.replace(';', ',')
    mge_df_plas = pd.DataFrame(columns = ['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'ARGs','resistance', 'e_val', 'identity', 'coverage'])
    mge_df_oth = pd.DataFrame(columns = ['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'ARGs','resistance', 'e_val', 'identity', 'coverage'])

    plas_file_map = pd.DataFrame(columns=['plasmid', 'file_id', 'contig_id'])
    plas_resistance_summary = pd.DataFrame(columns=['plasmid', 'file_id'])
    
    to_skip = []
    seqfeature_dict = {}
    for idx, row in full_df.iterrows():
        #If the row is a mobile element, reduce to ARGs in that node
        if row.type != 'ARG':
            #Reduce to node
            sub_df = full_df[full_df['contig_id'] == row.contig_id]

            #Reduce to args in node
            args_df = sub_df[sub_df['type'] == 'ARG']
            args = []
            resist = []

            #Reduce to same MGE type
            non_args_df = sub_df[sub_df['type'] == row['type']]

            #Aggregate all args and resistance if plasmid
            if row.type == 'plasmid':
                
                args = ",".join(set(list(args_df['name'])))
                resist = ",".join(set(list(args_df['resistance'])))
                row.start = '-'
                row.end = '-'
                row.strand = None
                name = ",".join([x for x in set(list(non_args_df['name'])) if not pd.isnull(x) ])
                mobility = ",".join([x for x in set(list(non_args_df['mobility'])) if not pd.isnull(x) ])
                inc_rep = ",".join([x for x in set(list(non_args_df['inc_rep'])) if not pd.isnull(x) ])
                identity = ",".join([str(x) for x in set(list(non_args_df['identity'])) if not pd.isnull(x) ])

                coverage = ",".join([str(x) for x in set(list(non_args_df['coverage'])) if not pd.isnull(x) ])
                e_val = ",".join([str(x) for x in set(list(non_args_df['e_val'])) if not pd.isnull(x) ])
                src = ",".join(set(list(non_args_df['src'])))

                #if ARGs found on plasmid, add to mge output and summary df
                if len(args) > 0:
                    
                    new_row = [row['file_id'], 
                               row.contig_id, 
                               row['type'],
                               name, 
                               row.start, 
                               row.end, 
                               row.strand, 
                               src,
                               inc_rep, 
                               mobility, 
                               args, 
                               resist, 
                               e_val, 
                               identity, 
                               coverage]
                    
                    mge_df_plas.loc[len(mge_df_plas)] = new_row
                    
                    #Only mobilizable plasmids
                    if 'non-mobilizable' not in mobility:
                        plas_resistance_summary.loc[len(plas_resistance_summary)] = [row['name'], row['file_id']]
                        seen = False
                        #If a common, mobilizable plasmid, create feature sets
                        if visuals > 0 and row['name'] in plas_list:
                            
                            #Create empty feature sets if file:contig not in dict
                            if row.file_id in seqfeature_dict.keys() and row.contig_id not in seqfeature_dict[row.file_id].keys():
                                    seqfeature_dict[row.file_id][row.contig_id] = [GDFeatureSet(name='ARGs'),
                                                                                   GDFeatureSet(name='Base'),
                                                                                   GDFeatureSet(name='MGEs')]
                            elif row.file_id not in seqfeature_dict.keys():
                                seqfeature_dict[row.file_id] = {}
                                seqfeature_dict[row.file_id][row.contig_id] = [GDFeatureSet(name='ARGs'),
                                                                               GDFeatureSet(name='Base'),
                                                                               GDFeatureSet(name='MGEs')]
                            else:
                                seen = True
                            if not seen:
                                #Remove repeat args
                                args_df = args_df.drop_duplicates(subset=['name', 'start', 'end'])
                                
                                #Add feature for each ARG 
                                for a_idx, a_row in args_df.iterrows():
                                    if a_row.strand == '+': 
                                        a_row.strand = 1
                                    elif a_row.strand == '-':
                                        a_row.strand = -1
                                        
                                    #Add to Arrow track
                                    seqfeature_dict[row.file_id][row.contig_id][0].add_feature(
                                        SeqFeature(
                                            FeatureLocation(
                                                int(a_row.start), int(a_row.end), a_row.strand
                                                )
                                            ),
                                        label=True, 
                                        label_position = 'start', 
                                        label_angle = 30, 
                                        label_strand = 1,
                                        name = a_row['name'], 
                                        sigil = "ARROW", 
                                        color = resist_map[args_resist[a_row['name']]], 
                                        arrowhead_length=0.25,
                                        arrowshaft_length=0.5
                                        )
                                    #Add to box track
                                    seqfeature_dict[row.file_id][row.contig_id][1].add_feature(
                                        SeqFeature(
                                            FeatureLocation(
                                                int(a_row.start), int(a_row.end), a_row.strand
                                                )
                                            ),
                                        label=True, 
                                        label_position = 'middle', 
                                        label_angle = 45, 
                                        name = a_row['name'], 
                                        sigil = "BOX", 
                                        color = resist_map[args_resist[a_row['name']]], 
                                        arrowhead_length=0.25
                                        )
        
                                #Identify MGEs on plasmid
                                other_mge_df = sub_df[sub_df['type']!= 'ARG']
                                other_mge_df = other_mge_df[other_mge_df['type']!= 'plasmid']
                                other_mge_df = other_mge_df.drop_duplicates(subset=['start', 'end'])
                                
                                #Add MGEs to featureset
                                for m_idx, m_row in other_mge_df.iterrows():
                                    seqfeature_dict[row.file_id][row.contig_id][2].add_feature(
                                        SeqFeature(
                                            FeatureLocation(
                                                int(m_row.start), int(m_row.end)
                                                )
                                            ),
                                        label=True, 
                                        label_position = 'left', 
                                        label_strand = 0,
                                        label_angle = 0, 
                                        name = m_row['type'] + ": " + m_row['name'], 
                                        color = colors.Color(1, 1, 0, alpha=0.3), 
                                        sigil = "BOX"                               
                                        )
                            plas_file_map.loc[len(plas_file_map)] = [row['name'], row.file_id, row.contig_id]

            elif row.type != 'chromosome' and idx not in to_skip:
                #For non-plasmids, find ARGS within 31kb
                name = []
                src = []
                coverage = []
                identity = []
                e_val = []

                for s_idx, s_row in args_df.iterrows():
                    if int(s_row.start) in range(int(row['start']) - 31000, int(row['start']) + 31000) or int(s_row['end']) in range(int(row['end']) - 31000, int(row['end']) + 31000):
                        [args.append(x) for x in s_row['name'].split(';')]
                        [resist.append(x) for x in s_row['resistance'].split(';')]

                #Also condense MGEs of same type with nonzero overlap
                for n_idx, n_row in non_args_df.iterrows():
                    if len(range(max(int(n_row['start']), int(row['start'])), min(int(n_row['end']), int(n_row['end']))+1)) > 0:
                        to_skip.append(n_idx)
                        name.append(n_row['name'])
                        src.append(n_row['src'])
                        row.start = min(n_row['start'], row['start'])
                        row.end = max(n_row['end'], n_row['end'])
                        coverage.append(n_row.coverage)
                        identity.append(n_row.identity)
                        e_val.append(n_row.e_val)
                if len(args) > 0:
                    new_row = [row['file_id'], row.contig_id, row['type'],",".join(set(name)), row.start, row.end, row.strand, ",".join(set(src)), row.inc_rep, row.mobility, ",".join(set(args)), ",".join(set(resist)), ",".join(e_val), ",".join(identity), ",".join(coverage)]
                    mge_df_oth.loc[len(mge_df_oth)] = new_row

    
    mge_df_plas = mge_df_plas.drop_duplicates(subset = ['file_id', 'contig_id'])
    
    #Save summary counts
    plas_resistance_summary = plas_resistance_summary.drop_duplicates()
    plas_resistance_summary = plas_resistance_summary['plasmid'].value_counts()
    plas_resistance_summary = plas_resistance_summary.to_frame()
    plas_resistance_summary.columns = ["num. strains"]
    plas_resistance_summary.index.rename('plasmid', inplace=True)
    plas_resistance_summary.to_csv(parent_dir + '/mobile_resistant_plasmids.csv')
    
    

    
    return pd.concat([mge_df_plas, mge_df_oth]).sort_values(by = ["file_id", "contig_id", "start"]), seqfeature_dict, plas_file_map, plas_resistance_summary

def make_plasmid_output(full_df):
    plas_df = pd.DataFrame(columns = ['file_id','contig_id', 'plasmid', 'inc_rep', 'ARGs', 'resistance', 'mobility','MGEs','src'])
    subsets = full_df.groupby(['file_id','contig_id']) 
    #iterate through contigs
    for idx, sub_df in subsets:
        #isolate other MGEs on each contig
        if 'plasmid' in sub_df['type'].to_list():
            plas_row = sub_df[sub_df['type'] == 'plasmid']
            sub_df = sub_df[sub_df['type'] != 'plasmid']
            assert(len(plas_row['name']) == 1)
            
            #Ignore non-mobilizable plasmids
            if(plas_row.iloc[0]["mobility"] != "non-mobilizable"):
                #If other MGEs exist, concat type, name, and location information for all into one column value
                if len(sub_df) > 0:
                    sub_df = sub_df.fillna('+')
                    sub_df['MGE_info'] = "(" + sub_df["type"] + " " + sub_df["start"].str.split('.').str[0] + ".." + sub_df["end"].str.split('.').str[0] + "(" + sub_df["strand"] + "): "+ sub_df["name"] + ")"
                    mge_info = str(list(sub_df.MGE_info))   
                else:
                    mge_info = '[]'
                new_row = [idx[0], idx[1], plas_row.iloc[0]['name'], plas_row.iloc[0]['inc_rep'], plas_row.iloc[0]['ARGs'], plas_row.iloc[0]['resistance'], plas_row.iloc[0]['mobility'], mge_info, plas_row.iloc[0]['src']]
                plas_df.loc[len(plas_df)] = new_row
    return plas_df

def make_strain_plas_summaries(plas_df):
    plas_sum = plas_df.groupby(['file_id', 'plasmid'], as_index = False).agg({'ARGs': set, 'resistance': set,  'contig_id': set, 'inc_rep': 'max', 'mobility': 'max', 'src': set})
    plas_sum = plas_sum[['file_id', 'plasmid', 'inc_rep', 'mobility','ARGs', 'resistance', 'contig_id',
            'src']]


    for idx, row in plas_sum.iterrows():
        row['ARGs'] = ",".join(set((",".join(list(row['ARGs']))).split(',')))
        row['resistance'] = ",".join(set((",".join(list(row['resistance']))).split(',')))
        row['contig_id'] = list(row['contig_id'])
        plas_sum.iloc[idx] = row
    return plas_sum

def make_visuals(input_dir, feature_dict, contig_df):
    #Reduce contig dataframe to plasmids occuring in >threshold strains
    

    count_dict = collections.defaultdict(lambda: 0)
    diagram_dict = {}
    
    #For each strain in contig_df, add tracks for all contigs
    for file in os.listdir(input_dir):

            #Extract file_id - CHANGE
            file_id = file.split('.')[0].split('_')[0]
            if file_id in contig_df['file_id'].to_list():
                print('Strain:', file_id)
                #Reduce df to relevant contigs for this strain
                strain_df = contig_df[contig_df['file_id'] == file_id]

                #Read in contigs - ONLY ONE ITERATOR PER FILE!!!
                strain_seq = SeqIO.parse(input_dir + "/" + file, "fasta")
                
                #Create tracks for relevant contigs
                for record in strain_seq:
                    if record.id in strain_df['contig_id'].to_list():
                        
                        #Idenitfy plasmid corresponding to contig
                        strain_plas_list = set(strain_df[strain_df['contig_id'] == record.id]['plasmid'].to_list())
                        for plasmid in strain_plas_list:
                            print('Current plasmid:', plasmid)
                            #Create contig tracks and add feature sets
                            
                            arg_track = GDTrack(name= file_id+ ":" + record.id,  
                                                greytrack=1, 
                                                greytrack_labels=1, 
                                                scale = 1, 
                                                axis_labels=1)
                            arg_track.add_set(feature_dict[file_id][record.id][0])
                            arg_track.add_set(feature_dict[file_id][record.id][2])
                            
                          #  mge_track = GDTrack(name= file_id+ ":" + record.id + "–MGEs" , start=0, end=len(record) )
                           # mge_track.add_set(feature_dict[file_id][record.id][2])
                            
                            #base_ft_track = GDTrack(name=file_id+ ":" + record.id, start=0, end=len(record), greytrack=1)
                            #base_ft_track.add_set(feature_dict[file_id][record.id][1])
                            
                            #Add contig tracks to plasmid diagram
                            top = count_dict[plasmid]
                            print('Adding to ' + str(top) + ' level')
                            if plasmid not in diagram_dict.keys():
                                diagram_dict[plasmid] = GDDiagram(name=plasmid)
                            
                        #    diagram_dict[plasmid].add_track(mge_track, track_level=top)
                            #diagram_dict[plasmid].add_track(base_ft_track, track_level=top+1)
                            diagram_dict[plasmid].add_track(arg_track, track_level=top)
                            count_dict[plasmid] += 2
                    
                
                #Add a blank space at top of each updated diagram
                for pls in pd.unique(strain_df['plasmid']):
                    count_dict[pls] += 2
                    
    for key, value in diagram_dict.items():

        value.draw(format="linear", pagesize=(count_dict[key] * cm,max(count_dict.values()) * cm), fragments=1, start=value.range()[0], end = value.range()[1])
        value.write(args.d + "/plasmid_diagrams/" + key + ".pdf", "PDF")
        
       
        
def main():
    print("Running")
    
    prog_list = os.listdir(parent_dir)
    full_df = pd.DataFrame(columns = ['file_id','contig_id', 'type', 'name','start', 'end', 'strand', 'src', 'inc_rep', 'mobility', 'resistance', 'e_val', 'identity', 'coverage'])
    cur_df = None
    for prog in prog_list:
        print(prog)
        if prog in ['card', 'ncbi', 'resfinder', 'plasmidfinder']:
            in_file = parent_dir + '/' + prog + '/' + prog + '_results.txt'
            cur_df = abricate_parser(in_file, prog)

        if prog in ['mefinder']:
            in_dir = parent_dir + '/' + prog + "/"
            cur_df = mefinder_parser(in_dir)

        if prog in ['mobtyper']:
            in_dir = parent_dir + '/' + prog + "/"
            cur_df = mobtyper_parser(in_dir)

        full_df = pd.concat([full_df, cur_df])
    
    #full_df= pd.read_csv(parent_dir + "/gbk_maker_input.csv")
    full_df = full_df.sort_values(by = ["file_id", "contig_id", "start"])
    full_df = full_df.drop_duplicates()
    full_df.to_csv(args.d + "/compiled_output/all_hits.csv")
    
    args_df = full_df[full_df['type'] == 'ARG']
    args_df = args_df.drop_duplicates()
    make_color_map(args_df)
    print(resist_map)
    #List of all mobilizable plasmids occuring in >= threshold strains
    plas_list = full_df[full_df['type'] == 'plasmid'].drop_duplicates(subset=['name', 'file_id'])
    plas_list = plas_list[plas_list['mobility'] != 'non-mobilizable'].groupby('name', as_index = False).agg({'file_id': list})
    plas_list = plas_list[plas_list['file_id'].map(len) >= threshold]
    plas_list = list(plas_list['name'])

    
    mge_df, feature_dict, contig_df, plas_summary = make_full_output(full_df, args.v, plas_list)
    mge_df.to_csv(args.d + "/compiled_output/all_mges.csv")
    
    plasmid_df = make_plasmid_output(mge_df)
    plasmid_df.to_csv(args.d + "/compiled_output/all_plasmids.csv")
        
    plasmids_by_strain = make_strain_plas_summaries(plasmid_df)
    plasmids_by_strain.to_csv(args.d + "/compiled_output/plasmids_by_strain.csv")
    
    
    if args.v > 0:
        vis_to_make = plas_summary[plas_summary['num. strains'] >= threshold].index.to_list()
        contig_df = contig_df[contig_df['plasmid'].isin(vis_to_make)]
        print(pd.unique(contig_df['plasmid']))

        make_visuals(input_dir, feature_dict, contig_df)

if __name__ == '__main__':
    main()
    