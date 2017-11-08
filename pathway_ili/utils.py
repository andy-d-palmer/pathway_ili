import io
import pandas as pd
import numpy as np
from pyMSpec.pyisocalc.tools import normalise_sf
from sm_annotation_utils.sm_annotation_utils import SMInstance
import sys
PY3 = sys.version_info[0] == 3

#import sm_annotation_utils

def generate_ili_csv(groups, map_id, database='HMDB', fdr=0.2):
    dataset_names, annotations = get_annotations(groups, database, fdr)
    ili_coords, kegg_compounds = get_kegg_compounds(map_id)
    annotations_df = format_annotations(kegg_compounds, annotations, dataset_names, ili_coords)
    ili_out = partition_by_group(annotations_df, groups, ili_coords)
    #if PY3:
    #    print('PY3')
    #    s_buf = io.StringIO()
    #else:
    s_buf = io.BytesIO()
    s_buf.write(str(ili_out.to_csv(encoding='utf-8')).encode())
    s_buf.seek(0)
    print('returning csv')
    return s_buf

def partition_by_group(annotations_df, groups, ili_coords):
    all_annotations = pd.DataFrame(annotations_df.sum(axis=1), columns=['all annotations', ])
    gg = []
    for ii, group1 in enumerate(groups):
        g1 = annotations_df[groups[group1]].mean(axis=1)
        gg.append(g1.rename(group1))
    for ii, group1 in enumerate(groups):
        g1 = annotations_df[groups[group1]].mean(axis=1)
        for jj, group2 in enumerate(groups):
            if not ii == jj:
                g2 = annotations_df[groups[group2]].mean(axis=1)
                gn = " - ".join([group1, group2])
                gg.append((g1 - g2).rename(gn))
    groups_df = pd.concat(gg, axis=1)
    ili_out = pd.concat([ili_coords, annotations_df, all_annotations, groups_df], axis=1)
    print (ili_out.head())
    return ili_out

def format_annotations(kegg_compounds, annotations, dataset_names, ili_coords):
    kegg_sf = [normalise_sf(sf) for sf in kegg_compounds['sf'].values]
    annotations_df = pd.DataFrame([[sf in ds_annotations for sf in kegg_sf] for ds_annotations in annotations],
                                  index=dataset_names, columns=kegg_compounds.index).transpose()
    print(annotations_df.index)
    print(ili_coords.index)
    annotations_df = annotations_df.loc[ili_coords.index].fillna(0).astype(int)
    return annotations_df

def get_kegg_compounds(map_id):
    kegg_compounds = pd.read_csv(open('/media/embl/palmer/databases_in/kegg_compounds_all.csv'), sep=',',
                                 names=['id', 'name', 'sf', 'mz', 'map'], index_col=0)
    kegg_compounds = kegg_compounds[kegg_compounds.mz > 0.].dropna()
    ili_coords = pd.read_csv(open('/media/embl/palmer/tmp/{}.csv'.format(map_id)), index_col=0)
    return ili_coords, kegg_compounds

def get_map_filename(mapid):
    return '/media/embl/palmer/tmp/{}.png'.format(mapid)

def get_annotations(groups, database, fdr):
    config = {}#'/home/palmer/Documents/ipython_notebooks/alpha_es_mb_config.json'
    remote_instance = SMInstance()
    annotations = []
    dataset_names = []
    for group in groups:
        for target_dataset in groups[group]:
            print('getting this target:', target_dataset, "\n")
            dataset = remote_instance.dataset(name=target_dataset)
            dataset_names.append(target_dataset)
            annotations.append([normalise_sf(sf_a[0]) for sf_a in dataset.annotations(database=database, fdr=fdr)])
            print ("-", target_dataset, "-", len(annotations[-1]))
        print (group, [len(annotations[ii]) for ii in range(len(annotations)-len(groups[group]), len(annotations))])
    return dataset_names, annotations


def get_FA_formula(chain_length, n_dbl_bonds, modifier=''):
    modifier_lookup = {'': (0, 0, 0), 'O-': (-1, +2, 0), 'dO-': (-2, +4, 0), 'tO-': (-3, +6, 0),
                       "P-": (-1, 0, 0),
                       "d": (0, 3, 1), "t": (1, 3, 1)}  # (O,H,N)
    n_oxygen = 1 + modifier_lookup[modifier][0]
    n_carbon = int(chain_length)
    n_hydrogen = n_carbon * 2 - 2 * int(n_dbl_bonds)
    n_nitrogen = 0 + modifier_lookup[modifier][2]
    sf ="".join("{}{}".format(e,n) for e,n in zip(["C", "H", "N", "O"], [n_carbon, n_hydrogen, n_nitrogen, n_oxygen]) if n>0)
    return sf


def getRGroups():
    r_groups = []
    for c in [14, 16]:#, 18, 20]:
        for dbl in [0,1]:#,2]:
            mf = get_FA_formula(c, dbl)
            r_groups.append(mf)
    return r_groups


def replaceR(mf):
    # Replace all instances of 'R' with templates in R_GROUPS
    import re
    mf_split = re.split('([A-Z][a-z]*)', mf)
    el = mf_split[1::2]
    num = [int(n) if n != "" else 1 for n in mf_split[2::2]]
    mf_out=[]
    try:
        if num[el.index('R')] > 4:
            raise ValueError('too many substitutions {}'.format(num[el.index('R')]))
        num[el.index('R')] -= 1
        mf_base = "".join('{}{}'.format(e,n) for e, n in zip(el,num) if n>0)
        mf_out.append( [replaceR(mf_base + r) for r in R_GROUPS])
    except:
        mf_out.append(normalise_sf("".join('{}{}'.format(e,n) for e, n in zip(el,num) if n>0)))
    return np.asarray(mf_out).flatten()


def parse_kegg(kegg_fn):
    from pyMSpec.pyisocalc.pyisocalc import InvalidFormulaError, ParseError, parseSumFormula
    kegg_compounds = pd.read_csv(kegg_fn, sep='\t', index_col=0)
    kegg_compounds = kegg_compounds.loc[kegg_compounds.FORMULA.isnull()==False] #drop entries without a mf
    id_to_sf = {}
    for row in kegg_compounds.iterrows():
        db_id = row[0]
        mf = row[1]['FORMULA']
        if ")n" in mf:
            continue
        mf = mf.replace('. ', "").replace("(R1)","R").replace("(R2)","R")
        if "(" in mf:
            print ('not parsing {}'.format(mf))
            continue
        try:
            parsed_mf = [parseSumFormula(mf).__unicode__(),]
        except (ValueError, ParseError, InvalidFormulaError) as e:
            try:
                parsed_mf = replaceR(mf)
            except ValueError as e:
                print (e, mf)
        id_to_sf[db_id] = parsed_mf
    sf_to_id = {}
    for k, v in id_to_sf.iteritems():
        for _v in v:
            sf_to_id.setdefault(_v, []).append(k)
    return id_to_sf, sf_to_id
