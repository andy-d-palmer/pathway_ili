import io
#import pandas as pd
#import numpy as np
#
#import sm_annotation_utils

def generate_ili_csv(groups, database='HMDB', fdr=0.2, map_id='main'):
    print ">>>>>", groups
    dataset_names, annotations = get_annotations(groups, database, fdr)
    ili_coords, kegg_compounds = get_kegg_compounds(map_id)
    annotations_df = format_annotations(kegg_compounds, annotations, dataset_names, ili_coords)
    ili_out = partition_by_group(annotations_df, groups, ili_coords)
    s_buf = io.BytesIO()
    ili_out.to_csv(s_buf)
    s_buf.seek(0)
    return s_buf

def partition_by_group(annotations_df, groups, ili_coords):
    import pandas as pd
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
    print ili_out.head()
    return ili_out

def format_annotations(kegg_compounds, annotations, dataset_names, ili_coords):
    from pyMSpec.pyisocalc.tools import normalise_sf
    import pandas as pd
    kegg_sf = [normalise_sf(sf) for sf in kegg_compounds['sf'].values]
    annotations_df = pd.DataFrame([[sf in ds_annotations for sf in kegg_sf] for ds_annotations in annotations],
                                  index=dataset_names, columns=kegg_compounds.index).transpose()
    annotations_df = annotations_df.loc[ili_coords.index].fillna(0).astype(int)
    return annotations_df

def get_kegg_compounds(map_id):
    import pandas as pd
    kegg_compounds = pd.read_csv(open('/media/embl/palmer/databases_in/kegg_compounds_all.csv'), sep=',',
                                 names=['id', 'name', 'sf', 'mz', 'map'], index_col=0)
    kegg_compounds = kegg_compounds[kegg_compounds.mz > 0.].dropna()
    ili_coords = pd.read_csv(open('/media/embl/palmer/tmp/KEGG_EC_metabolitemap_bw.csv'), index_col=0)
    return ili_coords, kegg_compounds

def get_annotations(groups, database, fdr):
    from sm_annotation_utils.sm_annotation_utils import SMInstance
    import numpy as np
    from pyMSpec.pyisocalc.tools import normalise_sf
    config = '/home/palmer/Documents/ipython_notebooks/alpha_es_mb_config.json'
    remote_instance = SMInstance(config)
    annotations = []
    dataset_names = []
    for group in groups:
        for target_dataset in groups[group]:
            dataset = remote_instance.dataset(dataset_name=target_dataset)
            dataset_names.append(target_dataset)
            annotations.append([normalise_sf(sf_a[0]) for sf_a in dataset.annotations(database=database, fdr=fdr)])
            print "-", target_dataset, "-", len(annotations[-1])

        print group, [len(annotations[ii]) for ii in range(len(annotations)-len(groups[group]), len(annotations))]
    return dataset_names, annotations
