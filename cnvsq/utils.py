import numpy as np


def filter_normal(cnarr, is_female):
    '''过滤掉正常的区间'''
    if is_female:
        pairsome = (cnarr['chromosome'] != 'chrY')
        #is_normal = ((pairsome & (cnarr['cn'] == 2)) | ((~pairsome) & (cnarr['cn'] == 0)))
        is_normal = ((pairsome & (cnarr['cn'] == 2)) | (~pairsome))
    else:
        autosome = (~cnarr['chromosome'].isin(['chrX', 'chrY']))
        is_normal = ((autosome & (cnarr['cn'] == 2)) | ((~autosome) & (cnarr['cn'] == 1)))

    return cnarr[~is_normal]

def add_cnv_type(cnarr, is_female):
    '''CNV是DEL还是DUP'''
    cnarr['type'] = np.where(cnarr['cn'] == 2, 'WILD', 'CNV')
    if is_female:
        cnarr['type'] = np.where(cnarr['cn'] > 2, 'DUP', 'DEL')
    else:
        cnarr['type'] = np.where(cnarr['cn'] > 2, 'DUP', 'DEL')
        for chrom in ('chrX', 'chrY'):
            cnarr.data.loc[(cnarr['chromosome'] == chrom) & (cnarr['cn'] > 1) ,'type'] = 'DUP'
            cnarr.data.loc[(cnarr['chromosome'] == chrom) & (cnarr['cn'] < 1) ,'type'] = 'DEL'

    return cnarr

