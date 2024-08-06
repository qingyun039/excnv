#!/usr/bin/env python
import os
import sys
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
import plotly.graph_objects as go
#from scipy.interpolate import make_interp_spline
# from cnvlib.cnary import CopyNumArray as CNA

def chromosome_collections(df, y_positions, height, **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict or float
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df = df.copy()
        df.loc[:, 'width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        #print(chrom)
        if isinstance(y_positions, dict):
            yrange = (y_positions[chrom], height)
        else:
            yrange = (y_positions, height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], edgecolor='black', **kwargs)
    if del_width:
        del df['width']

def load_ideo_data(cytoband):
    """
    载入核型信息，可以从 UCSC's Table Browser 中的 'cytoBandIdeo' 表获取.
    Parameters
    ----------
    cytoband : file
        加载的核型信息文件
    """
    ideo = pd.read_table(cytoband, skiprows=1, names=['chrom', 'start', 'end', 'name', 'gieStain'])
    color_lookup = {
        'gneg': (1., 1., 1.),
        'gpos25': (.6, .6, .6),
        'gpos50': (.4, .4, .4),
        'gpos75': (.2, .2, .2),
        'gpos100': (0., 0., 0.),
        'acen': (.8, .4, .4),
        'gvar': (.8, .8, .8),
        'stalk': (.9, .9, .9),
    }
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])
    return ideo

def genome_whole_view(ideo, cnarr, figout='genome_whole_view', min_cnv=5000000, title=None):
    """
    绘制基因组CNV的全景图
    Parameters
    ----------
    ideo : pandas.DataFrame
        代表基因组核型的数据框，必须包含['chrom', 'start', 'end', 'color']列
    cnarr : pandas.DataFrame
        代表拷贝数变异的数据框，必须包含['chromosome', 'start', 'end', 'type']列
    fig : str
        输出图片路径
    min_cnv : int
        最小绘制的CNV跨度
    """
    chrom_height = 1
    chrom_spacing = 4
    chromosome_list = [f'chr{i}' for i in list(range(1,23)) + ['X', 'Y']]

    ybase = 0
    chrom_ybase = {}
    chrom_centers = {}
    cnv_ybase = {}

    for chrom in chromosome_list[::-1]:
        chrom_ybase[chrom] = ybase
        chrom_centers[chrom] = ybase + chrom_height / 2.
        cnv_ybase[chrom] = ybase + chrom_height + chrom_spacing / 4.
        ybase += chrom_height + chrom_spacing

    figsize = (16, 16)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    # 绘制染色体
    ideo = ideo[ideo.chrom.isin(chromosome_list)]
    for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
        ax.add_collection(collection)

    # 绘制拷贝数变异
    for _, cnv in cnarr.iterrows():
        if cnv.end - cnv.start < min_cnv:
            continue
        color = 'red' if cnv.type.lower() == 'dup' else 'blue'
        ypos = cnv_ybase[cnv.chromosome]
        ax.plot([cnv.start, cnv.end], [ypos, ypos], color=color, linestyle='-', label=cnv.type)

    dup_line = Line2D([0], [0], label='dup', color='red')
    del_line = Line2D([0], [0], label='del', color='blue')
    ax.legend(handles=[dup_line, del_line], loc='lower right')

    ax.set_yticks([chrom_centers[i] for i in chromosome_list])
    ax.set_yticklabels(chromosome_list)
    ax.axis('tight')

    if title:
        fig.suptitle(title)
    fig.savefig(figout)
    plt.close()

def chromosome_cnv_html(cnr, cns=None, ref_copy=2, figout='genome.html'):
    """
    所有染色体上CNV的散点图，html可放大缩小

    Parameters
    ----------
    cnr : pandas.DataFrame
        染色体bin区间的拷贝数信息, 必须包含['chromosome', 'start', 'end', 'log2']列
    cns : pandas.DataFrame
        染色体bin区间的拷贝数信息, 必须包含['chromosome', 'start', 'end', 'log2']列
    """
    outdir = '.'
    figout = os.path.join(outdir, f'{figout}')
    f = open(figout, 'w')
    for chrom in cnr.chromosome.unique():
        chrom_cnr = cnr[cnr.chromosome == chrom]
        
        fig = go.Figure()
        copy_num = (np.exp2(chrom_cnr.log2) * ref_copy).round(2)
        copy_num[copy_num > 4] = 4
        fig.add_trace(go.Scatter(x=chrom_cnr.start, y=copy_num, mode='markers', name='bin'))
        if cns is not None:
            chrom_cns = cns[cns.chromosome == chrom]
            xs = chrom_cns[['start', 'end']].to_numpy().flatten()
            ys = (np.exp2(chrom_cns[['log2', 'log2']].to_numpy().flatten()) * ref_copy).round(2)
            ys[ys > 4] = 4
            fig.add_trace(go.Scatter(x=xs, y=ys, fillcolor='red', mode='lines'))

        f.write(f"<h3>{chrom}</h3>")
        fig.write_html(f)

    f.close()

def chromosome_cnv_view(ideo, cnr, cns=None, ref_copy=2, figout='chromosome_cnv', title=None):
    """
    所有BIN的CN值在染色体上的分布
    Parameters
    ----------
    ideo : pandas.DataFrame
        染色体核型信息
    cnr : pandas.DataFrame
        染色体bin区间的拷贝数信息, 必须包含['chromosome', 'start', 'end', 'log2']列
    cns : pandas.DataFrame
        染色体bin区间的拷贝数信息, 必须包含['chromosome', 'start', 'end', 'log2']列
    """
    figsize = (16, 5)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)

    ax.add_collection(next(chromosome_collections(ideo, -1, 0.5)))
    copy_num = (np.exp2(cnr.log2) * ref_copy).round(2)
    copy_num[copy_num > 4] = 4
    ax.scatter(cnr.start, copy_num)
    if cns is not None:
        xs = cns[['start', 'end']].to_numpy().flatten()
        ys = (np.exp2(cns[['log2', 'log2']].to_numpy().flatten()) * ref_copy).round(2)
        ys[ys > 4] = 4

        #xnew = np.linspace(xs.min(), xs.max(), 1000)
        #spl = make_interp_spline(xs, ys, k=3)
        #ynew = spl(xnew)
        ax.plot(xs, ys, color='red')

    if title:
        fig.suptitle(title)
    fig.savefig(figout)
    plt.close()

def gene_cnv_view(cnr, ref_copy=2, figout='gene_cnv', title=None):
    """
    绘制基因上的CNV分布图
    Parameters
    ----------
    cnr : pandas.DataFrame
        染色体bin区间的拷贝数信息, 必须包含['chromosome', 'start', 'end', 'log2']列
    figout : str
        输出图片路径
    """
    figsize = (8, 4)
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    if cnr.gene.str.contains('_').all():
        exons = cnr.gene.str.split('_', expand=True).iloc[:, -1]
    else:
        exons = [i+1 for i in range(len(cnr))]
    print(exons)
    ax.bar(exons, np.exp2(cnr.log2) * ref_copy)
    ax.axhline(y=2, color='black', linestyle='--')
    ax.axhline(y=1, color='red', linestyle='--')
    ax.axhline(y=3, color='red', linestyle='--')
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xticks(exons)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.bottom.set_visible(False)
    ax.tick_params(bottom=False)
    ax.set_xlabel('Exon Rank')
    ax.set_ylabel('Copy Number')
    if title:
        fig.suptitle(title)
    fig.savefig(figout)
    plt.close()

def test():
    cnr = pd.read_table(sys.argv[1])
    cns = pd.read_table(sys.argv[2])
    chromosome_cnv_html(cnr, cns, figout='genome.html')

def parse_args():
    parser = argparse.ArgumentParser(description='CNV检测结果绘图')
    parser.add_argument(
            '-i', '--ideo',
            required = True,
            help = '指定ideo文件'
            )
    parser.add_argument(
            '-o', '--outdir',
            help = '输出图像目录'
            )
    parser.add_argument(
            'cnrfile',
            help = '指定cnr文件'
            )
    parser.add_argument(
            'cnsfile',
            help = '指定cns文件'
            )
    parser.add_argument(
            'tsvfile',
            help = '指定tsv文件'
            )
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    sample = os.path.basename(args.tsvfile)[:-4]
    if args.outdir:
        outdir = args.outdir
    else:
        outdir = sample
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    ideo = load_ideo_data(args.ideo)
    cnr = pd.read_table(args.cnrfile)
    cns = pd.read_table(args.cnsfile)
    tsv = pd.read_table(args.tsvfile, skiprows=6)
    
    # genome
    figout = os.path.join(outdir, 'genome')
    genome_whole_view(ideo, tsv, min_cnv=1000, figout=figout)

    # chrom
    for chrom in cnr.chromosome.unique():
        chrom_ideo = ideo[ideo.chrom == chrom]
        chrom_cnr = cnr[cnr.chromosome == chrom]
        chrom_cns = cns[cns.chromosome == chrom]
        figout = os.path.join(outdir, f'{chrom}')
        chromosome_cnv_view(chrom_ideo, chrom_cnr, chrom_cns, figout=figout)

    # gene
    tsv = tsv[tsv.probes < 100]
    genes = []
    for genelist in tsv.gene.str.split(','):
        genes.extend(genelist)
    for gene in pd.Series(genes).str.replace('_.+?$', '', regex=True).unique():
        gene_cnr = cnr[cnr.gene.str.match(gene)]
        exons = []
        for genelist in gene_cnr.gene:
            genes = genelist.split(',')
            if len(genes) > 1:
                for exon in genes:
                    if exon.startswith(gene+'_'):
                        exons.append(exon)
                        break
                else:
                    exons.append(gene+'_0')
            else:
                exons.append(genes[0])
        gene_cnr = gene_cnr.copy()
        gene_cnr.loc[:, 'gene'] = exons

        figout = os.path.join(outdir, f'{gene}')
        gene_cnv_view(gene_cnr, figout=figout)



if __name__ == '__main__':
    test()
