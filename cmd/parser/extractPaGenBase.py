'''
Database_web_link = http://bioinf.xmu.edu.cn/PaGenBase/index.jsp
'''
import pandas as pd
import glob
import os


grch37 = "/Source/scripts/data/GRCh37.gene.bed"
file = "/Source/scripts/data/geneTpye/Blood.txt"
filedir = "/Source/scripts/data/geneTpye/"


# output_dir = os.path.dirname(file)
output_dir = '/Source/scripts/data/deal'


grch_df = pd.read_csv(grch37, sep="\t", header=None, names=[
                      'chr', 'start', 'end', 'strand', 'gene'])


def main():
    extactDir(filedir)
    # extractFile(file)


def extactDir(filedir):
    print("extracting files in directory: " + filedir)
    for file in glob.glob(os.path.join(filedir, "*.txt")):
        extractFile(file)


def extractFile(file, housekeeping=True, spm_min=0.5, spm_max=1.0):
    file_name = os.path.basename(file)
    filename_without_extension = os.path.splitext(file_name)[0]
    data_df = pd.read_csv(file, sep="\t", skiprows=7)

    if len(data_df.columns) == 5:
        data_df = filterSPM(data_df, spm_min, spm_max)
        print("filtering SPM of " + filename_without_extension)
        extractBaseToTSS(
            data_df, grch_df, output_dir, filename_without_extension)
    else:
        if housekeeping:
            print("filtering housekeeping of " + filename_without_extension)
            hk_df = filterGeneType(data_df, housekeeping=True)
            extractBaseToTSS(
                hk_df, grch_df, output_dir, filename_without_extension)
            low_df = filterDPM(data_df, dpm_min=0.7, dpm_max=1.0)
            extractBaseToTSS(
                low_df, grch_df, output_dir, filename_without_extension, low=True)


def extractBaseToTSS(data_df, grch_df, output_dir, filename, out_gene=False, low=False):
    gene_df = grch_df[grch_df['gene'].isin(data_df['Gene Symbol'])]
    if out_gene:
        gene_df.to_csv(os.path.join(output_dir, filename +
                       "_Gene.bed"), sep="\t", header=None, index=None)
    tss_bed = extractTSS(gene_df)
    if low:
        out_path = os.path.join(output_dir, filename + "_TSS_low.bed")
    else:
        out_path = os.path.join(output_dir, filename + "_TSS.bed")
    tss_bed.to_csv(out_path, sep="\t", header=None, index=None)
    print('save to ' + out_path)


def filterGeneType(df, specific=False, housekeeping=False, selective=False, repressed=False):
    df_new = df.copy()
    if specific:
        df_new = df_new[df_new['Specific'] == 'Y']
    if housekeeping:
        df_new = df_new[df_new['Housekeeping'] == 'Y']
    if selective:
        df_new = df_new[df_new['Selective'] == 'Y']
    if repressed:
        df_new = df_new[df_new['Repressed'] == 'Y']

    return df_new


def filterDPM(df, dpm_min=0, dpm_max=1.0):
    return df[(df['DPM/PubMed ID'] >= dpm_min) & (df['DPM/PubMed ID'] <= dpm_max)]


def filterSPM(df, spm_min=0.5, spm_max=1.0):
    return df[(df['SPM'] >= spm_min) & (df['SPM'] <= spm_max)]


def extractTSS(df_bed):
    for index in df_bed.query("strand=='+'").index:
        df_bed.at[index, 'start'] = df_bed.loc[index, 'start'] - 1000
        df_bed.at[index, 'end'] = df_bed.loc[index, 'start'] + 2000
    for index in df_bed.query("strand=='-'").index:
        df_bed.at[index, 'start'] = df_bed.loc[index, 'end'] - 1000
        df_bed.at[index, 'end'] = df_bed.loc[index, 'start'] + 2000
    df_bed = df_bed[(df_bed['strand'] == '+') | (df_bed['strand'] == '-')]
    return df_bed


if __name__ == "__main__":
    main()
