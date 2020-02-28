import pandas as pd
import gzip
import subprocess
import os

import sys
sys.path.append("/work/isabl/home/zhouy1/mergeSVvcf")
from mergesvvcf import mergedfile

def cat1_filtering_noSUPP(row):
    row=list(row)
    count=str(row[105:120]).count("yes")
    return count

def merge(sample):
    print(sample)
    vcf = f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/uk_all_sv_0227/{sample}/merged_svs.all.vcf.gz"
    df = pd.read_csv(f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/uk_all_sv_0227/{sample}/{sample}_vs_I-H-108298-N1-1-D1-1.annotated.flagged.output.tsv", sep='\t')
    df = df[df['AnnotSV type']=='full']
    df['Cat1_noSUPP']=df.apply(lambda row: cat1_filtering_noSUPP(row), axis=1)
    df=df[df['Cat1_noSUPP']==0]
    df['end1'] = df['end1'].astype(str)
    df['end2'] = df['end2'].astype(str)

    # *** any additional filters goes here ***
    with gzip.open(vcf, 'rb') as f:
        content = f.readlines()
        content = [i.decode('ascii') for i in content]
        header = [i for i in content if i.startswith('#')]
        content = [i for i in content if not i.startswith('#')]

        match = []
        for index, row in df.iterrows():
            for i in content:
                start = i.split('\t')[1]
                end = i.split('\t')[7].split(';')[0].strip('END=')
                if ((row['end1'] == str(start)) & (row['end2'] == str(end)))| ((row['end1'] == str(end)) & (row['end2'] == str(start))):
                    match.append(i)
                    df.loc[index, 'match'] = 1
                    continue

    with open(f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/brass_filtered_vcf/{sample}.brass.vcf", 'w') as f:
        f.write(''.join(header))
        for i in match:
            f.write(i)

    vcf = f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/brass_filtered_vcf/{sample}.brass.vcf"
    cmd = ['bcftools', 'sort', vcf, '-o', vcf+'.sorted.vcf']
    subprocess.check_call(cmd)
    cmd = ['bgzip', vcf+'.sorted.vcf']
    subprocess.check_call(cmd)
    cmd = ['tabix', '-p', 'vcf', vcf+'.sorted.vcf.gz']
    subprocess.check_call(cmd)
    
    
    svaba = f"/work/isabl/home/gutierj2/uk_all/svABA_unmatch/universal_normal/joe/{sample}.svaba.somatic.sv.vcf.gz"
    gridss = f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/gridss_somatic/{sample}_gridss_somatic.vcf"
    brass = f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/brass_filtered_vcf/{sample}.brass.vcf.sorted.vcf.gz"
    mergedfile.merge(
        filenames = [brass, svaba, gridss],
        programs = ['brass', 'svaba', 'gridss'],
        slop=100,
        verbose=False,
        min_num_callers=2,
        forceSV=True,
        outfile = f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/merged/universal_normal/{sample}.merged.vcf"
    )
    
    os.makedirs(f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/annotated/universal_normal/{sample}", exist_ok=True)
    cmd = [
        "singularity", "exec",
        "--bind", "/work:/work",
        "--bind", "/ifs:/ifs",
        "--bind", "/juno:/juno",
        "/work/isabl/home/zhouy1/images/singularity/toil_unmatched_cnvkit_0.4.1.simg",
        "AnnotSV.tcl",
        "-SVinputFile", f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/merged/universal_normal/{sample}.merged.vcf",
        "-bedtools", "/usr/local/bin/bedtools",
        "-outputFile", f"/work/isabl/home/zhouy1/unmatechedSVpipeline/test/unmatched_benchmark/sv_triplecaller_merge/annotated/universal_normal/{sample}/{sample}.merged.annotated.tsv",
        "-annotationFolder", "/work/isabl/home/zhouy1/annotsv/2.1_patch/AnnotSV/Annotations"
    ]
    subprocess.check_call(cmd)


def main():
    sample = sys.argv[1]
    merge(sample)

if __name__ == '__main__':
    main()
