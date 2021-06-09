import pandas as pd

target_name = config['target_path'].split('.txt')[-2].split('/')[-1]

rule all:
    input:
        f"output/{target_name}.hg38.txt"

def write_vcf(target, meta, out_path):
    with open(out_path, 'w') as write_obj:
        write_obj.write(meta)
        target.to_csv(write_obj, sep="\t", index=False)

def read_meta(meta_path):
    with open(meta_path, 'r') as read_obj:
        lines = read_obj.readlines()
    return ''.join(lines)

rule pseudo_vcf:
    input:
        target_path=config['target_path'],
        header_path="resources/vcf_meta_data.txt"

    output:
        temp(f"temp/{target_name}.vcf")

    run:
        target = pd.read_csv(str(input.target_path), sep='\t', header=None,names=['chrom:pos:ref:alt', 'rsid'])

        # liftover needs chr prefix
        target['#CHROM'] = target['chrom:pos:ref:alt'].str.split(':', expand=True)[0]
        target['#CHROM'] =  'chr' + target['#CHROM'].astype(str)

        target['POS'] = target['chrom:pos:ref:alt'].str.split(':',expand=True)[1]
        target['REF'] = target['chrom:pos:ref:alt'].str.split(':',expand=True)[2]
        target['ALT'] = target['chrom:pos:ref:alt'].str.split(':',expand=True)[3]
        target['ID'] = target['chrom:pos:ref:alt'] + ';' + target['rsid']

        target['QUAL'] = '.'
        target['FILTER'] = 'PASS'
        target['INFO'] = '.'
        target['FORMAT'] = 'GT'

        # make dummy sample
        target['dummy'] = '0|0'

        target_order = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'dummy']
        target_ordered = target[target_order]
        meta = read_meta(meta_path=str(input.header_path))
        write_vcf(target=target_ordered,meta=meta,out_path=str(output))

rule get_ref:
    output:
        ref_path="resources/hg38.fasta"

    params:
        zipped="resources/hg38.fa.gz"

    shell:
        """
        wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz -P resources
        zcat {params.zipped} > {output.ref_path}
        """

rule get_chain:
    output:
        chain_path="resources/hg19ToHg38.over.chain.gz"

    shell:
        """
        wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz -P resources
        """

rule seq_dict_path:
    input:
        ref_path="resources/hg38.fasta"

    output:
        dict_path="resources/hg38.fasta.dict"

    shell:
        """
        picard \
        CreateSequenceDictionary \
        R={input.ref_path} \
        O={output.dict_path} 
        """

rule liftover:
    input:
        target_vcf_path=f"temp/{target_name}.vcf",
        ref_path="resources/hg38.fasta",
        dict_path="resources/hg38.fasta.dict",
        chain_path="resources/hg19ToHg38.over.chain.gz"

    output:
        rejected_vcf_path=temp(f"output/{target_name}.rejected.vcf"),
        lifted_vcf_path=temp(f"output/{target_name}.hg38.vcf"),
        lifted_vcf_index_path=temp(f"output/{target_name}.hg38.vcf.idx")

    shell:
        """
         picard -Xmx32g LiftoverVcf \
         I={input.target_vcf_path} \
         O={output.lifted_vcf_path} \
         CHAIN={input.chain_path} \
         REJECT={output.rejected_vcf_path} \
         R={input.ref_path}
        """


def get_header_row(vcf_path):
    # detect header row, cna be inconsistent
    with open(vcf_path, 'r') as read_obj:
        for counter, line in enumerate(read_obj):
            line_clean = line.strip('\n').strip('\r').split('\t')
            if '#CHROM' in line_clean:
                return counter

def clean_lifted(lifted_raw):
    join_cols = ['#CHROM', 'POS', 'REF', 'ALT']
    lifted_raw['hg38_chrom_pos_ref_alt'] = lifted_raw[join_cols].agg(':'.join, axis=1)
    lifted_raw['hg19_chrom_pos_ref_alt'] = lifted_raw['ID'].str.split(';', expand=True)[0]
    lifted_raw['rsid'] = lifted_raw['ID'].str.split(';', expand=True)[1]
    lifted_raw['issue_flag'] = 'None'
    keep_cols = ['hg38_chrom_pos_ref_alt', 'hg19_chrom_pos_ref_alt', 'rsid', 'issue_flag']
    return lifted_raw[keep_cols]

def clean_rejected(rejected_raw):
    rejected_raw['hg38_chrom_pos_ref_alt'] = 'None'
    rejected_raw['hg19_chrom_pos_ref_alt'] = rejected_raw['ID'].str.split(';', expand=True)[0]
    rejected_raw['rsid'] = rejected_raw['ID'].str.split(';', expand=True)[1]
    rejected_raw['issue_flag'] = rejected_raw['FILTER'] + rejected_raw['INFO']
    keep_cols = ['hg38_chrom_pos_ref_alt', 'hg19_chrom_pos_ref_alt', 'rsid','issue_flag']
    return rejected_raw[keep_cols]

rule clean_results:
    input:
        rejected_vcf_path=f"output/{target_name}.rejected.vcf",
        lifted_vcf_path=f"output/{target_name}.hg38.vcf"

    output:
        lifted_text_path = f"output/{target_name}.hg38.txt"

    run:
        lifted_header_row = get_header_row(vcf_path=str(input.lifted_vcf_path))
        lifted_raw = pd.read_csv(str(input.lifted_vcf_path), sep='\t', skiprows=lifted_header_row, dtype={'POS' : 'str'})
        lifted_clean = clean_lifted(lifted_raw=lifted_raw)

        rejected_header_row = get_header_row(vcf_path=str(input.rejected_vcf_path))
        rejected_raw = pd.read_csv(str(input.rejected_vcf_path), sep='\t', skiprows=rejected_header_row, dtype={'POS' : 'str'})
        rejected_clean = clean_rejected(rejected_raw=rejected_raw)

        lifted_text = pd.concat([lifted_clean, rejected_clean])
        lifted_text.to_csv(str(output), sep='\t', index=None)




