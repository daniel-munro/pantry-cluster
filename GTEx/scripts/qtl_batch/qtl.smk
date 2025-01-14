def groups_arg(wildcards, input):
    """Pass the phenotype groups file as a Python arg if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return f'--groups {input.groups}'
    else:
        return ''

def groups_input(wildcards):
    """Include the phenotype groups file as an input if applicable"""
    if modalities[wildcards.modality]['grouped']:
        return pheno_dir / f'{wildcards.modality}.phenotype_groups.txt'
    else:
        return []

rule tensorqtl_perm_chr:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
    output:
        interm_dir / 'qtl' / '{modality}.cis_qtl.{chrom}.txt.gz',
    params:
        geno_prefix = geno_prefix,
        qtl_dir = interm_dir / 'qtl',
        groups_arg = groups_arg,
    conda: 'tensorqtl'
    resources:
        walltime = 24,
        mem_mb = 32000,
        gpu = '-l ngpus=1',
    retries: 5
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            {params.groups_arg} \
            --mode cis \
            --chrom {wildcards.chrom} \
            >> log.{wildcards.modality}.cis.{wildcards.chrom}.txt 2>&1
        """

rule tensorqtl_perm:
    """Map cis-QTLs, determining significance using permutations.
    Outputs the top association per phenotype.
    """
    input:
        expand(interm_dir / 'qtl' / '{{modality}}.cis_qtl.{chrom}.txt.gz', chrom=[f'chr{x}' for x in range(1, 23)]),
    output:
        output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    params:
        prefix = str(interm_dir / 'qtl' / '{modality}.cis_qtl'),
        qtl_dir = output_dir / 'qtl',
    conda: 'pantry'
    shell:
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/assemble_tensorqtl.py \
            {params.prefix} \
            {output} \
            --mode cis
        """

rule tensorqtl_independent_chr_batch:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        geno = multiext(geno_prefix, '.bed', '.bim', '.fam'),
        bed = pheno_dir / '{modality}.bed.gz',
        bedi = pheno_dir / '{modality}.bed.gz.tbi',
        covar = interm_dir / 'covar' / '{modality}.covar.tsv',
        groups = groups_input,
        cis = output_dir / 'qtl' / '{modality}.cis_qtl.txt.gz',
    output:
        interm_dir / 'qtl' / '{modality}.cis_independent_qtl.{chrom}' / 'batch{batch}.txt.gz',
    params:
        batch_dir = str(interm_dir / 'qtl' / '{modality}.cis_independent_qtl.{chrom}'),
        geno_prefix = geno_prefix,
        groups_arg = groups_arg,
    conda: 'tensorqtl'
    resources:
        walltime = 48, # Give extra time for combined-phenotype mapping
        mem_mb = 32000,
        gpu = '-l ngpus=1',
    retries: 5
    shell:
        ## Cluster environments may require cuda to be loaded, e.g.:
        # module load cuda
        """
        mkdir -p {params.batch_dir}
        python3 scripts/run_tensorqtl.py \
            {params.geno_prefix} \
            {input.bed} \
            {output} \
            --covariates {input.covar} \
            --cis_output {input.cis} \
            {params.groups_arg} \
            --mode cis_independent \
            --chrom {wildcards.chrom} \
            --batch {wildcards.batch} \
            >> log.{wildcards.modality}.cis_independent.{wildcards.chrom}.{wildcards.batch}.txt 2>&1
        """

rule tensorqtl_independent_chr:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        expand(interm_dir / 'qtl' / '{{modality}}.cis_independent_qtl.{{chrom}}' / 'batch{batch}.txt.gz', batch=range(10)),
    output:
        interm_dir / 'qtl' / '{modality}.cis_independent_qtl.{chrom,[A-Za-z0-9]+}.txt.gz',
    params:
        uncompressed = str(interm_dir / 'qtl' / '{modality}.cis_independent_qtl.{chrom}.txt'),
    conda: 'pantry'
    shell:
        """
        zcat {input[0]} | head -n1 > {params.uncompressed}
        for f in {input}; do
            zcat $f | tail -n+2 >> {params.uncompressed}
        done
        gzip {params.uncompressed}
        """

rule tensorqtl_independent:
    """Use stepwise regression to identify multiple conditionally independent cis-QTLs per phenotype."""
    input:
        # lambda w: expand(interm_dir / 'qtl' / f'{w.modality}.cis_independent_qtl.{{chrom}}.txt.gz', chrom=[f'chr{x}' for x in range(1, 23)]),
        expand(interm_dir / 'qtl' / '{{modality}}.cis_independent_qtl.{chrom}.txt.gz', chrom=[f'chr{x}' for x in range(1, 23)]),
    output:
        output_dir / 'qtl' / '{modality}.cis_independent_qtl.txt.gz',
    params:
        prefix = str(interm_dir / 'qtl' / '{modality}.cis_independent_qtl'),
        qtl_dir = output_dir / 'qtl',
    conda: 'pantry'
    shell:
        """
        mkdir -p {params.qtl_dir}
        python3 scripts/assemble_tensorqtl.py \
            {params.prefix} \
            {output} \
            --mode cis_independent
        """

