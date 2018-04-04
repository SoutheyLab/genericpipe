'''
Build the pipeline workflow by plumbing the stages together.
'''

from ruffus import Pipeline, suffix, formatter, add_inputs, output_from, regex
from stages import Stages
import glob
from utils import safe_make_dir

def make_pipeline_map(state):
    '''Build the pipeline by constructing stages and connecting them together'''
    # Build an empty pipeline
    pipeline = Pipeline(name='genericpipe')
    # Stages are dependent on the state
    stages = Stages(state)

    safe_make_dir('alignments')
    safe_make_dir('metrics')
    safe_make_dir('metrics/amplicon')
    safe_make_dir('metrics/summary')

    # The original FASTQ files
    fastq_files = glob.glob('fastqs/*')

    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.original_fastqs,
        name='original_fastqs',
        output=fastq_files)

    # Align paired end reads in FASTQ to the reference producing a BAM file
    pipeline.transform(
        task_func=stages.align_bwa,
        name='align_bwa',
        input=output_from('original_fastqs'),
        # Match the R1 (read 1) FASTQ file and grab the path and sample name.
        # This will be the first input to the stage.
        filter=formatter('fastqs/(?P<sample>[a-zA-Z0-9_-]+)_R1_(?P<lib>[a-zA-Z0-9-:]+).fastq.gz'),
        # Add one more inputs to the stage:
        #    1. The corresponding R2 FASTQ file
        add_inputs=add_inputs(
            'fastqs/{sample[0]}_R2_{lib[0]}.fastq.gz'),
        # Add an "extra" argument to the state (beyond the inputs and outputs)
        # which is the sample name. This is needed within the stage for finding out
        # sample specific configuration options
        extras=['{sample[0]}', '{lib[0]}'],
        # The output file name is the sample name with a .bam extension.
        output='alignments/{sample[0]}.sort.hq.bam')

    # generate mapping metrics.
    pipeline.transform(
        task_func=stages.generate_amplicon_metrics,
        name='generate_amplicon_metrics',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/amplicon/{sample[0]}.amplicon-metrics.txt',
        extras=['{sample[0]}'])

    pipeline.transform(
        task_func=stages.intersect_bed,
        name='intersect_bed',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.intersectbed.bam')

    pipeline.transform(
        task_func=stages.coverage_bed,
        name='coverage_bed',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.bedtools_hist_all.txt')

    pipeline.transform(
        task_func=stages.genome_reads,
        name='genome_reads',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.mapped_to_genome.txt')

    pipeline.transform(
        task_func=stages.target_reads,
        name='target_reads',
        input=output_from('intersect_bed'),
        filter=suffix('.intersectbed.bam'),
        output='.mapped_to_target.txt')

    pipeline.transform(
        task_func=stages.total_reads,
        name='total_reads',
        input=output_from('align_bwa'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9_-]+).sort.hq.bam'),
        output='metrics/summary/{sample[0]}.total_raw_reads.txt')

    pipeline.collate(
        task_func=stages.generate_stats,
        name='generate_stats',
        input=output_from('coverage_bed', 'genome_reads', 'target_reads', 'total_reads'),
        #filter=regex(r'.+/(.+BS\d{4,6}.+)\..+\.txt'),
        filter=regex(r'.+/(.+)\.(bedtools_hist_all|mapped_to_genome|mapped_to_target|total_raw_reads)\.txt'),
        output=r'metrics/summary/all_sample.summary.\1.txt',
        extras=[r'\1', 'all_sample.summary.txt'])
    
    summary_file = 'all_sample.summary.txt'	 

    (pipeline.originate(
        task_func=stages.grab_summary_file,
        name='grab_summary_file',
        output=summary_file)  
        .follows('generate_stats'))

    pipeline.transform(
        task_func=stages.filter_stats,
        name='filter_stats',
        input=output_from('grab_summary_file'),
        filter=suffix('.summary.txt'),
        output='.passed.summary.txt',
        extras=['all_sample.failed.summary.txt'])
 
    return pipeline
    
def make_pipeline_call(state):
    #this part of the pipeline will take the summary results of "map" and turn them into gatk and undr_rover vcfs
    pipeline = Pipeline(name='genericpipe')
    
    with open("all_sample.passed.summary.txt", 'r') as inputf:
        passed_files = inputf.read().split('\n')
    
        stages = Stages(state)

    safe_make_dir('variants')
    safe_make_dir('variants/gatk')
   
    pipeline.originate(
        task_func=stages.passed_filter_files,
        name='passed_filter_files',
        output=passed_files)

    ###### GATK VARIANT CALLING ######
    # Call variants using GATK
    pipeline.transform(
        task_func=stages.call_haplotypecaller_gatk,
        name='call_haplotypecaller_gatk',
        input=output_from('passed_filter_files'),
        filter=formatter('.+/(?P<sample>[a-zA-Z0-9-_]+).sort.hq.bam'),
        output='variants/gatk/{sample[0]}.g.vcf')
    
    return pipeline

def make_pipeline_process(state):
    # Define empty pipeline
    pipeline = Pipeline(name='genericpipe')
    # Get a list of paths to all the directories to be combined for variant calling
    run_directories = state.config.get_option('runs')
    #grab files from each of the processed directories in "runs"
    gatk_files = []
    
    for directory in run_directories:
        gatk_files.extend(glob.glob(directory + '/variants/gatk/*.g.vcf'))
    
    # Stages are dependent on the state
    stages = Stages(state)

    # This is a dummy stage. It is useful because it makes a node in the
    # pipeline graph, and gives the pipeline an obvious starting point.
    pipeline.originate(
        task_func=stages.glob_gatk,
        name='glob_gatk',
        output=gatk_files)
   
    safe_make_dir('variants')
    safe_make_dir('variants/gatk')

    # Combine G.VCF files for all samples using GATK
    pipeline.merge(
        task_func=stages.combine_gvcf_gatk,
        name='combine_gvcf_gatk',
        input=output_from('glob_gatk'),
        output='ALL.combined.vcf')

    # Genotype G.VCF files using GATK
    pipeline.transform(
        task_func=stages.genotype_gvcf_gatk,
        name='genotype_gvcf_gatk',
        input=output_from('combine_gvcf_gatk'),
        filter=suffix('.combined.vcf'),
        output='.raw.vcf')

    # Apply GT filters to genotyped vcf
    pipeline.transform(
        task_func=stages.genotype_filter_gatk,
        name='genotype_filter_gatk',
        input=output_from('genotype_gvcf_gatk'),
        filter=suffix('.raw.vcf'),
        output='.raw.gt-filter.vcf')

    # Decompose and normalise multiallelic sites
    pipeline.transform(
        task_func=stages.vt_decompose_normalise,
        name='vt_decompose_normalise',
        input=output_from('genotype_filter_gatk'),
        filter=suffix('.raw.gt-filter.vcf'),
        output='.raw.gt-filter.decomp.norm.vcf')

    # Annotate VCF file using GATK
    pipeline.transform(
        task_func=stages.variant_annotator_gatk,
        name='variant_annotator_gatk',
        input=output_from('vt_decompose_normalise'),
        filter=suffix('.raw.gt-filter.decomp.norm.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.vcf')

    # Filter vcf
    pipeline.transform(
        task_func=stages.gatk_filter,
        name='gatk_filter',
        input=output_from('variant_annotator_gatk'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vcf')

      
   #Apply VEP
    (pipeline.transform(
        task_func=stages.apply_vep,
        name='apply_vep',
        input=output_from('gatk_filter'),
        filter=suffix('.raw.gt-filter.decomp.norm.annotate.filter.vcf'),
        output='.raw.gt-filter.decomp.norm.annotate.filter.vep.vcf')
        .follows('index_final_vcf'))

    return pipeline
