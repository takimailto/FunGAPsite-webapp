import sys
from Bio import SeqIO
import os
import re
import shlex
from datetime import datetime
from subprocess import check_call


FunGAP_dir = '/home/mailto/FunGAP'  # according to user's Fungap installing
sys.path.append(FunGAP_dir)
from set_logging import set_logging

run_check_dependencies_path = os.path.join(FunGAP_dir, 'check_dependencies.py')
run_hisat2_path = os.path.join(FunGAP_dir, 'run_hisat2.py')
run_trinity_path = os.path.join(FunGAP_dir, 'run_trinity.py')
run_repeat_modeler_path = os.path.join(FunGAP_dir, 'run_repeat_modeler.py')

run_augustus_path = os.path.join(FunGAP_dir, 'run_augustus.py')
run_maker_path = os.path.join(FunGAP_dir, 'run_maker.py')
run_braker1_path = os.path.join(FunGAP_dir, 'run_braker1.py')

run_busco_path = os.path.join(FunGAP_dir, 'run_busco.py')
run_pfam_scan_path = os.path.join(FunGAP_dir, 'run_pfam_scan.py')
make_nr_prot_path = os.path.join(FunGAP_dir, 'make_nr_prot.py')
run_blastp_path = os.path.join(FunGAP_dir, 'run_blastp.py')
make_transcripts_path = os.path.join(FunGAP_dir, 'make_transcripts.py')
run_blastn_path = os.path.join(FunGAP_dir, 'run_blastn.py')
import_blast_path = os.path.join(FunGAP_dir, 'import_blastp.py')
import_busco_path = os.path.join(FunGAP_dir, 'import_busco.py')
import_pfam_path = os.path.join(FunGAP_dir, 'import_pfam.py')
import_blastn_path = os.path.join(FunGAP_dir, 'import_blastn.py')
catch_bad_genes_path = os.path.join(FunGAP_dir, 'catch_bad_genes.py')
filter_gff3s_path = os.path.join(FunGAP_dir, 'filter_gff3s.py')
gff3_postprocess_path = os.path.join(FunGAP_dir, 'gff3_postprocess.py')

copy_output_path = os.path.join(FunGAP_dir, 'copy_output.py')
create_markdown_path = os.path.join(FunGAP_dir, 'create_markdown.py')

download_sister_orgs_path = os.path.join(FunGAP_dir, 'download_sister_orgs.py')
get_augustus_species_path = os.path.join(FunGAP_dir, 'get_augustus_species.py')


def set_loggings(output_dir):
    create_dir(output_dir)
    log_file = os.path.join(output_dir, 'logs', 'fungap.log')
    global logger_time, logger_txt
    logger_time, logger_txt = set_logging(log_file)

    logger_txt.debug('\n============ New Run {} ============'.format(
        datetime.now())
    )


def create_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    log_dir = os.path.join(output_dir, 'logs')
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)


def run_hisat2(
        genome_assembly, trans_read_files, output_dir, num_cores, max_intron,
):
    if len(trans_read_files) == 1 and trans_read_files[0].endswith('.bam'):
        return trans_read_files

    hisat2_output_dir = os.path.join(output_dir, 'hisat2_out')
    log_dir = os.path.join(output_dir, 'logs')

    # run_hisat2.py -r <fastq1> <fastq2> <fastq3> ... \
    # -o <output_dir> -l <log_dir> -f <ref_fasta> -c <num_cores>
    # -m <max_intron>
    command = (
        'python {} --read_files {} --output_dir {} --log_dir {} --ref_fasta {} '
        '--num_cores {} --max_intron {}'.format(
            run_hisat2_path, ' '.join(trans_read_files), hisat2_output_dir,
            log_dir, genome_assembly, num_cores, max_intron
        ))
    logger_time.debug('START: wrapper_run_hisat2')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_hisat2\n')

    # Get output BAM file paths
    '''trans_bams = []
    for trans_read_file in trans_read_files:
        prefix = re.sub(r'_[12s]$', '',
                        os.path.basename(os.path.splitext(trans_read_file)[0])
                        )
        hisat2_output = os.path.join(hisat2_output_dir, '{}.bam'.format(prefix))
        trans_bams.append(hisat2_output)
    trans_bams2 = list(set(trans_bams))
    return trans_bams2
    '''


def run_trinity(
        trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron
):
    trinity_output_dir = os.path.join(output_dir, 'trinity_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_trinity.py -b <bam_files> -o <output_dir> -l <log_dir> -c <num_cores>
    # -m <max_intron> --jaccard_clip
    command = (
        'python {} --bam_files {} --output_dir {} --log_dir {} --num_cores {} '
        '--max_intron {} {}'.format(
            run_trinity_path, ' '.join(trans_bams), trinity_output_dir,
            log_dir, num_cores, max_intron,
            no_jaccard_clip
        )
    )
    logger_time.debug('START: wrapper_run_trinity')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_trinity\n')

    # Get output transcriptome assembly files
    '''trinity_asms = glob(os.path.join(
        output_dir, 'trinity_out', '*/Trinity_*.fasta')
    )
    return trinity_asms
    '''


def run_repeat_modeler(genome_assembly, output_dir, num_cores):
    # run_repeat_modeler.py -g <genome_assembly> -o <output_dir> -l <log_dir>
    # -c <num_cores>
    rm_output_dir = os.path.join(output_dir, 'repeat_modeler_out')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        'python {} --genome_assembly {} --output_dir {} --log_dir {} '
        '--num_cores {}'.format(
            run_repeat_modeler_path, genome_assembly, rm_output_dir, log_dir,
            num_cores
        )
    )
    logger_time.debug('START: wrapper_run_repeat_modeler')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_repeat_modeler\n')

    '''repeat_model_file = glob(
        os.path.join(rm_output_dir, 'RM*/consensi.fa.classified')
    )[0]
    return repeat_model_file
    '''


def run_maker(
        genome_assembly, output_dir, augustus_species, sister_proteome, num_cores,
        repeat_model_file, trinity_asms, no_genemark_fungus
):
    maker_out_dir = os.path.join(output_dir, 'maker_out')
    # run_maker.py -i <input_fasta> -a <augustus_species> -p <protein_db_fasta>
    # -R <repeat_model> -e <est_files> -o <output_dir> -c <num_cores>
    # -l <log_dir> --gmes_fungus
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        'python {} --input_fasta {} --augustus_species {} --protein_db_fasta {}'
        ' --repeat_model {} --est_files {} --output_dir {} --num_cores {} '
        '--log_dir {} {}'.format(
            run_maker_path, genome_assembly, augustus_species, sister_proteome,
            repeat_model_file, ' '.join(trinity_asms), maker_out_dir, num_cores,
            log_dir, no_genemark_fungus
        )
    )
    logger_time.debug('START: wrapper_run_maker')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_maker\n')

    '''maker_gff3s = glob(
        os.path.join(output_dir, 'maker_out', '*/maker_*.gff3')
    )
    maker_faas = glob(os.path.join(output_dir, 'maker_out', '*/maker_*.faa'))
    return maker_gff3s, maker_faas'''


def run_augustus(masked_assembly, output_dir, augustus_species):
    # run_augustus.py -m <masked_assembly> -s <species> -o <output_dir>
    # -l <log_dir>
    output_dir = os.path.join(output_dir, 'augustus_out')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        'python {} --masked_assembly {} --species {} --output_dir {} '
        '--log_dir {}'.format(
            run_augustus_path, masked_assembly, augustus_species, output_dir,
            log_dir
        )
    )
    logger_time.debug('START: wrapper_run_augustus')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_augustus\n')
    '''augustus_gff3 = os.path.join(output_dir, 'augustus.gff3')
    augustus_faa = os.path.join(output_dir, 'augustus.faa')
    return augustus_gff3, augustus_faa'''


def run_braker1(
        masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus
):
    braker1_output_dir = os.path.join(output_dir, 'braker1_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_braker1.py -m <masked_assembly> -b <bam_files> -o <output_dir>
    # -l <log_dir> -c <num_cores> --fungus
    command = (
        'python {} --masked_assembly {} --bam_files {} --output_dir {} '
        '--log_dir {} --num_cores {} {}'.format(
            run_braker1_path, masked_assembly, ' '.join(trans_bams),
            braker1_output_dir, log_dir, num_cores, no_braker_fungus
        )
    )
    logger_time.debug('START: wrapper_run_braker1')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_braker1\n')

    '''prefixes = [os.path.basename(os.path.splitext(x)[0]) for x in trans_bams]
    prefixes_u = list(set(prefixes))

    braker1_gff3s = []
    braker1_faas = []
    for prefix in prefixes_u:
        braker1_gff3 = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.gff3'.format(prefix)
        )
        braker1_gff3s.append(braker1_gff3)
        braker1_faa = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.faa'.format(prefix)
        )
        braker1_faas.append(braker1_faa)
    return braker1_gff3s, braker1_faas'''


def run_busco(input_faa, output_dir, num_cores):
    busco_output_dir = os.path.join(output_dir, 'busco_out')
    log_dir = os.path.join(output_dir, 'logs')
    # run_busco.py -i <input_fasta> -o <output_dir> -l <log_dir> -c <num_cores>
    command = (
        'python {} --input_fasta {} --output_dir {} --log_dir {} '
        '--num_cores {}'.format(
            run_busco_path, input_faa, busco_output_dir, log_dir, num_cores
        )
    )
    logger_time.debug('START: wrapper_run_busco')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_busco\n')


def run_busco_tatall(faa_files, output_dir, num_cores):
    for faa_file in faa_files:
        run_busco(faa_file, output_dir, num_cores)
    '''busco_out_dir = os.path.join(output_dir, 'busco_out')'''


def make_nr_prot(faa_files, output_dir,):
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    # make_nr_prot.py -i <faa_files> -o <output_dir>
    command = 'python {} --faa_files {} --output_dir {}'.format(
        make_nr_prot_path, ' '.join(faa_files), gene_filtering_dir
    )
    logger_time.debug('START: wrapper_make_nr_prot')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_make_nr_prot\n')

    '''nr_prot_file = os.path.join(gene_filtering_dir, 'nr_prot.faa')
    nr_prot_mapping_file = os.path.join(
        gene_filtering_dir, 'nr_prot_mapping.txt'
    )

    return nr_prot_file, nr_prot_mapping_file'''


def run_blastp(nr_prot_file, output_dir, sister_proteome, num_cores):
    # run_blastp.py -q <query_fasta> -d <db_fasta> -l <log_dir> -c <num_cores>
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        'python {} --query_fasta {} --db_fasta {} --log_dir {} '
        '--num_cores {}'.format(
            run_blastp_path, nr_prot_file, sister_proteome, log_dir,
            num_cores
        )
    )
    logger_time.debug('START: wrapper_run_blastp')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_blastp\n')

    '''blastp_output = os.path.join(output_dir, 'gene_filtering', 'nr_prot.blastp')

    return blastp_output'''


def run_pfam_scan(nr_prot_file, output_dir, num_cores):
    # run_pfam_scan.py -i <input_fasta> -l <log_dir> -c <num_cores>
    log_dir = os.path.join(output_dir, 'logs')
    command = 'python {} --input_fasta {} --log_dir {} --num_cores {}'.format(
        run_pfam_scan_path, nr_prot_file, log_dir, num_cores
    )
    logger_time.debug('START: wrapper_run_pfam_scan')
    logger_txt.debug('[Wapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_pfam_scan\n')

    '''pfam_scan_out = os.path.join(
        output_dir, 'gene_filtering', 'nr_prot.pfam_scan'
    )
    return pfam_scan_out'''


def concatenate_transcripts(trinity_asms, trinity_asm):
    set_loggings(os.path.join(trinity_asm, 'output'))
    command = 'cat {} > {}'.format(' '.join(trinity_asms), trinity_asm)
    logger_time.debug('Create transcript')
    logger_txt.debug('[Run] {}'.format(command))
    os.system(command)


def make_transcripts(genome_assembly, gff3_file):
    # make_transcripts.py -f <input_fasta> -g <input_gff3>
    command = 'python {} --input_fasta {} --input_gff3 {}'.format(
        make_transcripts_path, genome_assembly, gff3_file
    )
    logger_time.debug('START: wrapper_make_transcripts')
    logger_txt.debug('[Wapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_make_transcripts\n')

    gff3_base = os.path.splitext(gff3_file)[0]
    transcript_file = '{}_transcript.fna'.format(gff3_base)
    return transcript_file


def run_blastn(predicted_transcript, assembled_transcript, output_dir):
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    prefix = re.sub(
        r'_transcript\.fna', '', os.path.basename(predicted_transcript)
    )
    out_prefix = os.path.join(gene_filtering_dir, prefix)
    log_dir = os.path.join(output_dir, 'logs')
    # run_blastn.py -q <query_fasta> -d <db_fasta> -o <output_prefix>
    # -l <log_dir> -c <num_cores>
    command = (
        'python {} --query_fasta {} --db_fasta {} --output_prefix {} '
        '--log_dir {}'.format(
            run_blastn_path, predicted_transcript, assembled_transcript,
            out_prefix, log_dir
        )
    )
    logger_time.debug('START: wrapper_run_blastn')
    logger_txt.debug('[Wapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_blastn\n')

    blastn_out = '{}.blastn'.format(out_prefix)
    return blastn_out


def run_blastn_totall(gff3_files, genome_assembly, trinity_asm, output_dir):
    blastn_out_files = []
    for gff3_file in gff3_files:
        transcript_file = make_transcripts(genome_assembly, gff3_file)
        blastn_out_file = run_blastn(transcript_file, trinity_asm, output_dir)
        blastn_out_files.append(blastn_out_file)
    return blastn_out_files


def import_blastp(blastp_output, nr_prot_mapping_file):
    # import_blastp.py -b <blastp_out_file> -n <nr_prot_mapping>
    # blastp_out_dir = os.path.dirname(blastp_output)
    command = 'python {} --blastp_out_file {} --nr_prot_mapping {}'.format(
        import_blast_path, blastp_output, nr_prot_mapping_file
    )
    logger_time.debug('START: wrapper_import_blastp')
    logger_txt.debug('[Wrapper] {}'.format(command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_blastp\n')

    '''Return files
    blastp_dict = os.path.join(blastp_out_dir, 'blastp_score.p')

    return blastp_dict'''


def import_busco(busco_out_dir, output_dir):
    # import_busco.py -b <busco_dir> -o <output_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    command = 'python {} --busco_dir {} --output_dir {}'.format(
        import_busco_path, busco_out_dir, gene_filtering_dir
    )
    logger_time.debug('START: wrapper_import_busco')
    logger_txt.debug('[Wrapper] {}'.format(command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_busco\n')
    ''' # Return files
    busco_dict = os.path.join(gene_filtering_dir, 'busco_score.p')

    return busco_dict'''


def import_pfam(pfam_scan_out, nr_prot_mapping_file):
    # import_pfam.py -p <pfam_scan_out_file> -n <nr_prot_mapping>
    command = 'python {} --pfam_scan_out_file {} --nr_prot_mapping {}'.format(
        import_pfam_path, pfam_scan_out, nr_prot_mapping_file
    )
    logger_time.debug('START: wrapper_import_pfam')
    logger_txt.debug('[Wrapper] {}'.format(command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_pfam\n')

    # Return files
    '''pfam_scan_out_dir = os.path.dirname(pfam_scan_out)
    pfam_dict = os.path.join(pfam_scan_out_dir, 'pfam_score.p')
    return pfam_dict'''


def import_blastn(blastn_out_files, output_dir):
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    # import_blastn.py -b <blastn_out_files> -o <output_dir>
    command = 'python {} --blastn_out_files {} --output_dir {}'.format(
        import_blastn_path, ' '.join(blastn_out_files), gene_filtering_dir
    )
    logger_time.debug('START: wrapper_import_blastn')
    logger_txt.debug('[Wrapper] {}'.format(command))
    os.system(command)
    logger_time.debug('DONE : wrapper_import_blastn\n')
    '''blastn_dict = os.path.join(gene_filtering_dir, 'blastn_score.p')

    return blastn_dict'''


def catch_bad_genes(gff3_files, genome_assembly, output_dir):
    # catch_bad_genes.py -g <gff3_files> -a <genome_assembly> -o <output_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    command = (
        'python {} --gff3_files {} --genome_assembly {} --output_dir {}'.format(
            catch_bad_genes_path, ' '.join(gff3_files), genome_assembly,
            gene_filtering_dir
        )
    )
    logger_time.debug('START: wrapper_catch_bad_genes')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_catch_bad_genes\n')

    '''bad_dict = os.path.join(gene_filtering_dir, 'D_bad.p')
    return bad_dict'''


def filter_gff3s(
        genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict, blastn_dict,
        bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir
):
    # filter_gff3s.py -a <genome_assembly> -i <input_gff3s> -m <mapping_file>
    # -b <blastp_dict> -B <busco_dict> -p <pfam_dict> -N <blastn_dict> -g <bad_dict>
    # -n <nr_prot_file> -o <output_dir> -l <log_dir>
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    log_dir = os.path.join(output_dir, 'logs')
    command = (
        'python {} --genome_assembly {} --input_gff3s {} --mapping_file {} '
        '--blastp_dict {} --busco_dict {} --pfam_dict {} --blastn_dict {} '
        '--bad_dict {} --nr_prot_file {} --output_dir {} --log_dir {}'
    ).format(
        filter_gff3s_path, genome_assembly, ' '.join(gff3_files), nr_prot_mapping_file,
        blastp_dict, busco_dict, pfam_dict, blastn_dict, bad_dict, nr_prot_file,
        gene_filtering_dir, log_dir
    )
    logger_time.debug('START: wrapper_filter_gff3s')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_filter_gff3s\n')


def gff3_postprocess(genome_assembly, input_gff3):
    # gff3_postprocess.py -g <genome_assembly> -i <input_gff3> -o <output_gff3>
    output_gff3 = os.path.join(os.path.dirname(input_gff3), 'filtered_2.gff3')
    command = (
        'python {} --genome_assembly {} --input_gff3 {} --output_gff3 {}'
    ).format(
        gff3_postprocess_path, genome_assembly, input_gff3, output_gff3
    )
    logger_time.debug('START: wrapper_gff3_postprocess')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_gff3_postprocess\n')


def copy_output(output_dir):
    # copy_output.py -o <output_dir>
    command = 'python {} --output_dir {}'.format(copy_output_path, output_dir)
    logger_time.debug('START: wrapper_copy_output')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE: wrapper_copy_output\n')


def create_markdown(genome_assembly, output_dir, trans_bams, trinity_asms, fungap_gff3):
    # python create_markdown.py -f <input_fasta> -g <input_gff3>
    # -t <trinity_assembly> -b <bam_file> -o <output_dir>
    if fungap_gff3 == '':
        fungap_gff3 = os.path.join(output_dir, 'gene_filtering/filtered_2.gff3')
    trans_bam = trans_bams[0]
    trinity_asm = trinity_asms[0]
    markdown_out_dir = os.path.join(output_dir, 'fungap_out')

    command = (
        'python {} --input_fasta {} --input_gff3 {} --trinity_assembly {} '
        '--bam_file {} --output_dir {}'
    ).format(
        create_markdown_path, genome_assembly, fungap_gff3, trinity_asm,
        trans_bam, markdown_out_dir
    )
    logger_time.debug('START: wrapper_create_markdown')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE: wrapper_create_markdown\n')

    logger_time.debug('## DONE: FunGAP ##')


'''
Check if inputs are proper to run FunGAP
'''


def check_inputs(
        trans_read_1, trans_read_2, trans_read_single, trans_bam, genome_assembly,
        sister_proteome,
):
    print('Check input files...')

    if trans_read_1 and not os.path.exists(trans_read_1):
        return '[ERROR] No such file {}'.format(os.path.basename(trans_read_1))
    if trans_read_2 and not os.path.exists(trans_read_2):
        return '[ERROR] No such file {}'.format(os.path.basename(trans_read_2))
    if trans_read_single and not os.path.exists(trans_read_single):
        return '[ERROR] No such file {}'.format(os.path.basename(trans_read_single))
    trans_read_files = []
    if trans_read_1 and trans_read_2:
        # Check extension
        if (
                not trans_read_1.endswith('_1.fastq') and
                not trans_read_1.endswith('_1.fq')
        ):
            error_message = (
                '[ERROR] TRANS_READ_1 file name is incorrect. '
                'Should be <prefix>_1.fastq. You provided {}'.format(
                    trans_read_1
                )
            )
            return error_message

        elif (
                not trans_read_2.endswith('_2.fastq') and
                not trans_read_2.endswith('_2.fq')
        ):
            error_message2 = (
                '[ERROR] TRANS_READ_2 file name is incorrect. '
                'Should be <prefix>_2.fastq. You provided {}'.format(
                    trans_read_1)
            )
            return error_message2

        # Check prefix
        prefix_1 = os.path.basename(trans_read_1).replace('_1.fastq', '')
        prefix_2 = os.path.basename(trans_read_2).replace('_2.fastq', '')
        if prefix_1 != prefix_2:
            error_message3 = (
                '[ERROR] Two paired-end trans_read_files should have same'
                'prefix. <prefix>_1.fastq and <prefix>_2.fastq'
            )
            return error_message3

        trans_read_files = [trans_read_1, trans_read_2]

    elif trans_read_single:
        # Check extension
        if (
                not trans_read_single.endswith('_s.fastq') and
                not trans_read_single.endswith('_s.fq')
        ):
            error_message4 = (
                '[ERROR] TRANS_READ_SINGLE file name is incorrect. Should be '
                '<prefix>_s.fastq. You provided {}'.format(trans_read_single)
            )
            return error_message4
        trans_read_files = [trans_read_single]

    elif trans_bam:
        trans_read_files = [trans_bam]

    if not trans_read_files:
        error_message5 = (
            '[ERROR] You did not provide any transcriptome files: '
            '-1 and -2, -U, or -A should be provided'
        )
        return error_message5

    print('TRANS_READ_FILES is ok...')

    # Check geonme assembly in FASTA
    with open(genome_assembly, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            return '[ERROR] FASTA file is invalid: {}'.format(
                genome_assembly
            )

        for record in SeqIO.parse(handle, 'fasta'):
            if '|' in record.id:
                error_message6 = (
                    '[ERROR] FASTA defline contains "|" character, please '
                    'remove and re-run'
                )
                return error_message6

    print('GENOME_ASSEMBLY is ok...')

    # Check sister_proteome in FASTA
    with open(sister_proteome, 'r') as handle:
        fasta = SeqIO.parse(handle, 'fasta')
        if not any(fasta):
            return '[ERROR] FASTA file is invalid: {}'.format(
                genome_assembly
            )

    print('SISTER_PROTEOME is ok...')
    print('')

    return trans_read_files


def run_repeat_masker(genome_assembly, output_dir , repeat_model_file, num_cores):
    command = "FunGAP/external/RepeatMasker/RepeatMasker {} -dir{} -pa {} -lib {} -xsmall -e ncbi".format(
        genome_assembly, output_dir, num_cores, repeat_model_file)
    logger_time.debug('START: wrapper_run_repeat_masker')
    logger_txt.debug('[Wrapper] {}'.format(command))
    command_args = shlex.split(command)
    check_call(command_args)
    logger_time.debug('DONE : wrapper_run_repeat_masker\n')


def download_sister_orgs(taxon, email_address, num_sisters, output_dir):
    command = "python {} \
  --download_dir {} \
  --taxon '{}' \
  --num_sisters {} \
  --email_address {}".format(download_sister_orgs_path, output_dir, taxon, num_sisters, email_address)
    command_args = shlex.split(command)
    check_call(command_args)
    os.system("zcat {}/*faa.gz > {}/prot_db.faa".format(output_dir, output_dir))


def get_augustus_specie(genus_name, email_address):
    command = "python {} \
  --genus_name '{}' \
  --email_address {}".format(get_augustus_species_path, genus_name, email_address)

