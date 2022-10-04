import argparse as ap


def demux_arg_parser(subcommand: str) -> ap.ArgumentParser:
    parser = ap.ArgumentParser(add_help = False)
    default = parser.add_argument_group('default')

    if subcommand == 'single':
        default.add_argument(
            '--bctag',
            '-t',
            default = 'BC',
            help = 'BAM file tag carrying the barcode information'
        )

    if subcommand == 'paired':
        default.add_argument(
            '--bc1tag',
            '-t1',
            default = 'BC',
            help = 'BAM file tag carrying the barcode information for read 1'
        )
        default.add_argument(
            '--bc2tag',
            '-t2',
            default = 'B2',
            help = 'BAM file tag carrying the barcode information for read 2'
        )

    return parser


def generic_arg_parser() -> ap.ArgumentParser:
    parser = ap.ArgumentParser(add_help = False)
    required = parser.add_argument_group('Required arguments')
    required.add_argument(
        'input',
        metavar = 'input_bam',
        help = '''
            BAM formatted file containing either single or paired-end reads (the latter is assumed to be interleaved)
        '''
    )
    required.add_argument(
        '--barcodes',
        '-b',
        help = '''
            tab-separated file containing the name of the sample followed by the assigned barcode(s) for the reads 
            (i.e. single: barcode\tsample_name; paired: barcode1\tbarcode2\tsample_name)
        ''',
        required = True
    )
    required.add_argument(
        '--statistics',
        '-s',
        help = 'name of the file to write the demux statistics to',
        required = True
    )

    default = parser.add_argument_group('default')
    default.add_argument(
        '--output_prefix',
        '-o',
        default = '',
        help = 'prefix to use for the written output files. Can also contain directories.'
    )
    default.add_argument(
        '--processes',
        '-p',
        default = 1,
        type = int,
        help = 'number of concurrent processes to use for writing demultiplexed reads'
    )
    default.add_argument(
        '--mismatches',
        '-m',
        default = 0,
        type = int,
        help = 'number of allowed mismatches between the true and the sequenced barcode'
    )
    default.add_argument(
        '--gzip',
        '-gz',
        action = 'store_true',
        default = False,
        help = 'flag indicating if written fastq files should be gzipped. If set files are gzipped.'
    )
    default.add_argument(
        '--verbose',
        '-v',
        action = 'store_true',
        default = False,
        help = 'flag indicating if logs about the progress should be written to the console or not'
    )

    return parser


def main_arg_parser() -> ap.ArgumentParser:
    generic_args = generic_arg_parser()
    parser = ap.ArgumentParser(
        formatter_class = ap.RawDescriptionHelpFormatter,
        description = '''
        pydemux is an easy to use tool for demultiplexing sequencing data in BAM format
        for single- and paired-end reads. For detailed help on subcommands use either 
        
        pydemux single -h
        
        pydemux paired -h
        
    ''',
    epilog = '''
        example usages:
        pydemux single -s stats.tsv -b barcode.tsv single_end.bam
        pydemux paired -s stats.tsv -b barcode.tsv interleaved_paired_end.bam
    '''
    )
    subparsers = parser.add_subparsers(
        title = 'commands',
        dest = 'command',
        description = 'subcommands',
        help = 'subcommands',
        metavar = ''
    )

    subparsers.add_parser(
        'single',
        formatter_class = ap.ArgumentDefaultsHelpFormatter,
        parents = [
            generic_args,
            demux_arg_parser('single')
        ],
        help = '''
            demultiplex single end reads from BAM input. Use -t/--bctag to set the BAM tag in which the barcode is saved
            (default is BC). By default, no mismatches are allowed between true and sequenced barcode, you can use
            -m/--mismatches to change this. However, it is advised to set this argument to min(bc_hamming_dist) / 2 - 1
            in order to avoid erroneous assignment of reads to samples, where bc_hamming_dist is a set of pairwise
            hamming distances between all true barcodes. Reads that cannot be safely assigned according to this
            condition (i.e. number of mismatches between sequenced and true barcode is larger than --mismatches)
            are written to a separate ambiguous.fq file 
        ''',
        usage = '%(prog)s -s stats.tsv -b barcodes.tsv single_end.bam'
    )

    subparsers.add_parser(
        'paired',
        formatter_class = ap.ArgumentDefaultsHelpFormatter,
        parents = [
            generic_args,
            demux_arg_parser('paired')
        ],
        help = '''
            demultiplex paired end reads from BAM input. Use -t1/--bc1tag and -t2/--bc2tag to set the BAM tag 
            in which the barcodes of the first and second read are saved respectively (tags are searched for both reads)
            (default is -t1 BC - t2 B2). By default, no mismatches are allowed between true and sequenced barcode, you 
            can use -m/--mismatches to change this. However, it is advised to set this argument to 
            min(bc_hamming_dist) / 2 - 1 in order to avoid erroneous assignment of reads to samples, 
            where bc_hamming_dist is a set of pairwise hamming distances between all true barcodes. Reads that cannot 
            be safely assigned according to this condition (i.e. number of mismatches between sequenced and true 
            barcode is larger than --mismatches) are written to a separate ambiguous_[12].fq file 
        ''',
        usage = '%(prog)s -s stats.tsv -b barcodes.tsv single_end.bam'
    )

    return parser
