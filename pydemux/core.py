import logging
import regex
import pysam
import gzip
import os
import multiprocessing as mp
from typing import IO, Union
from io import StringIO


class KeyMatcher:
    def __init__(
        self,
        regular_expressions: list[Union[tuple[regex], tuple[regex, regex]]],
        patterns: dict[regex, str]
    ):
        self.regex = regular_expressions
        self.patterns = patterns

    def barcode_to_key(self, barcodes: Union[tuple[str], tuple[str, str]]) -> str:
        key = None
        for reg_exprs in self.regex:
            if all_match(barcodes, reg_exprs):
                key = '+'.join(self.patterns[re] for re in reg_exprs)
                break

        return key


def extract_barcodes(
    read1: pysam.AlignedSegment,
    read2: Union[pysam.AlignedSegment, None],
    bc1tag: str,
    bc2tag: Union[str, None]
) -> tuple[str, Union[str, None]]:
    bc1 = read1.get_tag(bc1tag)
    bc2 = None

    if read2:
        try:
            bc2 = read2.get_tag(bc2tag)

        except KeyError:
            bc2 = read1.get_tag(bc2tag)

    return bc1, bc2


def bam_to_fq(
    read: dict,
    rn: int,
    bc: str
):
    '''
    takes a bamread as dictionary (i.e. read.to_dict) and returns a complementary fastq entry

    :param read:    bam read as dictionary as returned by read.to_dict()
    :param rn:      read number
    :param bc:      barcode sequence

    :return:        bam read in fastq format
    '''

    name = '@' + read['name'] + ' ' + ':'.join([str(rn), 'N', '0', bc])
    seq = read['seq']
    qual = read['qual']

    return '\n'.join([name, seq, '+', qual]) + '\n'


def all_match(
    barcodes: Union[tuple[str], tuple[str, str]],
    reg_exprs: list[Union[tuple[regex], tuple[regex, regex]]]
) -> bool:
    return all(re.match(bc) for re, bc in zip(reg_exprs, barcodes))


def split_dir_and_prefix(outfile_prefix: str) -> tuple[str, str]:
    tmp = outfile_prefix.split('/')
    return '/'.join(tmp[:-1]) + '/', tmp[-1]


def instantiate_objects(
    fields: list[str],
    outfile_prefix: str,
    subcommand: str,
    mismatches: int,
    gzipped: bool,
    mock: bool
) -> tuple[list[str], Union[tuple[str], tuple[str, str]], tuple[regex], str]:
    out_name = fields[-1]
    out_dir, out_prefix = split_dir_and_prefix(outfile_prefix)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    join_list = [out_prefix] if out_prefix else []
    suffix = '.gz' if gzipped else ''
    if subcommand == 'single':
        if mock:
            barcodes = ['ambiguous']
            handles = (
                out_dir + '_'.join(join_list + ['ambiguous']) + '.fq' + suffix,
            )
            res = None

        else:
            barcodes = fields[:1]
            handles = (
                out_dir + '_'.join(join_list + [out_name]) + '.fq' + suffix,
            )
            res = (
                regex.compile('(' + fields[0] + '){s<=' + str(mismatches) + '}'),
            )

    elif subcommand == 'paired':
        if mock:
            barcodes = ['ambiguous']
            handles = (
                out_dir + '_'.join(join_list + ['ambiguous', '1.fq' + suffix]),
                out_dir + '_'.join(join_list + ['ambiguous', '2.fq' + suffix])
            )
            res = None

        else:
            barcodes = fields[:2]
            handles = (
                out_dir + '_'.join(join_list + [out_name, '1.fq' + suffix]),
                out_dir + '_'.join(join_list + [out_name, '2.fq' + suffix])
            )
            res = tuple(
                [regex.compile('(' + b + '){s<=' + str(mismatches) + '}') for b in fields[:2]]
            )

    return barcodes, handles, res, out_name


def initialize(
    barcode_file: str,
    subcommand: str,
    outfile_prefix: str,
    mismatches: int,
    gzipped: bool = False
) -> tuple[dict[str, int], dict[str, tuple[IO]], dict[str, str], dict[regex, str], list[tuple[regex]]]:
    stat_counter, names = {}, {}
    reg_exprs, out_handles, patterns = [], {}, {}
    with open(barcode_file, 'r') as file:
        for line in file:
            fields = line.rstrip().split('\t')
            barcodes, handles, res, out_name = instantiate_objects(
                fields,
                outfile_prefix,
                subcommand,
                mismatches,
                gzipped,
                False
            )

            key = '+'.join(barcodes)
            out_handles[key] = handles
            names[key] = out_name
            stat_counter[key] = 0
            reg_exprs.append(res)

            for re, bc in zip(res, barcodes):
                patterns[re] = bc

    barcodes, handles, _, out_name = instantiate_objects(
        ['ambiguous'],
        outfile_prefix,
        subcommand,
        mismatches,
        gzipped,
        True
    )
    # there is actually no barcodes here just the ambiguous reads file
    key = '+'.join(barcodes)
    out_handles[key] = handles
    names[key] = key
    stat_counter[key] = 0

    return stat_counter, out_handles, names, patterns, reg_exprs


def parallel_writer(
    out_handles: dict[str, Union[tuple[str], tuple[str, str]]],
    write_queue: mp.Queue,
    result_queue: mp.Queue,
    reg_exprs: list[tuple[regex]],
    patterns: dict[regex, str],
    lock: mp.Lock,
    process_number: int,
    gzipped: bool,
    verbose: bool
):
    counter = {k: 0 for k in out_handles.keys()}
    file_chunks = {k: (StringIO(), StringIO()) for k in out_handles.keys()}
    key_matcher = KeyMatcher(
        reg_exprs,
        patterns
    )

    if verbose:
        logging.info(
            'writer-{0} waiting for input'.format(process_number)
        )

    while True:
        read_list = write_queue.get()
        if read_list:

            if verbose:
                logging.info(
                    'writer-{0} processing input'.format(process_number)
                )

            while len(read_list) > 0:
                barcodes, reads = read_list.pop(-1)
                key = key_matcher.barcode_to_key(barcodes)

                if key:
                    for stream, read, read_number in zip(
                        file_chunks[key],
                        reads,
                        [1, 2]
                    ):
                        stream.write(
                            bam_to_fq(
                                read,
                                read_number,
                                key
                            )
                        )
                    counter[key] += 1

                else:
                    for stream, read, read_number in zip(
                        file_chunks['ambiguous'],
                        reads,
                        [1, 2]
                    ):
                        stream.write(
                            bam_to_fq(
                                read,
                                read_number,
                                'ambiguous'
                            )
                        )
                    counter['ambiguous'] += 1

                del barcodes, reads

            del read_list

            lock.acquire()

            if verbose:
                logging.info(
                    'writer-{0} writing output'.format(process_number)
                )

            for key, streams in file_chunks.items():
                for filename, stream in zip(
                    out_handles[key],
                    streams
                ):
                    try:
                        if gzipped:
                            filehandle = open(filename, 'ab')
                            block = gzip.compress(
                                bytes(stream.getvalue(), 'utf-8'),
                                compresslevel = 9
                            )

                        else:
                            filehandle = open(filename, 'a')
                            block = stream.getvalue()

                        filehandle.write(block)
                        # needed to empty the file buffer which would
                        # otherwise result in not all reads being written
                        # to files
                        filehandle.flush()
                        filehandle.close()

                    except OSError as e:
                        logging.info(repr(e))

                file_chunks[key] = (StringIO(), StringIO())

            lock.release()

        else:
            if verbose:
                logging.info(
                    'writer-{0} has nothing left to process'.format(process_number)
                )

            write_queue.put([])
            break

    for streams in file_chunks.values():
        for stream in streams:
            stream.close()

    result_queue.put(counter)

    if verbose:
        logging.info(
            'nothing left to write. writer-{0} terminating'.format(process_number)
        )


def sequential_writer(
    out_handles: dict[str, Union[tuple[str], tuple[str, str]]],
    reg_exprs: list[tuple[regex]],
    patterns: dict[regex, str],
    read_list: list[tuple[tuple[str, str], tuple[dict[str, str], dict[str, str]]]],
    gzipped: bool
):
    counter = {k: 0 for k in out_handles.keys()}
    file_chunks = {k: (StringIO(), StringIO()) for k in out_handles.keys()}
    key_matcher = KeyMatcher(
        reg_exprs,
        patterns
    )

    while len(read_list) > 0:
        barcodes, reads = read_list.pop(-1)
        key = key_matcher.barcode_to_key(barcodes)

        if key:
            for stream, read, read_number in zip(
                file_chunks[key],
                reads,
                [1, 2]
            ):
                stream.write(
                    bam_to_fq(
                        read,
                        read_number,
                        key
                    )
                )
            counter[key] += 1

        else:
            for stream, read, read_number in zip(
                file_chunks['ambiguous'],
                reads,
                [1, 2]
            ):
                stream.write(
                    bam_to_fq(
                        read,
                        read_number,
                        'ambiguous'
                    )
                )
            counter['ambiguous'] += 1

    for key, streams in file_chunks.items():
        for filename, stream in zip(
            out_handles[key],
            streams
        ):
            try:
                if gzipped:
                    filehandle = open(filename, 'ab')
                    block = gzip.compress(
                        bytes(stream.getvalue(), 'utf-8'),
                        compresslevel = 9
                    )

                else:
                    filehandle = open(filename, 'a')
                    block = stream.getvalue()

                filehandle.write(block)
                # needed to empty the file buffer which would
                # otherwise result in not all reads being written
                # to files
                filehandle.flush()

            except OSError as e:
                logging.info(repr(e))

        file_chunks[key] = (StringIO(), StringIO())

    for streams in file_chunks.values():
        for stream in streams:
            stream.close()

    return counter
