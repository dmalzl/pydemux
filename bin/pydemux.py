import pysam
import logging
import multiprocessing as mp
from pydemux.parser import main_arg_parser
from pydemux.core import (
    extract_barcodes,
    parallel_writer,
    sequential_writer,
    initialize
)
from pydemux.utils import *


def main():
    parser = main_arg_parser()
    args = parser.parse_args()

    logging.info('initializing...')
    stats_counter, out_handles, names, patterns, reg_exprs = initialize(
        args.barcodes,
        args.subparser_name,
        args.output_prefix,
        args.mismatches,
        args.gzip
    )

    if args.nprocesses > 1:
        write_lock = mp.Lock()
        write_queue = mp.Queue(
            maxsize = MAXQUEUESIZE(args.nprocessors)
        )
        result_queue = mp.Queue()

        write_processes = []
        for i in range(args.nprocesses):
            write_processes.append(
                mp.Process(
                    target = parallel_writer,
                    args = (
                        out_handles,
                        write_queue,
                        result_queue,
                        reg_exprs,
                        patterns,
                        write_lock,
                        i,
                        args.gzip,
                        args.verbose
                    )
                )
            )
            write_processes[-1].start()

    logging.info('demultiplexing...')

    with pysam.AlignmentFile(args.input, 'rb', check_sq = False) as infile:
        total_number_of_reads, reads_in_buffer = 0, 0
        read_list = []
        while True:
            try:
                r1 = infile.__next__()
                r2 = infile.__next__() if args.subparser_name == 'paired' else None

            except StopIteration:
                if args.nprocesses > 1:
                    write_queue.put(read_list)
                    write_queue.put([])

                else:
                    counter = sequential_writer(
                        out_handles,
                        reg_exprs,
                        patterns,
                        read_list,
                        args.gzip
                    )

                    for key, read_count in counter.items():
                        stats_counter[key] += read_count

                del read_list
                break

            barcodes = extract_barcodes(
                r1,
                r2,
                args.bc1tag,
                args.bc2tag if args.subparser_name == 'paired' else None
            )

            if args.subparser_name == 'single':
                reads = (r1.to_dict(),)
                barcodes = barcodes[:1]
            else:
                reads = (r1.to_dict(), r2.to_dict())

            for read in reads:
                read.pop('tags')

            read_list.append((barcodes, reads))
            total_number_of_reads += 1
            reads_in_buffer += 1

            if reads_in_buffer == BUFFERSIZE:
                if args.nprocesses > 1:
                    write_queue.put(read_list)

                else:
                    counter = sequential_writer(
                        out_handles,
                        reg_exprs,
                        patterns,
                        read_list,
                        args.gzip
                    )

                    for key, read_count in counter.items():
                        stats_counter[key] += read_count

                del read_list
                reads_in_buffer = 0
                read_list = []

            if total_number_of_reads % 1000000 == 0 and args.verbose:
                logging.info(
                    'processed %i read pairs' % total_number_of_reads
                )


    if args.nprocesses > 1:
        logging.info('waiting for writerprocess to finish')
        for process in write_processes:
            process.join()

        del write_queue, write_processes


        logging.info('retrieving demux statistics')
        result_queue.put({})
        while True:
            counter = result_queue.get()
            if counter:
                for key, read_count in counter.items():
                    stats_counter[key] += read_count

            else:
                break

    logging.info('writing statsfile...')
    with open(args.stats, 'w') as sfile:
        sfile.write('file\t#readpairs\n')

        for key, count in stats_counter.items():
            sfile.write('\t'.join([names[key], str(count)]) + '\n')

    for handles in out_handles.values():
        for handle in handles:
            handle.close()

    logging.info('done')


if __name__ == '__main__':
    main()
