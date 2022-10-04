import logging
# setting this to large will result in a deadlock of the program due to the thread waiting for an item in
# the queue which will never happen because the pipe used to feed data in the queue is full see also
# https://stackoverflow.com/questions/10028809/maximum-size-for-multiprocessing-queue-item/10029074
BUFFERSIZE = 200000

# maximum memory and average read package size in byte usage to compute maxsize of queue
# the size of one element in the readlist was computed with pympler.asizeof.asizeof((BC, B2, r1, r2))
# AVGSIZEOFQUEUEPACKAGE = 3368 * BUFFERSIZE
# MAXQUEUESIZE = (args.maxMemory*1000000000)//AVGSIZEOFQUEUEPACKAGE - ((args.writerThreads + 1) * 2)
def MAXQUEUESIZE(nproc: int) -> int:
    return nproc * 6


def set_log_level(verbose: bool):
    if verbose:
        logging.basicConfig(
            format = '%(asctime)s - %(message)s',
            level = logging.INFO
        )