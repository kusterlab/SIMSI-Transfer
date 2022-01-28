from __future__ import print_function

import sys
import signal
import warnings
import logging
from multiprocessing import Pool
import traceback

logger = logging.getLogger(__name__)


class JobPool:
  def __init__(self, processes = 1, warningFilter = "default"):
    self.warningFilter = warningFilter
    self.pool = Pool(processes, self.initWorker)
    self.results = []
    
  def applyAsync(self, f, args):
    r = self.pool.apply_async(f, args)
    self.results.append(r)
  
  def initWorker(self):
    return init_worker(self.warningFilter)
  
  def checkPool(self, printProgressEvery = -1):
    try:
      outputs = list()
      for res in self.results:
        outputs.append(res.get(timeout = 10000)) # 10000 seconds = ~3 hours
        if printProgressEvery > 0 and len(outputs) % printProgressEvery == 0:
          logger.info(f' {len(outputs)} / {len(self.results)} {"%.2f" % (float(len(outputs)) / len(self.results) * 100)}%')
      self.pool.close()
      self.pool.join()
      return outputs
    except (KeyboardInterrupt, SystemExit):
      logger.error("Caught KeyboardInterrupt, terminating workers")
      self.pool.terminate()
      self.pool.join()
      sys.exit()
    except Exception as e:
      logger.error("Caught Unknown exception, terminating workers")
      logger.error(traceback.print_exc())
      logger.error(e)
      self.pool.terminate()
      self.pool.join()
      sys.exit()


def init_worker(warningFilter):
  # set warningFilter for the child processes
  warnings.simplefilter(warningFilter)
  
  # causes child processes to ignore SIGINT signal and lets main process handle
  # interrupts instead (https://noswap.com/blog/python-multiprocessing-keyboardinterrupt)
  signal.signal(signal.SIGINT, signal.SIG_IGN)


def addOne(i):
  return i+1


def unitTest():
  pool = JobPool(4)
  for i in range(20):
    pool.applyAsync(addOne, [i])
  results = pool.checkPool()
  print(results)
