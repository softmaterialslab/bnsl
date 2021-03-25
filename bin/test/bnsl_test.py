import testbook
import os
from datetime import datetime
import subprocess

DEFAULT_E2K2RATIO = 1
DEFAULT_SALT = 0.175

@testbook.testbook('bnsl-frontend.ipynb', execute=True)
def test_runPreprocessor(tb):
    t0 = datetime.now()
    func = tb.ref("runPreprocessor")
    func(DEFAULT_E2K2RATIO, DEFAULT_SALT)
    logf = "./preprocessor.log"
    assert os.path.exists(logf) == 1
    t1 = datetime.fromtimestamp(os.path.getmtime(logf))
    assert (t1>t0)
    os.remove(logf)

'''
def test_backend():
    retcode1, retstr = subprocess.getstatusoutput("make run-preprocessor")
    assert (retcode1 == 0)
    retcode2, retstr = subprocess.getstatusoutput("time lmp < ./in.lammps")
    #assert (retcode2 == 0)
    print (retstr)
    mk1 = "Total # of neighbors"
    mk2 = "Ave neighs/atom"
    assert ((mk1 in retstr) and (mk2 in retstr))
'''
