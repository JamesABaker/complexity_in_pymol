
import subprocess
from subprocess import Popen,PIPE
results = subprocess.Popen(["perl", "TMSOC.pl", "sequence.fasta", "segments.txt"], stdout=PIPE, stderr=PIPE)
print results.stdout.read()
