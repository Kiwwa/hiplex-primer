import os

refdir              = ' --refdir ../exampledata/examplechromosome'
genes               = ' --genes ../exampledata/exampleregions.txt'
blocksize           = ' --blocksize 100'
maxprimersize       = ' --maxprimersize 30'
primervar           = ' --primervar 10'
splicebuffer        = ' --splicebuffer 8'
melt                = ' --melt 64'
log                 = ' --log ../exampledata/exampleoutput/output.log'
idtfile             = ' --idtfile ../exampledata/exampleoutput/output.csv'
roverfile           = ' --roverfile ../exampledata/exampleoutput/rover_input.tsv'
maxhairpinsize      = ' --maxhairpinsize 30'
blocksizevar        = ' --blocksizevar 5'
scale               = ' --scale 25nmole'
purification        = ' --purification HPLC'
senseheelseq        = ' --senseheelseq aaaa'
antisenseheelseq    = ' --antisenseheelseq gggg'

os.system('../primer_design/primer_design.py' + refdir + genes + blocksize + maxprimersize + 
									primervar + splicebuffer + melt + log + idtfile + 
									roverfile + blocksizevar + scale + purification + 
									senseheelseq + antisenseheelseq)