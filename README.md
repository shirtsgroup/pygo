pygo
====

Python implementation of mixed MC/MD for a Go-like model for proteins

Requirements:
* [Parallel Python](http://www.parallelpython.com)
* numpy
* matplotlib

````
$ simulateGO.py -h

usage: simulateGO.py [-h] [-f FILENAME] [-p PARAMFILE]
                     [-t TEMPRANGE TEMPRANGE] [--tfile TFILE] [-r NREPLICAS]
                     [-n TOTMOVES] [-s SAVE] [-k SWAP] [--nswap NSWAP] [-w]
                     [--id ID] [--freq x x x x x x x] [--md MD MD] [-o ODIR]
                     [--surf] [--surfparamfile SURFPARAMFILE] [--scale SCALE]
                     [-Z ZUMBRELLA] [--k_Zpin K_ZPIN] [-Q QFILE]
                     [--k_Qpin K_QPIN] [--cluster] [--restart] [--extend ID]

Run a simulation

optional arguments:
  -h, --help            show this help message and exit

Input files and parameters:
  -f FILENAME           protein GO_xxxx.pdb file
  -p PARAMFILE          protein GO_xxxx.param file
  -t TEMPRANGE TEMPRANGE
                        temperature range of replicas (default: 200, 400)
  --tfile TFILE         file of temperatures
  -r NREPLICAS, --nreplicas NREPLICAS
                        number of replicas (default: 8)
  -n TOTMOVES, --nmoves TOTMOVES
                        total number of moves (default: 10000)
  -s SAVE, --save SAVE  save interval in number of moves (default: 1000)
  -k SWAP, --swap SWAP  replica exchange interval in number of moves (default:
                        1000)
  --nswap NSWAP         number of attempted exchanges at each (default: 500)
  -w, --writetraj       flag to write out trajectory (default: False)
  --id ID               the simlog id number or umbrella id number (default:
                        0)
  --freq x x x x x x x  ratio of move frequencies (tr:ro:an:di:gc:pr:md)
                        (default: 0:0:1:3:3:3:10)
  --md MD MD            step size (fs) and number of steps for MD move
                        (default: 45 fs, 50 steps)
  -o ODIR, --odir ODIR  output directory (default: ./)

Surface simulation input files and parameters:
  --surf                surface simulation flag (default: False)
  --surfparamfile SURFPARAMFILE
                        surface param file (default: avgsurfparam.npy)
  --scale SCALE         scaling of surface attraction strength (default: 1)

Umbrella simulation input files and parameters:
  -Z ZUMBRELLA, --Zumbrella ZUMBRELLA
                        umbrella simulation flag and distance of z pin
                        (default=0)
  --k_Zpin K_ZPIN       Z umbrella spring constant (default: 1)
  -Q QFILE, --Qumbrella QFILE
                        Q umbrella simulation flag and file of Q_pins
  --k_Qpin K_QPIN       Q umbrella spring constant (default: 10)

Other specifications:
  --cluster             flag for running on cluster
  --restart             restart from checkpoint files
  --extend ID           id number of existing simulation to extend
````

See examples subdirectory for example usages.
