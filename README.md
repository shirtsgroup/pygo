pygo
====

Python implementation of mixed MC/MD for a Go-like model for proteins

Requirements:
* [Parallel Python](http://www.parallelpython.com)
* numpy
* matplotlib

````
-bash-4.1$ python simulateGO.py -h

Usage: simulateGO.py [options]

Options:
  -h, --help            show this help message and exit
  -f FILENAME, --files=FILENAME
                        protein .pdb file
  -p PARAMFILE, --paramfile=PARAMFILE
                        protein .param file
  -t TEMPRANGE, --temprange=TEMPRANGE
                        temperature range of replicas
  --tfile=TFILE         file of temperatures
  -r NREPLICAS, --nreplicas=NREPLICAS
                        number of replicas
  -n TOTMOVES, --moves=TOTMOVES
                        total number of moves
  -s SAVE, --save=SAVE  number of moves between save operations
  -k SWAP, --swapstep=SWAP
                        number of moves between swap moves
  -w, --writetraj       flag to write out trajectory (default: False)
  --swapsper=SWAPSPER   number of swaps to perform (default: 500)
  --id=ID               the simlog id number or umbrella id number
  --freq=FREQ           ratio of move frequencies (tr:ro:an:di:gc:pr:md)
  --md=MD               step size (fs) and nsteps for MD move
  --surf                surface simulation flag
  --surfparamfile=SURFPARAMFILE
                        surface param file
  --scale=SCALE         scaling of surface attraction strength
  --umbrella=UMBRELLA   umbrella simulation flag and distance of pinning
  --umbrellak=UMBRELLAK
                        umbrella spring constant (default: 1)
  -Q QFILE, --Qumbrella=QFILE
                        Q umbrella simulation flag and file of Q_pins
  --k_Qpin=K_QPIN       Q umbrella spring constant (default: 10)
  --cluster             flag for running on cluster
  --restart             restart from a checkpoint
  --extend=EXTEND       id number of existing simulation to start simulation
                        from
````
See examples subdirectory for example usages.
