# SSeqCpp
This is a cross platform C++ project to compute and maintain spectral sequences based on sqlite3 databases and algorithms on Groebner basis.

## AdamsE2
This is one of the projects that is more polished. Compile this project and in the `build/bin` folder run
```
./AdamsE2
```
It will compute a minimal resolution of the Steenrod algebra and save it in a database file `AdamsE2.db`. You can stop the program anytime and it can resume later from where it was stopped.

The computation of total degrees `t<=100` takes about 70 seconds with a single CPU thread on a PC. It can also run in parallel and utilize multiple CPU cores.
