# SSeqCpp
This is a cross platform C++ project to compute and maintain spectral sequences based on sqlite3 databases and algorithms on Groebner basis.

## Adams
This is the project that computes free resolutions over the Steenrod algebra. Compile this project and you can find the executable file `Adams` in the `build/bin` folder.
```bash
./Adams res S0 100
```
This computes a minimal resolution of F2 over the mod 2 Steenrod algebra up to total degree of 100. It outputs a database file `S0_Adams_res.db`.

You can stop the program anytime and it can resume later from where it was stopped.

This computation takes 7 seconds on my AMD Ryzen 9 3900X CPU.

```bash
./Adams prod 100
```
This computes some chain maps in order to obtain the products in the cohomology of the Steenrod algebra up to total degree of 100. `S0_Adams_res.db` is the input file and it outputs another database file `S0_Adams_res_prod.db`.

You can stop the program anytime and it can resume later from where it was stopped.

This computation takes 16 seconds on my AMD Ryzen 9 3900X CPU.

```bash
./Adams export 100
```
This generates the cohomology of the Steenrod algebra up to total degree of 100 according to the input file `S0_Adams_res_prod.db`. It outputs another database file `S0_AdamsSS.db`.

This computation takes 0.1 seconds on my AMD Ryzen 9 3900X CPU.

## Implementation limits of the program
- The maximum Adams filtration is $s_{max}=2^{12}-1=4095$.
- The maximum degree of the Steenrod algebra is $n_{max}=383$.
- The maximum dimension of a free module over the Steenrod algebra is $dim_{max}=2^{19}=524288$