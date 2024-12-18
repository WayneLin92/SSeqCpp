# SSeqCpp
This is a cross platform project to compute and manage databases of spectral sequences based. The main programming language is C++. Recently we used this project to solve [the Last Kervaire Invariant Problem](https://arxiv.org/abs/2412.10879). See [Zenodo](https://zenodo.org/record/14475507) for more detail on the programming part.

## Adams
This project computes E2-page data of Adams spectral sequences at the prime 2. Compile this project and you can find an executable file `Adams` in the `build/bin` folder.

Features: This program is very robust in design. The output of the program will remain valid if
* there is a power failure during the computation,
* the program is terminated by the user or the system,
* two instances of the program are trying to write to the same file at the same time.

If the program is terminated for some reason, you can always resume the calculation from a recent checkpoint by running the same command again.

### Ring spectra
Assume that $R$ is a ring spectrum. To compute a minimal resolution of $H^*(R)$ over the Steenrod algebra up to internal degree `t=100`:
```bash
./Adams res R 100
```
It outputs a database file `R_Adams_res.db`.

To compute the multiplicative structure of the E2-page:
```bash
./Adams prod R 100
```
It outputs a database file `R_Adams_res_prod.db`.

To export the E2-page (the original resolution data is big):
```bash
./Adams export R 100
```
It outputs a database file `R_AdamsSS.db`.

### Modules over ring spectra
Assume that $M$ is a module over a ring spectrum $R$. To compute a minimal resolution of $H^*(M)$ over the Steenrod algebra up to internal degree `t=100`:
```bash
./Adams res M 100
```
It outputs a database file `M_Adams_res.db`.

To compute the module structure of the E2-page of $M$ over the E2-page of $R$:
```bash
./Adams prod_mod M R 100
```
It outputs a database file `M_Adams_res_prod.db`.

To export the E2-page:
```bash
./Adams export_mod M R 100
```
It depends on both `M_Adams_res_prod.db` and `R_AdamsSS.db`. It outputs a database file `M_AdamsSS.db`.

### d2 differentials
Assume that $X$ is a CW spectrum.
To compute the $d_2$ differentials:
```bash
./Adams d2 X d2
```
It outputs a database file `X_Adams_d2.db`.

To export the $d_2$ differentials:
```bash
./Adams export X 100
```
It outputs the results in a new column `d2` in the table `X_AdamsE2_basis` in the database file `X_AdamsSS.db`.

### Maps
Assume that there is a map between CW spectra named by `X__Y`. To compute this map between the E2-pages up to internal degree `t=100`:
```bash
./Adams map_res X Y 100
```
It outputs a database file `map_Adams_res_X__Y.db`.

To export the map on E2-pages:
```bash
./Adams export_mod X Y 100
```
It outputs a database file `map_AdamsSS_X__Y.db`.

### Custom definitions
The program has built-in support for spectra `S0`, `tmf` and stunted projective spaces `RPm_n`, `CPm_n`, `HPm_n`. The user can define their own spectra and maps in a file `Adams.json`. It should be placed in the same folder as the executable file. The format of the file is as follows:
```json
{
  "CW_complexes": {
    "C2": {
      "cells": [0, 1],
      "cells_gen": [0],
      "operations": [
        [0, 1]
      ]
    },
    "C2_C2": {
      "cells": [0, 1, 1, 2],
      "cells_gen": [0, 1],
      "operations": [
        [0, [1, 1]], 
        [0, 2], 
        [[1, 0], 2]
      ]
    },
    "CW_sigma_2sigma": {
      "cells_gen": [0],
      "cells": [0, 8, 16],
      "operations": [
        [0, 8],
        [0, 16]
      ]
    }
  },
  "maps": {
    "C2__C2_C2": {
      "from": "C2",
      "to": "C2_C2",
      "images": [0, []],
      "sus": 0
    }
  }
}
```

For CW spectra, the value for `cells` should be the ordered list of dimensions of cells. `cells_gen` should be a minimal list of generators of the cohomology as a module over the Steenrod algebra. `operations` should be a list of all nontrivial $Sq^{2^k}$ operations on the cells. If there is only one cell in dimension $d$, then we use `d` to denote the cell. If there are more than one cell in dimension $d$, we use `[d, i]` to denote the $i$'th cell in dimension $d$.

For a map $X\to \sum^k Y$, the value for `images` should be a list of images of the generators of $H^*(Y)$ in $H^*(X)$. The value for `sus` should be the number $k$.

### Implementation limits of the program
- The maximum Adams filtration is $s_{max}=2^{12}-1=4095$.
- The maximum degree of the Steenrod algebra is $n_{max}=383$.
- The maximum dimension of a free module over the Steenrod algebra is $dim_{max}=2^{19}=524288$

## ss
This project propagates Adams differentials and extensions among CW spectra. Compile this project and you can find an executable file `ss` in the `build/bin` folder.

Features:
* The program can deduce Adams differentials and extensions and produced proofs for the results.
* If the program is terminated by the user or the system or due to a power failure, the database files including the logging of proofs will remain the state before this run. This feature is due to the design of the modern database system [sqlite3](https://www.sqlite.org/).
* The program can generate plots for Adams differentials and extensions along cofiber sequences. The plots are actually webpages written in javascript. These webpages can generate latex tikz code for the plots if the user take a screenshot via the webpage.

After the user calculates a collection of Adams $E_2$ data using the `Adams` program, the user should place all the database files (containing `AdamsSS` in the name) in a folder next to the executable file `ss`. The folder name should not contain any space. We call this folder a *category* if it contains a configuration file `ss.json` describing the objects and maps that `ss` should consider.

A configuration file `ss.json` inside a category looks like the following.
```json
{
  "rings": [
    { "name": "S0", "path": "S0_AdamsSS.db", "deduce": "on", "plot_factors": [[[0, 1, 0], [1, 1, 0], [3, 1, 0]], [[7, 1, 0], [8, 3, 0], [9, 5, 0], [11, 5, 0], [14, 4, 0], [20, 4, 0]]] }
  ],
  "modules": [
    { "name": "C2", "path": "C2_AdamsSS.db", "over": "S0", "deduce": "on" },
    { "name": "Ceta", "path": "Ceta_AdamsSS.db", "over": "S0", "deduce": "on" }
  ],
  "maps": [
    { "name": "C2__S0", "path": "map_AdamsSS_C2__S0.db", "from": "C2", "to": "S0", "sus": 1, "t_max": 200 },
    { "name": "Ceta__S0", "path": "map_AdamsSS_Ceta__S0.db", "from": "Ceta", "to": "S0", "sus": 2, "t_max": 200 }
  ],
  "maps_v2": [
    { "name": "S0__S0_by_2", "factor": [0, 1, 0], "from": "S0", "to": "S0", "t_max": 260 },
    { "name": "S0__S0_by_eta", "factor": [1, 1, 0], "from": "S0", "to": "S0", "t_max": 259 },
    { "name": "S0__C2", "factor": [0, 0, 0], "from": "S0", "to": "C2", "t_max": 200 },
    { "name": "S0__Ceta", "factor": [0, 0, 0], "from": "S0", "to": "Ceta", "t_max": 200 }
  ],
  "cofseqs": [
    { "name": "S0__C2__S0", "path": "cofseq_S0__C2__S0.db", "i": "S0__C2", "q": "C2__S0", "d": "S0__S0_by_2" },
    { "name": "S0__Ceta__S0", "path": "cofseq_S0__Ceta__S0.db", "i": "S0__Ceta", "q": "Ceta__S0", "d": "S0__S0_by_eta" }
  ],
  "dir_plot": "mix"
}
```

We explain the fields in the configuration file. 
* `rings` is a list of ring spectra. 
* `modules` is a list of modules over the ring spectra.
* `maps` is a list of maps between spectra calculated by `./Adams`.
* `maps_v2` is a list of maps that is implied by the multiplicative structures of $E_2$-pages. `factor`: [stem,s,i] is the `i`'th basis element in degree (stem,s). We need these maps for defining cofiber sequences.
* `cofseqs` is a list of cofiber sequences. Usually `i` should be an inclusion of CW complexes. `q` should be a quotient map of CW complex. `d` is the connecting map.
* `dir_plot` is the folder name for the output plots. `plot_factors` in rings is the list of multiplications that we want the plots to support.

After the category in folder `category_name` is defined we should initialize the category by running the following command.
```bash
./ss reset category_name
```

To manually add an Adams differential:
```bash
./ss add_diff cw_name stem s r x dx category_name
```
For example, the following is $d_2h_4=h_0h_3^2$ in the Adams spectral sequence of `S0`.
```bash
./ss add_diff S0 15 1 2 0 0 category_name
```

To automatically compute Adams differentials and extensions in the category:
```bash
./ss deduce auto category_name num=3
```
Here `num` is the number of iterations. Each iteration will scan everywhere to try to deduce differentials.

By raising some flags, we can push the program harder in finding contradictions when assumptions are being made. There are several useful flags: `xy`, `syn`, `zero`, `pullback`. The program has a feature that can remember contradictions that happened before. Therefore it is efficient to apply zero or one flag at a time. To raise a flag (`zero` for example), run the following command.
```bash
./ss deduce auto category_name flags=zero num=3
```
The program will make records of human readable proofs of these results in the database file `log.db` in the category folder.

Before making plots for Adams differentials and extensions, the user should place a configuration file `ss.json` along side the executable file `ss`. You can think of this file as global configurations while each `ss.json` inside a category is a setting only for that category. The user should specify a path to hold webpages for the plots in `ss.json`. The webpage files can be found on [Zenodo](https://zenodo.org/record/14475507). For example, assuming that the webpage files are in a folder named `webpages`, the content of `ss.json` should be
```json
{
  "dir_website_ss": "webpages"
}
```

After this setup, the user can generate plots for the Adams differentials by running the following command.
```bash
./ss plot_ss category_name
```
To generate plots for the extensions along cofiber sequences:
```bash
./ss plot_cofseq category_name
```
To view the plots, open the file `plot.html` in the folder `dir_website_ss` by a browser. You can add parameters to the url to control the display of the plots. For example, `plot.html?data=Ceta` will show the plots for the Adams spectral sequence of `Ceta`, while `plot.html?data=Ceta__S0` will show the map $C\eta\to S^2$.

The reader can find plots we generate for [the Last Kervaire Invariant Problem](https://arxiv.org/abs/2412.10879) on this [webpage](https://waynelin92.github.io/ss/kervaire.html).

## Adamsp
This project aims to compute the Adams spectral sequences at odd primes $p$. It adapts some algorithms from the `Adams` project. The main contributor to this project is [Yu Zhang](https://github.com/yuzhangmath).
