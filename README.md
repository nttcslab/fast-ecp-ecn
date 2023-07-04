# A Fast and Exact Evaluation Algorithm for the Expected Number of Connected Nodes: an Enhanced Network Reliability Measure

This repository includes the codes to reproduce the results of experiments in the paper "A Fast and Exact Evaluation Algorithm for the Expected Number of Connected Nodes: an Enhanced Network Reliability Measure."

## Requirements

We use [TdZdd](https://github.com/kunisura/TdZdd) for implementing HH and HH+ methods. We also use [SAPPOROBDD](https://github.com/Shin-ichi-Minato/SAPPOROBDD) for implementing HH+ method. Before building our code, you must place the header files of [TdZdd](https://github.com/kunisura/TdZdd) into this directory. More specifically, all header files of [TdZdd/include/tdzdd](https://github.com/kunisura/TdZdd/tree/master/include/tdzdd) must be placed on `tdzdd/` directory; e.g. `tdzdd/DdEval.hpp`, `tdzdd/DdSpec.hpp`, etc. If you want to build HH+ method, you should build [SAPPOROBDD](https://github.com/Shin-ichi-Minato/SAPPOROBDD); after that, all header files must be placed on `SAPPOROBDD/include/` and the library file `BDD64.a` must be placed on `SAPPOROBDD/lib/`.

The codes for computing criticality (augmentability) of links also use [Adept](http://www.met.reading.ac.uk/clouds/adept/) C++ autodiff library. Before building them, you should install Adept with a higher level of optimization such as `-g -O3`. For the installation of Adept, please see [Adept's documentation](http://www.met.reading.ac.uk/clouds/adept/documentation.html).

After that, if your environment has CMake version >=3.8, you can build all codes with the following commands:

```shell
(moving to src/ directory)
mkdir release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

After running this, all binaries are generated in `release/` directory. Instead of `make`, you can build individual binary by the following commands:
```shell
make main          // Building the proposed method
make niny          // Building NINY method
make hh            // Building HH method; TdZdd is needed
make factoring     // Building factoring method
make hhplus        // Building HH+ method; TdZdd and SAPPOROBDD is needed
make crit-ecp      // Building the binary for computing ECP criticality; Adept is needed
make crit-allrel   // Building the binary for computing All-NR criticality; Adept is needed
```

### Verified environments

We verified that the building process of our codes and the commands presented below worked fine in the following macOS and Linux environments:

- macOS Big Sur 11.2.1 + Apple clang 12.0.0
- CentOS 7.9 + gcc 4.8.5

## Limitation

In our experiment, we include the algorithm for determining the link ordering by [the beam-search based heuristics of Inoue and Minato](https://www-alg.ist.hokudai.ac.jp/~thomas/TCSTR/TR/ABSTRACTS/abstr_tcstr16_80.html). This is due to the following reasons: (i) we included the elapsed time for computing the link ordering, and (ii) we computed the link ordering after performing some preprocessings such as removing degree 1 nodes. However, since the code for it is currently not publicly available, we excluded the codes for determining the link ordering. Instead, we applied this algorithm for the original graph (before performing preprocessings) and wrote the computed link ordering as the edgelist. This may produce similar experimental results, but not exactly the same.

## How to reproduce experimental results

All data used in our experiments are in `data.tar.gz`. After extracting, `*.txt` describes the graph files (as an edgelist), and `*.txt.prob` specifies each link's working probability. In addition, `src/*.txt.<n>.src` specifies the nodes in $T^\prime$ when $|T^\prime|=$`<n>`.

### Section VI-A

The proposed method can be executed by the following command:

```shell
./main [graph_file] [probability_file] [order_file] [T_file] <weight_file>
```

`[graph_file]`, `[probability_file]`, and `[T_file]` specify the path to the file describing the edgelist of the graph, each link's availability, and the list of the nodes in $T^\prime$, respectively. When $T^\prime=\emptyset$, `src/*.txt.0.src` should be specified for `[T_file]`. `[order_file]` specifies the path to the file describing the ordering among links. Note that we already computed the link order by the path-decomposition heuristics and wrote it on `*.txt`, so `[order_file]` can be specified with the same path as `[graph_file]` when reproducing the experimental results. `<weight_file>` is an optional augment specifying the path to the file describing the weight $w_i$ of nodes.

After execution, the program computes $S(\{v\})$ value for every node $v$.

_Running example:_ To compute every $S(\{v\})$ of Interoute topology, run:

```shell
./main ../data/0146-real-Interoute.edgelist.txt ../data/0146-real-Interoute.edgelist.txt.prob ../data/0186-real-UsCarrier.edgelist.txt ../data/src/0146-real-Interoute.edgelist.txt.0.src
```

The NINY and HH methods can also be executed by the following commands:

```shell
./niny [graph_file] [probability_file] [order_file] <weight_file>
./hh [graph_file] [probability_file] [order_file] <weight_file>
```

Note that since they only works when $T^\prime=\emptyset$, `[T_file]` is not needed.

### Section VI-B

Even if $|T^\prime|\neq\emptyset$, the proposed method can be executed by the same command:

```shell
./main [graph_file] [probability_file] [order_file] [T_file] <weight_file>
```

_Running example:_ To compute every $S(T^\prime\cup\{v\})$ of Interoute topology when $|T^\prime|=2$, run:

```shell
./main ../data/0146-real-Interoute.edgelist.txt ../data/0146-real-Interoute.edgelist.txt.prob ../data/0186-real-UsCarrier.edgelist.txt ../data/src/0146-real-Interoute.edgelist.txt.2.src
```

The factoring and HH+ methods can be executed by the following commands:

```shell
./niny [graph_file] [probability_file] [T_file] <weight_file>
./hh [graph_file] [probability_file] [order_file] [T_file] <weight_file>
```

### Section VI-C

The criticality measures with respect to ECP can be computed by the following command:

```shell
./crit-ecp [graph_file] [probability_file] [order_file] [mode] <weight_file>
```

`[graph_file]`, `[probability_file]`, and `[order_file]` are the same as above. The computed criticality measure can be switched by `[mode]`; ECP essentiality when `[mode]=1`, ECP augmentability when `[mode]=2`, and the contribution to ECP when `[mode]=3`.

_Running example:_ To compute ECP augmenrability for Interoute topology, run:

```shell
./crit-ecp ../data/0146-real-Interoute.edgelist.txt ../data/0146-real-Interoute.edgelist.txt.prob ../data/0186-real-UsCarrier.edgelist.txt 2
```

The criticality measures with respect to All-NR can also be computed by the following command:

```shell
./crit-allrel [graph_file] [probability_file] [order_file] [mode]
```

## License

This software is released under the NTT license, see `LICENSE.pdf`.
