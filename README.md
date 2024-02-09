Hello! Welcome to the Lionmap package, we're so excited to have you here!

This package is meant to enable easy, consistent and robust cell type mapping using a reference dataset. Optional features include using a custom hierarchy that describes the relationship between cell types, and setting different confidence thresholds for individual mappings.

There are a few different packages you have to install to use this package that are not hosted on the CRAN package directory. For most users, that will mean you have to manually install them (click the links below for instructions to install the following packages). You should also install the remotes library to install this package from Github.

Steps to have a successful install:

1)  Install the BPCells package (for more information [see here](https://bnprks.github.io/BPCells/))

2)  Install the [treeio package](https://bioconductor.org/packages/release/bioc/html/treeio.html)

3)  Install the Lionmap package: remotes::install_github("jonathan-columbiau/lionmap")

The package only works using BPCells matrices as input. BPCells is a package that enables fast gene expression transformations and IO operations while requiring little memory on your computing environment. For more info on BPCells, check out the package [here](https://bnprks.github.io/BPCells/).
