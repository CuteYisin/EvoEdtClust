# EdtClust: A fast homologous protein sequences clustering method based on edit distance

`EdtClust` is a software suite to cluster huge protein sequence sets. `EdtClust` is open source GPL-licensed software implemented in C++ for Linux , MacOS, and Windows.

## Quick start

1. Installation can be made with the following command line

```shell
git clone https://github.com/YX-Xiang/EdtClust.git
```

2. Put your sequences in one FASTA file (in.fa) and replace the input and output filenames in the script file with the following command.

```shell
cd EvoEdtClust
vim run.sh
```

Finally just run the script `bash run.sh` to achieve clustering.

Call `EdtClust` without any parameters to see the usage (or use '-h' option):

```c++
Usage: <option(s)> input_fasta_file output
Options:
    -h,  Show this help message
    -t,  threshold of similarity (default: 0.5)
    -p,  number of threads used (default: 1)
```

The `output` argument may be : a path to a folder

## Citation

When referencing, please cite "Xiang Y, Gu J, Zhou J. EdtClust: A fast homologous protein sequences clustering method based on edit distance[C]//2023 IEEE International Conference on Bioinformatics and Biomedicine (BIBM). IEEE, 2023: 690-697."
[PubMed](https://ieeexplore.ieee.org/abstract/document/10385950/ "EdtClust@PubMed")

---
Please send bug reports to Yisin \(yisinx@mail.nankai.edu.cn\).<br />
EdtClust is freely available at
[https://github.com/YX-Xiang/EdtClust](https://github.com/YX-Xiang/EdtClust "EdtClust@Github").