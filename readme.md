
# IBDkin

**IBDkin** is a software for **IBD**-based **kin**ship estimation. 
IBDkin scales to hundreds of billions of IBD segments detected in hundreds of thousand individuals.

If you use this software in your published analysis, please cite:

> Ying Zhou, Sharon R Browning, Brian L Browning, IBDkin: fast estimation of kinship coefficients from identity by descent segments, Bioinformatics, btaa569, https://doi.org/10.1093/bioinformatics/btaa569
Last update: June 18, 2020,
by Ying Zhou, yz001(at)uw(dot)edu



Content
--------

- [1 Installation](#1-installation)
- [2 Running IBDkin](#2-running-ibdkin)
    - [2.1 Example analysis](#21-example-analysis)
    - [2.2 Required Parameters](#22-required-parameters)
    - [2.3 Optional Parameters](#23-optional-parameters)
    - [2.4 Flags](#24-flags)
- [3 Output Files](#3-output-files)
    - [3.1 Kinship Coefficients](#31-kinship-coefficients)
    - [3.2 IBD Coverage](#32-ibd-coverage)
    - [3.3 Masked Regions](#33-masked-regions)
- [4 License](#4-license)



# 1 Installation

The following commands download the source code, change the working directory to the source code folder "IBDkin/src-v2.8.7.7/", and create the executable file, `IBDkin`:

```bash

git clone https://github.com/YingZhou001/IBDkin.git

cd IBDkin/src-v2.8.7.7/
make
```
**IBDkin** is compiled under linux CentOS 7.5 using gcc 4.8.5. If you encounter any problems compiling the program, please contact the author for assistance. 


[\[top\]](#ibdkin)

# 2 Running IBDkin

To run some example IBDkin analyses, change the working directory to "IBDkin/example.pub/run", copy the executable file "IBDkin" into this folder, and enter `sh run.sh` to run the examples in sections 2.1.1 to 2.1.3.

```bash
cd IBDkin/example.pub/run
cp ../../src-v2.8.7.7/IBDkin ./
sh run.sh
```
Alternatively, you can enter the commands in sections 2.1.1 to 2.1.3 to run the examples.

[\[top\]](#ibdkin)

## 2.1 Example Analysis

We first change the working directory to the [IBDkin/example.pub/run](example.pub/run) folder with the command:
```bash
cd IBDkin/example.pub/run
```
This folder contains four required **IBDkin** input files: "ibd.txt" is the list of files containing IBD segments, "ind.txt" is the list of individuals, "range.txt" is the ranges of markers, and "plink.map" is the genetic map (see [Required Parameters](#22-required-parameters)).
Before running the following scripts, we need to copy the executable file `IBDkin` to the working directory with command:

```bash
cp ../../src-v2.8.7.7/IBDkin ./
```

### 2.1.1 Computing kinship coefficients
We first run **IBDkin** with five threads (`--nthreads 5`), 
and output kinship coefficients between pairs of individuals who share at least one 4cM IBD segment (see the '--cutcm' parameter in the [Optional Parameters](#23-optional-parameters) section). 
The output file prefix is specified with the `-out` flag.
The output kinship coefficients are written to the gzip-compressed file "example-1.kinship.gz":

```bash
./IBDkin --ibdfile ibd.txt --map plink.map --ind ind.txt --range range.txt --nthreads 5 --out example-1
```

### 2.1.2 Computing IBD coverage and masked regions
In certain situations, we may want to output the masked regions (`--outmask` flag) and the IBD coverage across the genome (`--outcoverage` flag). 
The output masked regions will be written to the gzip-compressed file "example-2.mask.gz", and the output IBD coverage will be written to the gzip-compressed file "example-2.coverage.gz". 

```bash
./IBDkin --ibdfile ibd.txt --map plink.map --ind ind.txt --range range.txt --nthreads 5 --out example-2 --outmask --outcoverage
```

### 2.1.3 Distributed analysis
If there is not enough memory to run an analysis on a huge data set, we can use the option `--part` to distribute computation across multiple compute nodes.
In the following script, we divide the analysis into five parts and run each part separately by assigning integers, from 1 to 5, to the variable `${part}`.


```bash
part=1 # can be set as 2, 3, 4, and 5
./IBDkin --ibdfile ibd.txt --map plink.map --ind ind.txt --range range.txt --nthreads 5 --out example-3.${part} --part 5 ${part}
```
**Caution!** Specify a different output file prefix for each part when use the `--part` option, or you will overwrite your output files.

For advanced parameters setting, please refer to our [manuscript(add link to our manuscript)](#ibdkin) and check the following sections: [Required Parameters](#22-required-parameters) and [Optional Parameters](#23-optional-parameters) sections.

[\[top\]](#ibdkin)


## 2.2 Required Parameters

* **--ibdfile [file]** #\<string\> the [file] contains the pathnames of files that list the IBD segments on each chromosome (one pathname per line, and one line per chromosome). The IBD segments must be stored in gzip-compressed [hap-IBD format](https://github.com/browning-lab/hap-ibd).
* **--map [file]** #\<string\> the [file] is a genetic map with cM distances in [PLINK format](http://zzz.bwh.harvard.edu/plink/data.shtml#map), including all chromosomes having IBD segments. The chromosome identifiers must be consistent in the IBD inputs, the genetic map, and the range file.
* **--ind [file]** #\<string\> the [file] includes a list of individuals to be analyzed (one individual per line).
* **--range [file]** #\<string\> the [file] includes three columns: the chromosome identifier, starting bp position, and ending bp position (one chromosome per line). 

[\[top\]](#ibdkin)


## 2.3 Optional Parameters

Each option is followed by its default value


* **--out ./** # \<string\> output prefix.
* **--nthreads 2** #\<int\> number of threads.
* **--kinship 0.0** #\<float\> minimum output kinship coefficient. 
* **--binkb 1000.0** #\<float\> bin size in kbp to calculate IBD coverage.
* **--fold 4.0** #\<float\> max permitted fold deviation from the genome-wide median IBD coverage.  Regions having greater deviation will be excluded in the kinship estimation.
* **--part 1 1**   #\<int\> \<int\>total partitions and current partition. The first integer is the total number of partitions that will be analyzed, and the second integer (starting from 1) determines the partition that will be analyzed in this run. This option enables distributed analysis across multiple computation nodes. Please see the [Distributed analysis](#213-distributed-analysis) section for an example.
* **--cutcm 4.0 2.0** #\<float\> \<float\> The minimum long and short IBD segment cM lengths. The first float is the minimum long IBD segment length, and the second float is the minimum short IBD segment length. A kinship coefficient is estimated for each pair of individuals having at least one long IBD segment. All short and long IBD segments are included when estimating the kinship coefficient in pairs of individuals having more than one long IBD segment. 
* **--merge 5.0 20.0**   #\<float\> \<float\> max cM merge lengths for IBD1 and IBD2 regions. The first float is the maximum length of the IBD0 region between two IBD1 regions that will be merged, the second float is the maximum length of the non-IBD2 region between two IBD2 regions that will be merged.

[\[top\]](#ibdkin)

## 2.4 Flags
* **--outmask** # output masked regions.
* **--outcoverage** # output IBD coverage.

[\[top\]](#ibdkin)


# 3 Output Files

**IBDkin** has three outputs: kinship coefficients, IBD coverage, and masked regions.
These outputs are reported in three gzip-compressed files.

## 3.1 Kinship Coefficients

The kinship coefficient output file (ending with ".kinship.gz") has eight tab-delimited columns:

1. Sample identifier-1
2. Sample identifier-2
3. Number of IBD segments used for calculation
4. IBD0 proportion
5. IBD1 proportion
6. IBD2 proportion
7. Kinship coefficient
8. Degree of relationship (values: 0, 1, 2, 3, >3)


[\[top\]](#ibdkin)

## 3.2 IBD Coverage

The IBD coverage output file (ending in ".coverage.gz") has four tab-delimited columns:

1. Chromosome
2. Start position
3. End position
4. IBD coverage

The IBD coverage, the forth column in the output, is the total number of IBD segments intersecting the chromosome interval, with each segment weighted by the proportion of the interval that it covers.


[\[top\]](#ibdkin)

## 3.3 Masked Regions

The masked regions output file (ending in "mask.gz") has three tab-delimited columns:

1. Chromosome
2. Start position of each masked region
3. End position of each masked region

[\[top\]](#ibdkin)

# 4 License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

[\[top\]](#ibdkin)
