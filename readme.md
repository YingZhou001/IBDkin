
# IBDkin

**IBDkin** is a software for **IBD**-based **kin**ship estimation. 
IBDkin scales to hundreds of billions of IBD segments detected in hundreds of thousand individuals.

If you use this software in your published analysis, please cite:

> IBDkin: fast estimation of kinship coefficients from identity-by-descent segments

Last update: \today


Content
--------

- [1 Installation](#1-installation)
- [2 Running IBDkin](#2-running-ibdkin)
    - [2.1 Example analysis](#2.1-example-analysis)
    - [2.2 Required Parameters](#2.2-required-parameters)
    - [2.3 Optional Parameters](#2.3-optional-parameters)
    - [2.4 Flags](#2.4-flags)
- [3 Output Files](#3-output-files)
    - [3.1 Kinship coefficients](#3.1-kinship-coefficients)
    - [3.2 IBD coverage](#3.2-ibd-coverage)
    - [3.3 Masked regions](#3.3-masked-regions)
- [4 License](#4-license)



# 1 Installation

Download the source code `cd` into the source code folder "IBDkin/src-v2.8.7.3/", then type `make`, and the executable file `IBDkin` will be generated and ready to use.

```bash
git clone https://github.com/YingZhou001/IBDkin.git

cd IBDkin/src-v2.8.7.3/
make
```
This software is developed under linux CentOS 7.5.


[\[top\]](#ibdkin)

# 2 Running IBDkin

To run following examples, user can `cd` into the folder "IBDkin/example.pub/run", `cp` the executable file "IBDkin" into this folder, and type `sh run.sh` to see the outputs of following examples.

```bash
cd IBDkin/example.pub/run
cp ../../src-v2.8.7.3/IBDkin ./
sh run.sh
```
Otherwise, new users can follow the section 2.1.1 to 2.1.3 to run this software.

[\[top\]](#ibdkin)

## 2.1 Example Analysis

Let's take the first run of **IBDkin** in the folder [IBDkin/example.pub/run](example.pub/run).
In this folder, we have three required input files: "ibd.txt" is the IBD segments, "ind.txt" is the individual list, and "plink.map" is the recombination map. 
Detail description can be found in the section [Required Parameters](#required-parameters). 
Before running the following scripts, we need to copy the executable file `IBDkin` to this folder with command:

```bash
cp ../../src-v2.8.7.3/IBDkin ./
```

### 2.1.1 Computing kinship coefficients
We first run **IBDkin** with five threads (`-nthreads 5`), 
and output kinship coefficients between individual pairs whose relationship is up to the 9th degree (the default). 
The output file prefix is specified with the `-out` flag.
The output kinship coefficients are in the gzip-compressed file "example-1.kinship.gz":

```bash
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-1
```

### 2.1.2 Computing IBD coverage and masked regions
In certain situations, we may want to output the masked regions (flag `--outmask`) and the IBD coverage across the genome (flag `--outcoverage`). 
If we do not need to perform kinship estimation, we can make the analysis faster by adding the `--nokinship` flag.  When the `--nokinship` flag is specified, the program makes only one pass through the IBD data to calculate the mask regions and coverages.
The output mask filename will end in ".mask.gz" and the output coverage filename will end in ".coverage.gz". 
Both files are gzip-compressed. 
The command for this analysis is:

```bash
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-2 --outmask --outcoverage --nokinship
```

### 2.1.3 Distributed analysis
If there is not enough memory to run an analysis on a huge data set, we can use the option `-part` to distribute computation to different nodes.
In the following script, we assign the analysis into five partitions and run each partition separately by assigning integers, from 1 to 5, to the variable `${part}`.


```bash
part=1 # can be set as 2, 3, 4, and 5
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-3.${part} -part 5 ${part}
```
**Caution!**: specify different output prefix for each part when use the `-part` option, or you will overwrite your estimations.

For advanced parameters setting, please read our [manuscript(add link to our manuscript)](#ibdkin) and check the following sections: [Required Parameters](#required-parameters), [Optional Parameters](#optional-parameters).

[\[top\]](#ibdkin)


## 2.2 Required Parameters

* **-ibdfile [file]** #\<string\> > the [file] contains the pathnames of files of the IBD segments on each chromosome. One pathname per line, and one line per chromosome. The IBD segments are stored in gzip-compressed [hap-IBD format](https://github.com/browning-lab/hap-ibd).
* **-map [file]** #\<string\> the [file] is a [genetic map with cM distances in PLINK format](http://zzz.bwh.harvard.edu/plink/data.shtml), including all target chromosomes used for kinship estimation. 
The chromosome identifier in each hap-IBD output file must match the chromosome identifier in the genetic map file.
* **-ind [file]** #\<string\> the [file] includes a list of individuals to be analyzed.

[\[top\]](#ibdkin)


## 2.3 Optional Parameters

(default value(s) follow each option)

* **-out ./** # \<string\> output prefix.
* **-nthreads 2** #\<int\> number of threads.
* **-degree 9** #\<int\> max relationship degree for output kinship coefficients.
* **-binkb 1000** #\<float\> bin size in kbp to calculate IBD coverage.
* **-fold 4** #\<float\> max fold permitted deviation from the genome-wide median IBD coverage.  Regions having greater deviation will be excluded in the kinship estimation.
* **-part 1 1**   #\<int\> \<int\>total partitions and current partition. The first integer is the total number of partitions we are going to analyze separately, and the second integer (starting from 1) determines which partition to analyze in this run. This option allows us to distribute calculation to different computation nodes. Please check the section [Distributed analysis](#distributed-analysis) for an working example.
* **-cutcm 4 2** #\<float\> \<float\> The minimum long and short IBD segment lengths in cM. A kinship coefficient is estimated for each pair of individuals having at least one long IBD segment.  All short and long IBD segments are included when estimating the kinship coefficient in pairs of individuals having more than one long IBD segment. 
* **-merge 5 20**   #\<float\> \<float\> max cM merge lengths for IBD1 and IBD2 regions. The first float is the maximum length of the IBD0 region between IBD1 regions to be merged, the second float is the maximum length of the non-IBD2 region between IBD2 regions to be merged.

[\[top\]](#ibdkin)

## 2.4 Flags
* **--nokinship** # do not output kinship coefficients.
* **--outmask** # output genome mask.
* **--outcoverage** # output IBD coverage.

[\[top\]](#ibdkin)


# 3 Output Files

**IBDkin** has three outputs: kinship coefficients, IBD coverage, and masked regions.
These outputs are reported in three gzip-compressed files.

## 3.1 Kinship coefficients

The kinship coefficient output file (ending with ".kinship.gz") has eight columns:

1. Sample identifier-1
2. Sample identifier-2
3. Number of IBD segments used for calculation
4. IBD0 proportion
5. IBD1 proportion
6. IBD2 proportion
7. Kinship coefficient
8. Degree of relationship


[\[top\]](#ibdkin)

## 3.2 IBD coverage

The IBD coverage output file (ending in ".coverage.gz") has four columns:

1. Chromosome
2. Start position
3. End position
4. IBD coverage

The IBD coverage, the forth column in the output, is the total number of IBD segments intersecting the chromosome interval, with each segment weighted by the proportion of the interval that it covers.


[\[top\]](#ibdkin)

## 3.3 Masked regions

The masked output file (ending in "masked.gz") has three columns:

1. Chromosome
2. Start position of each mask region
3. End position of each mask region

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
