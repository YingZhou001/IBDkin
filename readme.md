
# IBDkin

**IBDkin** is a software for **IBD** based **kin**ship estimation.  
IBDkin scales to hundreds of billions of IBD segments detected in hundreds of thousand individuals.

If you use this software in your published analysis, please cite:

> IBDkin: fast estimation of kinship coefficients from identity-by-descent segments

Last update: \today

# Content

-   [Installation](#installation)
-   [Running IBDkin](#running-ibdkin)
    -   [Example analysis](#example-analysis)
    -   [Required Parameters](#required-parameters)
    -   [Optional Parameters](#optional-parameters)
    -   [Flags](#flags)
-   [Output](#output)
    -   [Kinship](#kinship)
    -   [Coverage](#coverage)
    -   [Mask region](#mask-region)
-   [License](#license)


# 1 Installation

Download the source code folder `src` and `cd` into it, then type `make`, and the executable file `IBDkin` will be generated and ready to use.

```bash
#commands to download and install this software
#will be adjust once the software is uploaded
git clone https://github.com/YingZhou001/IBDkin.git

cd IBDkin/src-v2.8.7.3/
make
```
This software is developed under linux CentOS 7.5.


[\[top\]](#content)

# 2 Running IBDkin

To run following examples, user can `cd` into the folder `IBDkin/example.pub/run`, `cp` the executable file `IBDkin` into this folder, and type `sh run.sh` to see the outputs of following examples.

```bash
cd IBDkin/example.pub/run
cp ../../src-v2.8.7.3/IBDkin ./
sh run.sh
```
Otherwise, new users can follow the section 2.1.1 to 2.1.3 to run this software.

## 2.1 Example Analysis

Let's take the first run of **IBDkin** in the folder [IBDkin/example.pub/run](example.pub/run).
In this folder, we have three required input files: `ibd.txt` is for IBD segments, `ind.txt` is for the individual list, and `plink.map` is the recombination map. Clear description can be found in the section [Required Parameters](#required-parameters). 
Before running the following scripts, we need to copy the executable file `IBDkin` to this folder with command:

```bash
cp ../../src-v2.8.7.3/IBDkin ./
```

### 2.1.1 Computing kinship coefficients
In this first script, we will run **IBDkin** with five threads (`-nthreads 5`), 
and output kinship coefficients between individual pairs whose relationship is up to the 9th degree (default). 
The output file prefix is specified with the `-out` flag.
The output kinship coefficients are in the gzip-compressed file "example-1.kinship.gz":

```bash
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-1
```

### 2.1.2 Computing IBD coverage and masked regions
In certain situations, we may want to output the masked regions (flag `--outmask`) and the IBD coverage across the genome (flag `--outcoverage`) without performing kinship estimation (flag `--nokinship`). 
If we do not need to perform kinship estimation, we can make the analysis more efficient by adding the `--nokinship` flag.  When the `--nokinship` flag is specified, the program only needs to make one pass through the IBD data to calculate the mask regions and coverages.
This option will be very helpful when the IBD input size is huge, in which file reading consumes the major analysis time.
The output mask filename will end in ".mask.gz" and the output coverage filename will end in ".coverage.gz". 
Both files are gzip-compressed. 
The command for this analysis is:

```bash
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-2 --outmask --outcoverage --nokinship
```

### 2.1.3 Distributed analysis
If there is not enough memory to run an analysis on a huge data set, we can use the option `-part` to distribute computation to different nodes or to allow several times running on the same node. 
In the following script, we assign the analysis into five partitions and run each partition separately.


```bash
for part in {1..5}
do
./IBDkin -ibdfile ibd.txt -map plink.map -ind ind.txt -nthreads 5 -out example-3.${part} -part 5 ${part}
done
```
**Caution!**: specify different output prefix for each part when use the `-part` option, or you will overwrite your estimations.

For advanced parameters setting, please read our [manuscript(add link to our manuscript)](#ibdkin) and check the following sections: [Required Parameters](#required-parameters), [Optional Parameters](#optional-parameters).

[\[top\]](#content)


## 2.2 Required Parameters

* **-ibdfile [file]** #\<string\> > the [file] contains the pathnames of files of the IBD segments on each chromosome. One pathname per line, and one line per chromosome. The IBD segments are stored in compressed [hap-IBD format](https://github.com/browning-lab/hap-ibd).
* **-map [file]** #\<string\> the [file] is a [genetic map in plink format](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), including all target chromosomes used for kinship estimation. The chromosome identifiers in the genetic map and the hap-IBD output file for the same chromosome must have identical chromosome identifiers.
* **-ind [file]** #\<string\> the [file] includes a list of individuals to be analyzed.

[\[top\]](#content)


## 2.3 Optional Parameters

(default value(s) are followed after each option)

* **-out ./** # \<string\> output prefix.
* **-nthreads 2** #\<int\> number of threads.
* **-degree 9** #\<int\> max relationship degree for output kinship coefficients.
* **-binkb 1000** #\<float\> bin size in kbp to calculate IBD coverage.
* **-fold 4** #\<float\> max fold permitted deviation from the genome-wide median IBD coverage.  Regions having greater deviation will be excluded in the kinship estimation.
* **-part 1 1**   #\<int\> \<int\>total partitions and current partition. The first integer is the total number of partitions we are going to analyze separately, and the second integer (starting from 1) determines which partition to analyze in this run. This option allows us to distribute calculation to different computation nodes. Please check the section [Distributed analysis](#distributed-analysis) for an working example.
* **-cutcm 4 2** #\<float\> \<float\> The minimum long and short IBD segment lengths in cM. A kinship coefficient is estimated for each pair of individuals having at least one long IBD segment.  All short and long IBD segments are included when estimating the kinship coefficient in pairs of individuals having more than one long IBD segment. Please check the manuscript for how we use short and long segments to calculate kinship coefficients.
* **-merge 5 20**   #\<float\> \<float\> max cM merge lengths for IBD1 and IBD2 regions. The first float is the maximum length of the IBD0 region between IBD1 regions to be merged, the second float is the maximum length of the non-IBD2 region between IBD2 regions to be merged.

[\[top\]](#content)

## 2.4 Flags
* **--nokinship** # do not output kinship coefficients.
* **--outmask** # output genome mask.
* **--outcoverage** # output IBD coverage.

[\[top\]](#content)


# 3 Output

**IBDkin** has three outputs: kinship coefficients, IBD coverage, and masked regions.
All these three files are gzip-compressed, ended with "kinship.gz", "coverage.gz", and "mask.gz"" separately.

## 3.1 Kinship

Kinship coefficient output includes eight columns:

1. sample identifier-1
2. sample identifier-2
3. number of IBD segments used for calculation
4. IBD0 proportion
5. IBD1 proportion
6. IBD2 proportion
7. kinship coefficient
8. degree of relationship


[\[top\]](#content)

## 3.2 Coverage

IBD coverage output includes four columns:

1. chromosome
2. start position
3. end position
4. IBD coverage

The IBD coverage, the forth column in the output, is the total number of IBD segments intersecting the chromosome interval, with each segment weighted by the proportion of the interval that it covers.


[\[top\]](#content)

## 3.3 Mask region

Mask output includes three columns:

1. chromosome
2. start position of each mask region
3. end position of each mask region

[\[top\]](#content)

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

[\[top\]](#content)
