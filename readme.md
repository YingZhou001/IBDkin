
# IBDkin
**IBDkin** is an **IBD** based **kin**ship estimation software, which is scaled to hundreds of billions of IBD segments detected among hundreds of thousand individuals.

If you use this software in your published analysis, please cite:

> TBD

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


# Installation

Download the source code folder `src` and `cd` into it, then type `make`, and the executable file `IBDkin` will be generated and ready to use.

```bash
#commands to download and install this software
#will be adjust once the software is uploaded
wget ..
tar xvf ..
make
```
This software is developed under linux CentOS 7.5.


[\[top\]](#content)

# Running IBDkin

## Example analysis

Let's take the first bite of **IBDkin** with the data in the [example.pub](example.pub/).

```bash
# initialize the inputs
map="plink.map"
ind="ind.txt"
ibdfile="ibdhead.txt"
```
In this script, we run **IBDkin** with five threads (`-nthreads 5`), 
and output kinship coefficients between individual pairs whose relatioship is closer than the 10th degree (default). 
The output will saved in gzip compressed file "first-bite.kinship.gz" as the prefix is specified by `-out first-bite`.

```bash
./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out first-bite --outmask --outcoverage
```

Under a certain situation, we want to output the mask regions (add flag `--outmask`), the IBD coverages (add flag `--outcoverage`), and without the kinship estimation (add flag `--nokinship`).
The flag `--nokinship` will make the computation more efficient becasue it only takes one pass of the data to calculate the mask regions and coverages.
The mask file would be "first-bite.mask.gz" and the coverage file would be "first-bite.coverage.gz", both of them are gzip compressed.

```bash
./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out first-bite --outmask --outcoverage --nokinship
```
When the memory is not enough for running analysis on a huge data set, we can use the option `-part` to distribute the computation to the different nodes or to allow several times running on the same node. 
In the bellowing command line, we assign the analysis into five partitions and run each partition separately. 


```bash
for part in {1..5}
do
./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out first-bite.${part} -part 5 ${part}
done
```
**Caution!**: specify different prefix names when use the option `-part`, or you will miss your estimations.

For advanced parameters setting, please read our [manuscript(add link to our manuscript)](#ibdkin) and check the following sections: [Required Parameters](#required-parameters), [Optional Parameters](#optional-parameters).

[\[top\]](#content)


## Required Parameters

(to Joe, compresssion input option?)

* **-ibdfile [file]** #\<string\> the [file] of pathnames, line by line, to each input IBD segments called from different chromosomes. The IBD segments are stored in compressed [hap-IBD format](https://github.com/browning-lab/hap-ibd).
* **-map [file]** #\<string\> the [file] is a [genetic map in plink format](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/), including all target chromsomes used for kinship estimation.(to Brian, how do you deal with the different chromosome identifiers? for example chr1 vs 1?)
* **-ind [file]** #\<string\> the [file] includes a list of individuals to be analyzed.

[\[top\]](#content)


## Optional Parameters

(default value(s) are followed after the options)

* **-out ./** # \<string\> output prefix
* **-nthreads 2** #\<int\> number of threads
* **-degree 9** #\<int\> degree cutoff to output kinship coefficients
* **-binkb 1000** #\<float\> bin size in kbp to calculate IBD coverage
* **-fold 4** #\<float\> total partitions and current partition. Both the regions covered by too much IBD segments or too little IBD segments are defined as the *masked region*. 
* **-part 1 1**   #\<int\> \<int\> total partitions and current partition. The first integer is the total number of partitions we are going to analyze separately, and the second integer (starting from 1) determines which partition to analyze in this run. This option allows us to distribute calculation to different computation nodes. See more details in the manuscript.
* **-cutcm 4 2** #\<float\> \<float\> IBD length filters in cM. The first float determines which individual pairs to be analyzed and the second float determines which IBD segments to be used for kinship calculation.
* **-gapfill 5 20**   #\<float\> \<float\> max cM length for IBD gap filling. The first float is the maximun length of the gap between IBD1 segments to be filled, the second float is the maximum length of the gap between IBD2 segments to be filled.

[\[top\]](#content)

## Flags
* **--nokinship** # do not output kinship coefficients.
* **--outmask** # output genome mask.
* **--outcoverage** # output IBD coverage.

[\[top\]](#content)


# Output

**IBDkin** has three outputs: kinship coefficients, IBD coverage, and masked regions.
All these three files are gzip-compressed, ended with "kinship.gz", "coverage.gz", and "mask.gz"" separately.

## Kinship

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

## Coverage

IBD coverage output includes four columns:

1. chromosome
2. start position
3. end position
4. IBD coverage

[\[top\]](#content)

## Mask region

Mask output includes three columns:

1. chromosome
2. start position of each mask region
3. end position of each mask region

[\[top\]](#content)

# License

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
