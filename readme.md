# IBDkin
**IBDkin** is an **IBD** based **kin**ship estimation software, which is scaled to billions of IBD segments detected among hundreds of thousand individuals.
IBDkin is able to handle more samples by applying '-part' to conduct analysis on separate individual pairs groups.

If you use this software in your published analysis, please cite:

> TBD

# Contents

-   [Installation](#installation)
-   [Running IBDkin](#running-ibdkin)
    -   [A running example](#a-running-example)
    -   [Help Documents](#help-documents)
    -   [Required Parameters](#required-parameters)
    -   [Optional Parameters](#optional-parameters)
    -   [Output Flags](#output-flags)
-   [Output](#output)
    -   [Kinship](#kinship)
    -   [Coverage](#coverage)
    -   [Mask region](#mask-region)
-   [License](#license)


# Installation

Download the source code folder `src` and `cd` into it, then type 'make', and the executable file `IBDkin` will be generated.
This software is developed under linux CentOS 7.5, it has not been tested under any other systems.


[\[top\]](#content)

# Running IBDkin

## A running example

An example can be found in [run.sh](example.pub/run/run.sh).

Here list the bash script:
```bash
map="plink.map"
ind="ind.txt"
ibdfile="ibdhead.txt"

./IBDkin -ibdfile ${ibdfile} -map ${map} -ind ${ind} -nthreads 5 -out test --outmask --outcoverage
```

[\[top\]](#content)


## Help Document
Type `./IBDkin` or input error would invoke the help document:

```bash
Usage: IBDkin [options] parameters

     (Required inputs:)
     	-ibdfile [file]	# <string> a list of pathnames to IBD segments
     	-map [file]	# <string> genetic map in plink format
     	-ind [file]	# <string> individuals to be analyzed

     (Optional parameters:)
     	-out ./		# <string> output prefix
     	-nthreads 2 	#<int> number of threads
     	-degree 9 	#<int> degree cutoff to output kinship coefficients
     	-binkb 1000 	#<float> bin size in kbp to calculate IBD coverage
     	-fold 4 	#<float> fold number to determine the mask region
     	-part 1 1 	#<int> numbers of total partitions and the running partition
     	-cutcm 4 2 	#<float> IBD length filters in cM. 
     	-gapfill 5 20	#<float> length cutoffs for IBD gap filling. 

     (Other flags:)
     	--nokinship	# do not output kinship coefficients
     	--outmask	# output genome mask
     	--outcoverage	# output IBD coverage
```

[\[top\]](#content)

## Required Parameters

* **-ibdfile [file]** #<string> the [file] includes a list of pathnames to IBD segments in [hap-IBD format](https://github.com/browning-lab/hap-ibd). Each line is a pathname to the IBD segments called from one chromosome, every line represents different chromosome.
* **-ind [file]** #<string> the [file] includes a list of individuals to be analyzed.
* **-map [file]** #<string> the [file] is a [genetic map in plink format](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).

[\[top\]](#content)


## Optional Parameters

* **-out ./** #<string> output prefix
* **-nthreads 2** #<int> number of threads
* **-degree 9** #<int> degree cutoff to output kinship coefficients
* **-binkb 1000** #<float> bin size in kbp to calculate IBD coverage
* **-fold 4** #<float> fold number of the median coverage to determine the mask region. Both the regions covered by too much IBD segments, higher than 4 times of the genome-wide median, and the regions covered by rare IBD segments, lower than 1/4 of the genome-wide median, will be regarded as the *masked region* in further analysis.
* **-part 1 1**   #<int> the first integer determines the number of partitions that individual pairs are divided into, and the second integer (starting from 1) determines which partition to analyze in this run. This option allows us to distribute calculation to different computation nodes. See more details in the manuscript.
* **-cutcm 4 2** #<float> IBD length filters in cM. The individuals pairs sharing segments longer than the first float value will be analyzed and IBD segments longer than the second float will be used for kinship calculation. The way we calculate the kinship coefficients is different between the pairs sharing one long IBD segment and the pairs sharing multiple long IBD segments. See more details in the manuscript.
* **-gapfill 5 20**   #<float> length cutoffs for IBD gap filling. The first float is the cutoff for filling the gap between IBD1 segments and the second float is the cutoff for filling the gap between IBD2 segments.

[\[top\]](#content)

## Output Flags
* **--nokinship** # do not output kinship coefficients.
* **--outmask** # output genome mask.
* **--outcoverage** # output IBD coverage.

[\[top\]](#content)


# Output

**IBDkin** has three output files, kinship estimation, IBD coverage, and masked regions.
All these three files are gzip-compressed, ended with ".kinship.gz", ".coverage.gz", and ".mask.bed.gz"" separately.

## Kinship

Kinship coefficient output is includes eight columns:

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

Mask output records the regions to be excluded in the estimation, it includes three columns:

1. chromosome
2. start position
3. end position

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
