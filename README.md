# ChIP-seq tutorial

## Introduction

So you have ChIP-seq data. Peaks of transcription factors or histone modifications. You might also have genome accessibility data, such as ATAC-seq. The data span multiple conditions, timepoints, developmental stages. Now what?

This tutorial covers downstream analysis of ChIP-seq data in a multi-sample setting. It will assume you have already mapped this data and performed peak calling. There are many tutorials available that illustrate how to do this. See for instance:

* [EBI ChIP-seq practical](https://www.ebi.ac.uk/training/online/course/ebi-next-generation-sequencing-practical-course/chip-seq-analysis/chip-seq-practical)
* [Hands-on introduction to ChIP-seq analysis - VIB Training](http://www.biologie.ens.fr/~mthomas/other/chip-seq-training/)

### Learning outcomes

After finishing this tutorial you:
* Can perform k-means clustering of ChIP-seq data and visualize the results as a heatmap.
* Can interpret a heatmap of clustered ChIP-seq data.
* Can identify biologically relevant transcription factor motifs from ChIP-seq and/or chromatin accessibility data.


## Before you begin

### Software prerequisites

This tutorial uses [genomepy](https://github.com/simonvh/genomepy), [fluff](https://github.com/simonvh/fluff) and [GimmeMotifs](https://gimmemotifs.readthedocs.io). You can install these tools using [bioconda](http://bioconda.github.io/). You can use bioconda on Linux, Mac OSX and Windows 10.


### Installation instructions for Windows 10 only

If you use Windows 10 you have to install the Windows Subsystem for Linux. Follow these instructions: [https://docs.microsoft.com/en-us/windows/wsl/install-win10](https://docs.microsoft.com/en-us/windows/wsl/install-win10). Install the Ubuntu distribution. Now you can start Ubuntu and follow the rest of the installation instructions. You are now running Linux, so follow the Linux installation for Miniconda!

### Install software using conda

If you have not used bioconda before, read and follow the installation instructions [here](http://bioconda.github.io/).  
Pay attention to the bioconda channel order if you have used it before (it has been changed recently). According to the current bioconda documentation it should be:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Now you can install all required software with the following command:

```
conda create -n tutorial python=3 biofluff gimmemotifs=0.13
```

You can activate the environment you have created with the following command:

```
conda activate tutorial
```

### Obtain data

#### Use your own

You can use your own data if you wish. You will need the following:

* A BED file with regulatory element locations (enhancers). In this tutorial we will obtain them from DNase I MACS2 narrowPeak files. Ideally, these peaks are based on a source that can identify the location of regulatory elements with reasonable accuracy. This includes, for instance, ATAC-seq, DNase I or p300 ChIP-seq. While histone marks such as H3K27ac or H3K4me1 mark enhancers it is very hard to accurately identify the correct location of the regulatory elements. The enhancer ChromHMM tracks from ROADMAP or ENCODE are also not accurate enough.
* At least two BAM files with ChIP-seq data. For this tutorial we will use H3K27ac ChIP-seq which correlates with enhancer activity.

#### Download tutorial data

As example data in this tutorial we will use ROADMAP data from the Bing Ren lab described in this paper: [Xie et al., 2013](https://doi.org/10.1016/j.cell.2013.04.022). 

![Figure 1: Xie et al. Epigenome Reference Maps for hESCs and Four hESC-Derived Lineages](https://i.imgur.com/evrXF69.jpg)

We will use DNase I and H3K27ac ChIP-seq data of human embryonic stem cells (hESCs) differentiated *in vitro* to four different lineages: Mesenchymal Stem Cells (MSC), Neural Progenitor Cells (NPC), Trophoblast-like Cell (TBL) and Mesendoderm (ME).

Download the tutorial data from Zenodo:

[https://doi.org/10.5281/zenodo.1560640](https://doi.org/10.5281/zenodo.1560640)

**Please note:** that this archive contains only a small selection of the full data, to be able to run everything in reasonable time.

### Install genome data

We will use [genomepy](https://github.com/simonvh/genomepy) to manage the genome data that we need. The example data has been mapped to human genome hg19. The following command will download this genome (and gene annotation)from UCSC:

```
genomepy install hg19 UCSC --annotation
```

You can also install other genomes from UCSC, Ensembl or NCBI. For instance, to download the current zebrafish genome from Ensembl (not needed for this tutorial):

```
genomepy install GRCz11 Ensembl --annotation
```

## Clustering and visualization

### Preparing the data

First extract the file `tutorial_data.tgz` that you downloaded from Zenodo.

```
tar xvzf tutorial_data.tgz
```

Let's start with creating a peak set to use for clustering. We will create a set of putative enhancers by overlapping the DNaseI peaks with the H3K27ac peaks. For these types of analysis it is important to use a combined peak set, that contains all peak regions of any one of the experiments. 

Let's start by combining the DNase I peak summits from different experiments:

```
combine_peaks -i data/peaks/GSM*summits.bed -g hg19 -w 200 > DNase_summits.bed
```

The script `combine_peaks` (included with GimmeMotifs >= 0.13) creates a merged peak summit file from a combination of MACS2 summit files. In case peaks overlap (with the peak size determined by the `-w` window parameter) the highest summit is selected. This means that the the resulting BED file has peaks that are always centered on a summit of one of the input files.


Now, let's select only the DNase I peaks that have H3K27ac signal, as a way of identifying putative enhancers. We will use the H3K27ac peaks as provided by ENCODE.

```
# Create a merged set of H3K27ac peaks
zcat data/peaks/ENC*H3K27ac*bed.gz | sort -k1,1 -k2g,2 | bedtools merge > H3K27ac.peaks.merged.bed
# Overlap DNase I with these peaks
bedtools intersect -a DNase_summits.bed -b H3K27ac.peaks.merged.bed > DNase_summits_with_H3K27ac.bed
```

Our final file `DNase_summits_with_H3K27ac.bed` contains the summits of DNase I peaks that overlap with H3K27ac peaks.

### Running fluff

You can create a heatmap of clustered ChIP-seq signals using [fluff](https://doi.org/10.7717/peerj.2209). 

For experimentation and explorative analysis it is often good to select a smaller subset of your data, so that everything runs quicker. Let's create a file with 5,000 randomly selected peaks:

```
shuf DNase_summits_with_H3K27ac.bed | head -n 5000 > DNase.H3K27ac.5k.peaks.bed
``` 

You can create a heatmap with the following command:

```
fluff heatmap -f DNase.H3K27ac.5k.peaks.bed -d data/bam/*bam -C kmeans -k 5 \
    -o DNase.H3K27ac.5k.heatmap.kmeans.5.png
```

This will select all reads in a 10kb window centered at your peaks and cluster this using kmeans clustering with 5 clusters. The resulting heatmap should look like this:

![](https://github.com/simonvh/bioinfosummer/raw/master/img/DNase.H3K27ac.5k.heatmap.kmeans.5.png)

Let's create nicer names for the BAM files:

```  
ln -s data/bam/H1_H3K27ac.10M.bam H1.bam
ln -s data/bam/H1_H3K27ac.10M.bam.bai H1.bam.bai
ln -s data/bam/H1_mesenchymal_H3K27ac.10M.bam MSC.bam
ln -s data/bam/H1_mesenchymal_H3K27ac.10M.bam.bai MSC.bam.bai
ln -s data/bam/H1_neuronal_progenitor_H3K27ac.10M.bam NPC.bam
ln -s data/bam/H1_neuronal_progenitor_H3K27ac.10M.bam.bai NPC.bam.bai
ln -s data/bam/H1_trophoblast_H3K27ac.10M.bam TBL.bam
ln -s data/bam/H1_trophoblast_H3K27ac.10M.bam.bai TBL.bam.bai
ln -s data/bam/H1_mesendoderm_H3K27ac.10M.bam ME.bam
ln -s data/bam/H1_mesendoderm_H3K27ac.10M.bam.bai ME.bam.bai
``` 

By using the linked files we get shorter names:

```
fluff heatmap -f DNase.H3K27ac.5k.peaks.bed -d H1.bam MSC.bam ME.bam NPC.bam TBL.bam  -C kmeans -k 5  -o DNase.H3K27ac.5k.heatmap.kmeans.5.example2.png
```

![](https://github.com/simonvh/bioinfosummer/raw/master/img/DNase.H3K27ac.5k.heatmap.kmeans.5.example2.png)


You can change the colors: 

```
fluff heatmap -f DNase.H3K27ac.5k.peaks.bed -d H1.bam MSC.bam ME.bam NPC.bam TBL.bam -C kmeans -k 5 -o DNase.H3K27ac.5k.heatmap.kmeans.5.example3.png -c green,green,red,#e41a1c,#ffff33 -b white,black,white,black,#984ea3
```

![](https://github.com/simonvh/bioinfosummer/raw/master/img/DNase.H3K27ac.5k.heatmap.kmeans.5.example3.png)


And also change the saturation / scaling:

```
fluff heatmap -f DNase.H3K27ac.5k.peaks.bed -d H1.bam MSC.bam ME.bam NPC.bam TBL.bam -C kmeans -k 5 -s 99% \
-o DNase.H3K27ac.5k.heatmap.kmeans.5.example4.png 
```

![](https://github.com/simonvh/bioinfosummer/raw/master/img/DNase.H3K27ac.5k.heatmap.kmeans.5.example4.png)


Let's go back to the intital clustering. 

**Q:** What kind of patterns, or clusters, would you expect to see? Do you see clusters with dynamic difference between the different cell types? 

This is what happens when ChIP-seq data, summarized in bins, is clustered using a the Euclidean distance metric. The resulting clusters are often *spatial*, i.e. ChIP-seq enrichment to the left and to the right of the center, or broad versus narrow peaks. While this may be interesting, often we are interested in dynamic patterns between different conditions. In `fluff` there is a way to specify dynamic clustering with a combination of two options: `-g` will only use one 1kb bin for each sample for clustering and `-M Pearson` will use the Pearson correlation distance as distance metric.

Let's have a look at the differences. 

```
fluff heatmap -f DNase.H3K27ac.5k.peaks.bed -d *bam -C kmeans -k 5 \
    -M Pearson -g -o DNase.H3K27ac.5k.heatmap.kmeans.5.dynamic.example5.png 
```

![](https://github.com/simonvh/bioinfosummer/raw/master/img/DNase.H3K27ac.5k.heatmap.kmeans.5.dynamic.example5.png)


**Q:** Do you see different clusters as compared to example 1? 

**Do:** Try running `fluff heatmap` with different cluster numbers, controlled by the `-k` parameter. Picking the correct number of clusters is not straightforward, and will depend on your data. Once you are satisified you can run `fluff heatmap` on the full dataset. Depending on the number of peaks, this might take some time.

## Motif Analysis using GimmeMotifs

### Motif enrichment in a single data set

Here, we will identify enriched motifs in (putative) enhancers of H1 embryonic stem cells. Feel free to use another data set of your choice!


#### Set up the genome variable

For some of the next questions we will need the location of the genome file. We will use the exports provided by `genomepy`. Find the config file of `genomepy` on your system:

```
genomepy config file
```

In my case this will show the following:

```
/home/simon/.config/genomepy/genomepy.yaml
```

Now replace `genomepy.yaml` with `exports.txt` and source this file:


```
source /home/simon/.config/genomepy/exports.txt
```

You now should have a shell variable `$HG19` (the genome name in upper-case) that points to the genome FASTA file. 

```
echo $HG19
```

```
/home/simon/.local/share/genomes/hg19/hg19.fa
```

If you want, you can add the `source` command above to your `~/.bashrc`. Then you will have a shell variable for every genome installed with `genomepy`.

#### Identify putative enhancers

From the DNase I peaks we have selected all regions that have H3K27ac. However, this will also include promoters. Let's focus on putative enhancers and remove all promoters from this set of peaks.

In the directory where the genome is located, you will also find an BED file containing gene annotation:

```
ls `dirname $HG19`
```

```
README.txt              hg19.annotation.gtf.gz  hg19.fa.fai    hg19.gaps.bed
hg19.annotation.bed.gz  hg19.fa                 hg19.fa.sizes
```

Let's use this to create a BED file of promoters (1kb around the TSS):

```
bedtools slop -l 1 -r 0 -s -i /data/genomes/hg19/hg19.annotation.bed.gz -g $HG19.sizes | \
    bedtools slop -b 500 -g $HG19.sizes | \
    cut -f1-6 > hg19.promoter.1kb.bed
``` 

Note the use of the `$HG19` shell variable in combination with `.sizes`. All genomes installed with `genomepy` will have this file, which contains a list of all chromosomes or scaffolds with their length. This can be used as a genome file for the `-g` option of all `bedtools` commands.

For our first motif analysis we will select peaks from H1 ES cells specifically. Using the next commmand we will:

* Select all peaks do not overlap with a promoter
* Select the peaks that overlap with H3K27ac in H1 cells
* Take the highest 1000 peaks (based on the value in the ENCODE peak file)
* Create 200bp sized regions

```
bedtools intersect -a DNase_summits_with_H3K27ac.bed -b hg19.promoter.1kb.bed -v | \
    bedtools intersect -a - -b data/peaks/ENCFF335JSA.H3K27ac.H1.bed.gz  -wo | 
    sort -k12gr  | head -n 1000 | cut -f1-3 | \
    bedtools slop -b 100 -g $HG19.sizes > top_enhancers.H1.bed
```

#### Run the motif analysis

For every motif analysis you need a background! One option is to select random regions from the genome. We can do this with `gimme background`. Generate a FASTA file with 5,000 sequences of length 200 from the hg19 genome:

```
gimme background background.random.w200.fa genomic -n 5000 -l 200 -f FASTA -g hg19
```

Now let's run our first motif analysis!

```
gimme roc top_enhancers.H1.bed background.random.w200.fa -g hg19 -r roc.report.random 
```

Note that you can use BED files to scan for motifs, as long as you supply a genome (GimmeMotifs needs to know where to get the sequences). Once the command has finished, you can open the html report in the output directory using your browser. If you ran this analysis remotely, retrieve the whole output directory!

You can sort the output by sorting on the columns header. Sort by 'Enrichmenty at 1% FPR'. 

**Q:** Notice something about these motifs?

Here you see an example of how the choice of background sequences will influence a motif analysis. Even though we removed all promoter peaks, which include CpG-island promoters, these putative enhancers apparently have a high GC%. Let's change backgrounds. The `gimme background` tool can also select sequences from the genome that match the GC% of your input sequences. For that, you will need to convert the peak BED file to FASTA using `bedtools`.

```
# Create FASTA file
bedtools getfasta -fi $HG19 -bed top_enhancers.H1.bed -fo top_enhancers.H1.fa
# Create GC matched background
gimme background background.gc.w200.fa gc -i top_enhancers.H1.bed -n 5000 -l 200 -f FASTA -g hg19
# Run the motif analysis again
gimme roc top_enhancers.H1.bed background.gc.w200.fa -g hg19 -r roc.report.gc 
```

**Q:** What are the top motifs now? Do these results make sense?

### Comparison of motifs in multiple data sets using maelstrom

We will create a tab-separated table with H3K27ac read counts (log-transformed and scaled). For H3K27ac ChIP-seq data we will use relatively broad regions (2kb). If you have other data that gives sharp peaks, such as p300 ChIP-seq, you would use a smaller window. For this tutorial we select only top 1000 peaks (as determined by the variance), as the analysis would otherwise take too long. For a normal analysis you could use more peaks, although this would influence the running time of `gimme maelstrom`. In my experience ~10,000 peaks would usually already give good results. If you have some other way of pre-selecting informative, dynamic or "interesting" peaks you can use that. 

```
coverage_table -p  DNase.H3K27ac.5k.peaks.bed -d *bam -w 2000 -l -s -t 1000 -T var > coverage_table.txt
``` 

Now we can analyse this file using `gimme maelstrom`. Normally, you would not have to specify the `-m` parameter, but one of the methods takes much longer to complete. By specifying these three methods, you will get results within a reasonable timeframe. With the `-N` argument you can specify the number of cores to use. Make this not higher than the number of cores in your system. There is a certain memory usage, if `gimme maelstrom` starts to use too much memory you would have to reduce the number of cores.
With four cores this command should take ~9 minutes.
```
gimme maelstrom coverage_table.txt hg19 maelstrom.out -N 4 -m bayesianridge,lightningregressor,xgboost
```

Once again the output directory `maelstrom.out` will contain an `html` file that you can open in your browser `maelstrom.out/gimme.maelstrom.report.html`.

You can change the ordering of the motifs by clicking on the header. 

**Q:** Check a few cell types. Do the motifs make sense (search the literature if you are not familiar with these cell types)?

Unfortunately, there is no easy way to check if the motifs are relevant. Usually domain knowledge or extensive literature search is necessary to follow up. This computational analysis would serve as a starting point for further analysis. Maybe the way forward would be some experiments?

# License

This material is made available under the [Creative Commons BY 4.0 license](https://creativecommons.org/licenses/by/4.0/). Feel free to share, remix and reuse!



