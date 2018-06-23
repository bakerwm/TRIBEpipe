.. _reading:


Reading paper
==============


FAQ: samtools mpileup output
------------------------------

Question: In some case, the depth calculated by the ``samtools mpileup`` is lower than ``bedtools genomecov -bz``, what's the reason?

Answer-1: There are some default filtering steps in samtools:

- Samtools mpileup will ignore bases with Q < 13

Use ``-Q 0`` to suppress the filtering

see: 

https://www.biostars.org/p/60841/#86579, 

https://www.biostars.org/p/60841/#60847

https://www.biostars.org/p/67579/#67791

- Samtools mpileup will remove duplicates or secondary reads:

Use ``--ff 4`` to request ONLY filtering of unmapped reads

https://github.com/samtools/samtools/issues/408


Question: What is the mean of characters in column-5 of ``samtools mpileup``

see samtools_ and wiki_  https://en.wikipedia.org/wiki/Pileup_format

.. _samtools: http://samtools.sourceforge.net/pileup.shtml

.. _wiki: https://en.wikipedia.org/wiki/Pileup_format


chr4    15685   a       28      >>>>>>>>..>...............>>    IJ:C<JCJJJJJJJGJIJIJ9JHIJJJJ    "!!!"~~!!!!!!!!!!!!!!!!!!!!!

from Wikipedia

	The columns
	Each line consists of 5 (or optionally 6) tab-separated columns:

	- 1 Sequence identifier  
	- 2 Position in sequence (starting from 1)  
	- 3 Reference nucleotide at that position  
	- 4 Number of aligned reads covering that position (depth of coverage)  
	- 5 Bases at that position from aligned reads  
	- 6 quality of those bases (OPTIONAL)  

	Column 5: The bases string

	- . (dot) means a base that matched the reference on the forward strand
	- , (comma) means a base that matched the reference on the reverse strand
	- </> (less-/greater-than sign) denotes a reference skip. This occurs, for example, if a base in the reference genome is intronic and a read maps to two flanking exons. If quality scores are given in a sixth column, they refer to the quality of the read and not the specific base.
	- AGTCN denotes a base that did not match the reference on the forward strand
	- agtcn denotes a base that did not match the reference on the reverse strand
	- A sequence matching the regular expression \+[0-9]+[ACGTNacgtn]+ denotes an insertion of one or more bases starting from the next position
	- A sequence matching the regular expression -[0-9]+[ACGTNacgtn]+ denotes a deletion of one or more bases starting from the next position
	- ^ (caret) marks the start of a read segment and the ASCII of the character following `^' minus 33 gives the mapping quality
	- $ (dollar) marks the end of a read segment
	- * (asterisk) is a placeholder for a deleted base in a multiple basepair deletion that was mentioned in a previous line by the -[0-9]+[ACGTNacgtn]+ notation
	
	Column 6: The base quality string

	This is an optional column. If present, the ASCII value of the character minus 33 gives the mapping Phred quality of each of the bases in the previous column 5. This is similar to quality encoding in the FASTQ format.



from SAMTools:

	where each line consists of chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases and base qualities. At the read base column, a dot stands for a match to the reference base on the forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch on the forward strand and `acgtn' for a mismatch on the reverse strand. A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. Here is an example of 2bp insertions on three reads:




Databases for RNA editing
~~~~~~~~~~~~~~~~~~~~~~~~~~~

+ RADAR [#]_ : http://rnaedit.com/

RADARa rigorously annotated database of A-to-I RNA editing.


+ REDIportal [#]_ : http://srv00.recas.ba.infn.it/atlas/

REDIportal, the largest and comprehensive collection of RNA editing in humans including more than 4.5 millions of A-to-I events detected in 55 body sites from thousands of RNAseq experiments.



**Reference**

.. [#] Gokul Ramaswami, Jin Billy Li; RADAR: a rigorously annotated database of A-to-I RNA editing, Nucleic Acids Research, Volume 42, Issue D1, 1 January 2014, Pages D109D113, https://doi.org/10.1093/nar/gkt996

.. [#] Ernesto Picardi, Anna Maria D'Erchia, Claudio Lo Giudice, Graziano Pesole; REDIportal: a comprehensive database of A-to-I RNA editing events in humans, Nucleic Acids Research, Volume 45, Issue D1, 4 January 2017, Pages D750D757, https://doi.org/10.1093/nar/gkw767



Tools
~~~~~~

+ REDItools [#]_

python scripts for RNA editing detection by RNA-Seq data

sourceforge: https://sourceforge.net/projects/reditools/files/?source=navbar

Github: https://github.com/ProfSmiles/reditools-updated


.. [#] Ernesto Picardi, Graziano Pesole; REDItools: high-throughput RNA editing detection made easy, Bioinformatics, Volume 29, Issue 14, 15 July 2013, Pages 18131814, https://doi.org/10.1093/bioinformatics/btt287


+ RES-Scanner [#]_ 

Here we present RES-Scanner, a flexible and efficient software package that detects and annotates RNA-editing sites using matching RNA-seq and DNA-seq data from the same individuals or samples.

.. [#] Zongji Wang, Jinmin Lian, Qiye Li, Pei Zhang, Yang Zhou, Xiaoyu Zhan, Guojie Zhang; RES-Scanner: a software package for genome-wide identification of RNA-editing sites, GigaScience, Volume 5, Issue 1, 1 December 2016, Pages 19, https://doi.org/10.1186/s13742-016-0143-4
















RNA-editing analysis
----------------------


- depth >= 20 in each replicate  

- A% >= 80% and G% = 0 in genomic DNA

- G% >= 10% in mRNA




From HyperTRIBE [#]_


    The criteria for RNA-editing events were as follows: (i) The nucleotide is covered by a minimum of 20 reads in each replicate; (ii) more than 80% of genomic DNA reads at this nucleotide are A with zero G (use the reverse complement if annotated gene is in the reverse strand); (iii) a minimum of 10% G is observed at this site in mRNA (or C for the reverse strand). Genomic DNA of S2 cells and background fly strain were sequenced to identify and exclude possible polymorphisms on the DNA level. RNA sequencing data were analyzed as previously described (Rodriguez et al. 2012; McMahon et al. 2016), with minor modifications. Background editing sites found in samples expressing Hyper-ADARcd alone were subtracted from the TRIBE identified editing sites both in S2 cells and in fly neurons. Overlap of editing sites from two data sets was identified using “bedtools intersect” with parameters “-f 0.9 -r”.


Reference:

.. [#] Joseph Rodriguez, Jerome S. Menet, Michael Rosbash, Nascent-Seq Indicates Widespread Cotranscriptional RNA Editing in Drosophila, Molecular Cell, Volume 47, Issue 1, 2012, [url1_]

.. _url1: https://www.sciencedirect.com/science/article/pii/S1097276512003541?via%3Dihub#sec4


From Rodriguez., 2012 Mol Cell [#]_

    RNA Editing Site Identification
    Base frequencies were calculated within exons and introns of UCSC annotation. Genes with multiple isoforms were flattened, where overlapping exons generate one exon. Base positions with one or more Gs in the nascent data sets and zero Gs in the sequenced genomic DNA were identified. We required that editing sites occur in at least five of six samples within each set of replicate time points, for a total of ten of 12 independent occurrences for each site. To avoid potential mismapping of reads at splice junctions by Tophat, we required that edited sites occur in at least one of the two middle quadrants of at least one read. Intronic sites that occurred within ten bases of an annotated splice site were also discarded. We also ranked our data sets by using the following g test log likelihood metric:

    When considering only zero G bases in the genomic, the equation simplifies to

    This metric was calculated for each sample and summed over all 12 replicates with 12 degrees of freedom. A one sided chi squared p value was generated and used to separate editing sites into high ranking with a 1 × 10−6 cutoff and low ranking with a 0.05 cutoff. We applied a similar approach to the yw pA-seq data, with the exception that we required that the editing site occur in both samples (two of two).

    Editing levels of our identified Nascent sites were also calculated in six paired-end 36 bp sequenced mRNA samples, as well as six single-end 72b p sequenced mRNA samples of the fly strain Cs. Editing levels were then calculated for the sites found in the nascent analysis.

    Determination of Editing Level
    For reproducibility of the nascent level, we pooled the editing counts for each replicate set of six time points together and calculated the percent editing level. The final percent editing level was determined by pooling the editing counts for all 12 samples and dividing by the total pooled counts of the 12 samples. Editing was similarly calculated for the yw pA-seq data and the small-scale Nascent-seq data by pooling of both replicate samples. The editing level for the Cs pA-seq data was calculated by pooling of all 12 samples. A two-tailed paired t test was used to test the significance of the observed editing levels between the nascent and yw mRNA data (p = 1 × 10−21).

.. [#] Joseph Rodriguez, Jerome S. Menet, Michael Rosbash, Nascent-Seq Indicates Widespread Cotranscriptional RNA Editing in Drosophila, Molecular Cell, Volume 47, Issue 1, 2012, Pages 27-37, ISSN 1097-2765, [url2_] 

.. _url2: https://doi.org/10.1016/j.molcel.2012.05.002.
