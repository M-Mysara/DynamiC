
				||||||||||||||||||||||||||||||||||||||||||||||||
				||              Welcome To DynamiC	          ||
				|| Software for constructing a lookup table   ||
				||		for more accurate OTU clustering	  ||
				||    Copyright (C) 2016  <M.Mysara et al>    ||
				||||||||||||||||||||||||||||||||||||||||||||||||

				 DynamiC version 1, Copyright (C) 2016, M.Mysara et al
				 DynamiC comes with ABSOLUTELY NO WARRANTY.
				 This is free software, and you are welcome to redistribute it under
				 certain conditions; please refer to \'COPYING\' for details.;
				 The software also includes \"mothur\" also under GNU Copyright
###Introducing Dynamic Cutoff (**D**ynami**C**)
####Adjustable cutoff for OTU clustering
#####The development of high-throughput sequencing technologies has revolutionized the field of microbial ecology via 16S RNA gene amplicon sequencing approaches. Clustering those amplicon sequencing data into Operational Taxonomic Units (OTUs) is one of the most commonly used approaches to approximate a bacterial species. Since a 97% 16S rRNA sequence similarity has been widely used in bacterial taxonomy as one of the criteria to delineate species, this cut-off is often applied when clustering amplicon reads into OTUs. However, where this cut-off is derived based on full-length 16S rRNA genes, the amplicons obtained with current high-throughput sequencing approaches in general only rely on one or two variable regions within this 16S rRNA gene. Therefore, within this work we assess the paradigm that applying a clustering step using a sequence similarity cut-off of 97% would lead to OTUs accurately corresponding to species. We show that the robustness of this species cut-off is questionable when applied to short amplicons that are only representing a small part of the full 16S rRNA gene. Indeed, the selected amplicon might be evolutionary more conserved for a specific taxonomic lineage, leading to the merging of different species at the OTU level. Based on our observations we claim that integrating the differential evolutional rates of taxonomic lineages by defining a taxonomic dependent OTU cut-off score, provides a more accurate correspondence between OTUs and species. 
#####In this context, we have developed a new tool named Dynamic Cutoff (DynamiC) capable of building a dynamic cut-off table "lookup table" for clustering that is taxonomic dependent region specific. In addition, DynamiC is also able to utilize this lookup table to further cluster 16S rRNA sequencing data into a more accurate clusters of operational taxonomic units (OTUs). DynamiC is perl based and freely available.  
#Installation Requirement:
Perl, R and mothur need to be installed in order to be able to run the software, they can be installed from https://www.perl.org/, https://www.r-project.org/ and http://www.mothur.org/ respectively. (Mothur needs to be installed in the same directory)


#Syntax:
####	perl DynamiC.pl {options}
#####make sure you use an underscore (not a hyphen) to specify the options!
#####make sure you use the complete path when describing files!


There are two modes to run DynamiC, Training or testing mode, they can be specified via the "_m" option. The training mode is used to build up your lookup table that will be utilized in later on when running the testing mode. The lookup table is a table of specific clustering cutoff that depend on the taxonomic rank specified (default is taxonomic family) and the location (and length) of the used amplicon within the 16S rRNA gene. 
For instance:

    Taxonomic Family      Position(0-300)   Position(50-350)                          
    Staphylococcaceae     0.027             0.021
    Streptococcaceae      0.012             0.03
    Pseudomonadaceae      0.025             0.01

Yet. the produced lookup table would calculate appropriate cutoffs for the entire 16S rRNA gene regions (with a length and frequency dependent on the options used see below). This lookup table will be then fed to as one of the input for the testing mode, where DynamiC would utilize these cutoffs to cluster the input sequences into OTUs with potentially closer correspondence to the existing species with the sample, and lesser variability upon comparing results from different region of the 16S rRNA gene of the same sample(s). 


##Training Mode:
Used to calculate a lookup table to be used in the testing mode
to specify training mode use the following option 
#####	perl DynamiC.pl _m train

####Mandatory Options:	
	_f	Fasta file
		Fasta file of the full length 16S rRNA gene sequences. 
		this file can be downloaded from various databases such as Silva, greengene and RDP. 
		Here we applied the type straing Living Tree Project fasta file.

	_t	Two column file with the database IDs and the Taxonomic Rank.
		We recommend using the Family taxonomic rank. The file should look like this: 
		ID		Taxonomic_Rank
		X78017		Acidaminococcaceae
		AF473835	Acidaminococcaceae
		X72865		Acidaminococcaceae
		AB490811	Acidaminococcaceae
		X81037		Spirochaetaceae
		AJ0012380	Spirochaetaceae

	_r	a pre-aligned Reference dataset (that will be used for alignment)
		Here we apply mothur version of silva database
		It can be downloaded from http://www.mothur.org/wiki/Silva_reference_files

####Non-mandatory Options:
	_w	Average window size (default 500)
		This option allows specifying the length of your amplicon. for instance: 
		If MiSeq with completely overlapping reads is used with a length of 300 bases, this option is set to 300

	_s	Window frame shift size (default 50)
		This option allows specifying the resolution level of your lookup table. It can range from 1 (slowest) to 1600 (fastest).
		It should be used to easily locate that cutoff of your amplicon region within the 16S rRNA gene.
		Applying the default value (with _w 500) would calculate the cutoffs for positions ranging from (1-500), (50-550), (100-600) etc.

 
	_c	Cut-off stringent level (between 0 to 1, default= 0.1)
		The distances between sequences within the selected taxonomic rank (e.g. Family rank) would be sorted.
		0 would select the minimum value, and 1 would select the maximum value.
		For Sanger sequenced full length 16S rRNA gene, a value of 0.025 is acceptable to remove outlayers.
		For NGs sequencing data a value of 0.1 is acceptable to count for sequencing errors.

	_z	Minimum size of accepted number of sequences with the selected taxonomic rank  (default 3)
		Used to determine the minimum number of sequences Within a taxonomic rank sufficient to get reliable cut-off.
 		Otherwise, this taxonomic rank wont be included in the produced lookup table.
 		In such case reads belonging to this taxonomic rank would be clustered using the default cut-off (see testing mode)

	_u	 Upper cut-off value for OTUs clustering (default= 0.01)
 		This is used to assign a minimum cut-off in order not to have too stringent cut-off 
 		This can be important in order not two split species with various paralogous 16S rRNA gene or various strains

	_h	Lower cut-off value for OTUs clustering (default= 0.03) 
 		This is used to assign a maximum cut-off in order not to have too loose cut-off 
 		This can be important in order not to overmerge species into the same OTU.

##Testing Mode:
Used to cluster 16s rRNA gene amplicon sequencing data with our proposed lookup table obtained via training mode
To specify testing mode use the following option 
#####	perl DynamiC.pl _m test

####Mandatory Options:
	_w	Look Up table (created from the training mode)
	_f	Fasta file of your unique (dereplicated) fasta IDs and sequences within your sample.

	_n	Name file of your sample
		Obtain using mothur (unique.seqs) command
		Two column file "tab separated" 
			First column: IDs of unique sequences within the fasta file
			Second column: all IDs with sequences  matching the corresponding ID "comma separated"
		For example:
			ID1	ID1,ID5,ID6,ID7
			ID2	ID2,ID3
			ID4	ID4

	_t	Taxonomy file of your sample
		A file with the sequence ID on one coloumn and taxonomic classification on the other column
		Created using mothur (classify.seqs) command, for example:
		ID1	Bacteria(100);Deinococcus-Thermus(100);Deinococci(100);Thermales(100);Thermaceae(100);unclassified;
		ID2	Bacteria(100);Aquificae(100);Aquificae(100);Aquificales(100);Aquificaceae(100);unclassified;

	_k	location within the 16S rRNA gene
		According to the window length and frame shift (specified within the training mode of the lookup table)
		For instance:
			If your amplicon  targeting V3-4 region within the 16S rRNA gene with a length of 500 bases
			The _k should be equal 6 which corresponds to position 250 to 750 within the 16S rRNA gene
			As: 1 (1-500), 2 (50-550). 3 (100-600), 4 (150-650), 5 (200-700), 6 (250-750)
			Assuming _s 50 and _w 500 parameter and their values were used in training


####Non-mandatory Options:
		_e mothur clustering algorithm default is "average"
		Options are:
			Average Neighborhood  => average
			Nearest Neighborhood  => nearest 
			Furthest Neighborhood => furthest
##General Options:
	_p number of processors (default =1)
	_o Output path (Mandatory)

#Output Files
The DynamiC program generates different text output files distributed over two folders "Final" and " Temp". Inside each of them another folder can be found, having different output depending on the mode.

###Training Mode
####Final
Contains the final lookup table file named LookUp_table.
####Temp
#####Families folder:
Containing several fasta file each including all sequences with the same taxonomic rank inserted.
#####Distance_window_ "size":
The word "size" would change depending on the size specified via _w option
Containing a subfolder named according for each taxonomic rank within each:

	Fasta:		With sequences trimmed for each frame shift (specified by _s)
	Distance:	With the corresponding distance of the sequences within each fasta file
	Summary:	Produced via "R" with the distance distribution (minimum, median, max) distances among others

###Testing Mode
####Final
Contains the final list file ending with (.list) [commonly refered to as mapping file]
####Temp
Containing a subfolder named according to the existing taxonomic rank. within each:

	Fasta:		With sequences belonging to that taxonomic rank
	Distance:	With the corresponding distance of the sequences within each fasta file
	List:		Produced via clustering each sequencing using the assigned cutoff from the lookup table

#Testing
###Example command (training): 
	perl DynamiC.pl _m train _f /PATH/LTP.fasta _t /PATH/LTP.table _r /PATH/silva.bacteria.fasta
###Example command (testing):
	perl DynamiC.pl _m test _w /PATH/Test_Lookup _f /PATH/test.fasta _n /PATH/test.names _t /PATH/test.Tax.wang.taxonomy _k 6

The different input file (LTP.fasta, LTP.table, Test_Lookup, test.fasta, test.names, test.Tax.wang.taxonomy) are included in the software. It will produce in the output path a file containing the results i.e. two files named Test_Lookup and test.lookup.an.list. The Test dataset is a part of the data published in:

#####M. Schirmer, U. Ijaz, R. D'Amore, Neil Hall, W. Sloan, and C. Quince (2015) Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform, Nucl. Acids Res. 

#Citing
If you are going to use DynamiC, please cite it with the included software (mothur):

#####Mysara M., P. Vandamme, N. Leys, J. Raes and P. Monsieurs, 2016, Reconciliation between Operational Taxonomic Units and Species Boundaries, in preparation.
#####Schloss PD, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, et al. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology 75:7537–41.

#Contact us
For questions, bugs and suggestions, please refer to mohamed.mysara@gmail.com & pieter.monsieurs@sckcen.be

Developed by M.Mysara et al. 2016
