#! /usr/local/bin/perl
#.....................License.................................
#	 Software for constructing a lookup table for more accurate OTU clustering 
#    Copyright (C) 2015  <M.Mysara et al>

#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
	
#......................Packages Used..........................
use strict;
use Getopt::Std;
#####################
my %opts;
getopt('ftrwsczuhmlnkepo',\%opts);
my @options=('f','t','r','w','s','c','z','u','h','m','l','n','k','e','p','o'); 
foreach my $value(@options){
	for(my $a=0;$a<30;$a=$a+2){
		my $temp_Value = '_'.$value;
		if($temp_Value eq $ARGV[$a]){$opts{$value}=$ARGV[$a+1];}
	}
}
#assigning variables
my $reference=$opts{r};
my $output=$opts{o};
my $processors=$opts{p};if(!$processors){$processors=1;}
my $window=$opts{w};if(!$window){$window=500;}
my $shift=$opts{s};if(!$shift){$shift=50;}
my $cutoff=$opts{c};if(!$cutoff){$cutoff=0.10;}
my $family_size=$opts{z}-1;if(!$family_size){$family_size=2;}$family_size=$family_size*2;
my $upper_threshold=$opts{u};if(!$upper_threshold){$upper_threshold=0.01;}
my $lower_threshold=$opts{h};if(!$lower_threshold){$lower_threshold=0.03;}
my $mode=$opts{m};
my $log=$opts{o}."/All.logfile"; unlink($log);
if($mode eq "train"){
	
	if(-e $opts{r}){}else{print "Please insert reference database for aligning as silva\n ";exit;}
	if(-e $opts{t}){}else{print "Please insert two column file with the database IDs and the Taxonomic families.\n";exit;}
	if(-e $opts{f}){}else{print "Please insert fasta file of database (e.g. type straing Living Tree Project database) \n";exit;}
	if(-e $opts{o}){system"rm -r $output/Temp";mkdir("$output/Temp");system"rm -r $output/Final";mkdir("$output/Final"); }else{print "Please assign output path \n";exit;}
	#####################
	mkdir("./families");
	mkdir("$output/Temp/families");
	#preparing the R file
	open FH,">",$output."/Temp/cutoff_stat.R";print FH 'habal<-read.table("./habal.dist", header=TRUE,  sep=",")
	#summary(habal)
	habal<- as.matrix(habal)
	#sd(habal, na.rm = FALSE)

	#quantile(habal, c(0.1))
	#quantile(habal,c(0,0.025,0.05,0.1,0.25,0.5,0.75,1))
	quantile(habal,c('.$cutoff.'))
	';



	#getting the unique list of families (e.g. LTP.table)
	system "cut -f2 $opts{t} | sort | uniq | sed '/\^\\\W\*\$/d'  > $output/Temp/All.families";

	open FH,"$output/Temp/All.families";my @families=<FH>;close FH;
	#Extracting fasta for each family, named family.all.fasta

	for(my $i=0;$i<scalar(@families);$i++){
		my $family=$families[$i];
		chomp ($family);
		my $family_accnos=$family.".all.accnos";
		my $family_fasta=$family.".all.fasta";
		system "grep \"$family\" $opts{t} | cut -f1 > $output/Temp/families/$family_accnos";
		system "./mothur \"\#get.seqs(accnos=$output/Temp/families/$family_accnos,fasta=$opts{f})\">> $log";
		my $pick_fasta=$opts{f};
		$pick_fasta=~s/\.(\w+)$/.pick.$1/;
		system "cp $pick_fasta $output/Temp/families/$family_fasta";
		unlink("$pick_fasta");
	}

	#Align all of the sequence to silva with E_coli
	system "cat $output/Temp/families/*all.fasta ./E_coli.fasta > $output/Temp/Families.fasta";
	system "echo Z832051 > accnos";
	system "./mothur \"\#align.seqs(processors=$processors,fasta=$output/Temp/Families.fasta,flip=T,reference=$reference);get.seqs(accnos=accnos,fasta=current)\">> $log";
	unlink("accnos");
	system "cp $output/Temp/Families.pick.align $output/Temp/E_coli.align";

	#Getting the cutting the start and end position in correspondance to E_Coli 
	system("rm -r $output/Temp/Distance_window_$window");
	mkdir("$output/Temp/Distance_window_$window");
	for (my $o=0;$o<1500-$window;$o=$o+$shift){
		open FH,"$output/Temp/E_coli.align";my @ecoli=<FH>;close FH;
		my @ecoli_arr=split("",$ecoli[1]);
		my $start=0;
		my $end=0;
		my $tracker=0;
		for(my $i=0;$i<scalar(@ecoli_arr);$i++){
			if($ecoli_arr[$i]=~/[A,G,C,T,N,U]/){$tracker++;}
			if($tracker eq $o+1){$start=$i;}
			if($tracker eq $o+$window){$end=$i;}
		}
		system "./mothur \"\#pcr.seqs(processors=$processors,fasta=$output/Temp/Families.align,start=$start,end=$end);filter.seqs(vertical=T)\">> $log";
		system("cp $output/Temp/Families.pcr.filter.fasta $output/Temp/Distance_window_$window/$o.fasta");
	}

	#Building the OTU_table

	for(my $i=0;$i<scalar(@families);$i++){
		my $family=$families[$i];
		chomp ($family);
		my $fasta="$output/Temp/families/".$family.".all.fasta";
		open FH,$fasta;my @temp_arr=<FH>;close FH;
		if(scalar(@temp_arr)>$family_size){
					mkdir("$output/Temp/Distance_window_$window/$family");
					#for (my $o=0;$o<450-$window;$o=$o+100){
					for (my $o=0;$o<1500-$window;$o=$o+$shift){
						system "./mothur \"\#list.seqs(fasta=$fasta);get.seqs(fasta=$output/Temp/Distance_window_$window/$o.fasta,accnos=current);dist.seqs(processors=$processors,cutoff=1)\">> $log";			
						system ("cp $output/Temp/Distance_window_$window/$o.pick.fasta $output/Temp/Distance_window_$window/$family/$o.fasta");	
						system ("cp $output/Temp/Distance_window_$window/$o.pick.dist $output/Temp/Distance_window_$window/$family/$o.dist");
						#system ("echo \"Dist\" > ./habal.dist");
						system ("cut -f3 -d \" \" $output/Temp/Distance_window_$window/$o.pick.dist >> ./habal.dist");
						if(-e "./habal.dist"){    
							system "R CMD BATCH $output/Temp/cutoff_stat.R";#habal.R
							system "cp ./cutoff_stat.Rout $output/Temp/Distance_window_$window/$family/$o._summary";
							unlink("./cutoff_stat.Rout");
							unlink("./habal.dist");
						}
					}

					
			
			
		} 
	}
	system "rm *logfile";
	#extracting the lookuptable
	my @CutOff;


	for(my $i=0;$i<scalar(@families);$i++){
			my $family=$families[$i];
			chomp ($family);
		my $fasta="$output/Temp/families/".$family.".all.fasta";
		open FH,$fasta;my @temp_arr=<FH>;close FH;
			if(scalar(@temp_arr)>$family_size){
			$CutOff[$i]=$family."\t";

			mkdir("$output/Temp/Distance_window_$window/$family");
					#for (my $o=0;$o<450-$window;$o=$o+100){
					for (my $o=0;$o<1500-$window;$o=$o+$shift){
				open FH,"$output/Temp/Distance_window_$window/$family/$o._summary";my @sum=<FH>;close FH;
				$sum[30]=~/([\d]+\.?[\d]*)/;
				$sum[30]=$1;
				chomp($sum[30]);
				if($sum[30]<$upper_threshold){$sum[30]=$upper_threshold;}
				if($sum[30]>$lower_threshold){$sum[30]=$lower_threshold;}
				$CutOff[$i]=$CutOff[$i].$sum[30]."\t";chomp($CutOff[$i]);
		
			}	
		}
	}
	open FH,">","$output/Final/LookUp_table";print FH join("\n", @CutOff);close FH;
	system ("rm $output/Temp/Distance_window_$window/*fasta");
	system ("rm $output/Temp/Distance_window_$window/*dist");
}
elsif($mode eq "test"){
#testing on real dataset

#defining variables
	my $window=$opts{k};if(!$window){$window=1;}
	my $lookup=$opts{w};if(!$lookup){print "Please insert a valid Lookuptable, or train one using the training mode\n"; exit;}
	my $fasta=$opts{f};if(!$fasta){print "Please insert a valid fasta file\n";exit;}
	my $names=$opts{n};if(!$names){print "Please insert a valid names file\n";exit;}
	my $tax=$opts{t};if(!$tax){print "Please insert a valid Taxonomy file\n";exit;}
	my $method=$opts{e};if(!$method){$method="average";}
	if(-e $opts{o}){system"rm -r $output/Temp";mkdir("$output/Temp");system"rm -r $output/Final";mkdir("$output/Final"); }else{print "Please assign output path \n";exit;}
	
	$method=~/^(\w)/;my $clst_method=$1."n";
	my $Final_list=$opts{f};$Final_list=~s/\.\w+$/.lookup.$clst_method.list/;
	my $F_list="";
	my $F_count=0;
	
#Opening the output list file

	open OUT,">","$output/Final/".$Final_list;
#generating random ID for running
	my $i=rand();
	mkdir("./Data_results/$i");



	system "cut -f5 -d \";\" $tax | cut -f 1 -d \"\(\" | sort | uniq | sed '/\^\\\W\*\$/d'  > $output/Temp/All.families";

	open FH,"$output/Temp/All.families";my @families=<FH>;close FH;
	
	for(my $j=0;$j<scalar(@families);$j++){
		my $family=$families[$j];
		chomp $family;
		mkdir("$output/Temp/$family");
		unlink("$output/Temp/family.accnos");
		system "cut -f1,5 -d \";\" $tax | cut -f 1-3 -d \"\(\" | grep \";$family\" | cut -f1 > $output/Temp/family.accnos";

		open FHH,"$output/Temp/family.accnos";my @FM=<FHH>;close FHH;
		if(scalar(@FM)<2){
#cant not calculate distance
			chomp($FM[0]);
			my $temp_id=$FM[0];
			system "grep \"$temp_id\" $names | cut -f2 > $output/Temp/temp_names";
			open TN, "$output/Temp/temp_names";my @Temp_NM=<TN>;close TN;
			chomp($Temp_NM[0]);
			#my $c = () = $Temp_NM[0] =~ /\,/g;  # $c is now 3
			$F_count=$F_count+1;
			$F_list=$F_list.$Temp_NM[0]."\t";
			
			system("mv $output/Temp/temp_names $output/Temp/$family/$family.names");
		}
		else{
		
		
			my $list=$fasta;
			$list=~s/fasta/pick.$clst_method.list/;
			unlink($list);
			system "./mothur \"\#get.seqs(fasta=$fasta,name=$names,accnos=$output/Temp/family.accnos);dist.seqs(processors=$processors,cutoff=0.1);cluster(method=$method,name=current,precision=1000,column=current,hard=f)\">> $log";

#Getting the cut-off
			system("grep \"$family\" $lookup > $output/Temp/LOOKUP.temp");
			open FH,"$output/Temp/LOOKUP.temp";my @LOOKUP_=<FH>;close FH;  
			my @LOOKUP=split("\t",$LOOKUP_[0]);
			my $cutoff=$LOOKUP[$window];
            if(!$cutoff){$cutoff=0.03;}
			print $family.":	".$cutoff; if($cutoff<0.01){$cutoff=0.01;}elsif($cutoff>0.04){$cutoff=0.04;}else{}
			print "\t".$cutoff."\n";		

#Getting the closest cluster to the cut-off selected
			my $cutoff_=$cutoff;
			system "tail -n +2 $list | cut -f1 > $output/Temp/cutoff_ls";
			open Fhh, "$output/Temp/cutoff_ls";my @cutoff_ls=<Fhh>;close Fhh;
			for(my $k=0;$k<scalar(@cutoff_ls);$k++){
				chomp($cutoff_ls[$k]);
				if($cutoff_ls[$k]<=$cutoff_){$cutoff=$cutoff_ls[$k]};
			}
#Adding to the Final list file
			system "grep \"$cutoff\" $list > $output/Temp/cutoff_ls";
			open Fhh, "$output/Temp/cutoff_ls";my @cutoff_ls=<Fhh>;close Fhh;
			chomp($cutoff_ls[0]);
			$cutoff_ls[0]=~/$cutoff\t(\d+)\t([\w\W]+)$/;
			$F_count=$F_count+$1;$F_list=$F_list.$2;
			my $files=$list;$files=~s/\.\w\w\.list/\*/;
			system("mv $files $output/Temp/$family/");
		}
		

	}
	my $header="label	numOtus	";	
	for(my $k=1;$k<$F_count+1;$k++){
		if($k<10){$header=$header."Otu0000".$k."\t";}
		elsif($k<100){$header=$header."Otu000".$k."\t";}
		elsif($k<1000){$header=$header."Otu00".$k."\t";}
		elsif($k<10000){$header=$header."Otu0".$k."\t";}
		elsif($k<100000){$header=$header."Otu".$k."\t";}
		else{$header=$header."Otu".$k;}		
	}
	$header=$header."\n0.03	$F_count	";
	print OUT $header.$F_list;
	close OUT;
	system("rm *logfile");
	
}
else{
	print "Please select either the training or testing mode using _m train/test";
	usage();
}
sub usage{
print
"	
	||||||||||||||||||||||||||||||||||||||||||||||||
	||              Welcome To DynamiC	    	  ||
	|| Software for constructing a lookup table   ||
	||		for more accurate OTU clustering	  ||
	||    Copyright (C) 2015  <M.Mysara et al>    ||
	||||||||||||||||||||||||||||||||||||||||||||||||

 DynamiC version 1, Copyright (C) 2015, M.Mysara et al
 DynamiC comes with ABSOLUTELY NO WARRANTY.
 This is free software, and you are welcome to redistribute it under
 certain conditions; please refer to \'COPYING\' for details.;

 The software also includes \"mothur\" also under GNU Copyright
 
 
Command Syntax:
(perl) DynamiC.pl {options} 	

 There are two modes to run DynamiC, Training or testing mode:
	_m mode (either train or test).

Training Mode:
Used to calculate a lookup table to be used in the testing mode

	#Mandatory Options:
	 _m train  
	 _f Fasta file of database (e.g. type straing Living Tree Project database) 
	 _t Two column file with the database IDs and the Taxonomic families.
	 _r Reference database for aligning as silva.
	 
	 #Non-mandatory Options:
	 _w average window size (resembling the size of the amplicons, default= 500)
	 _s window fram shift size (default= 50)
	 _c cut-off stringent level (between 0 to 1, default= 10)
	 _z Minimum size of accepted families with (default 3)
	 _u Upper cut-off value for OTUs clustering (default= 0.01)
	 _h Lower cut-off value for OTUs clustering (default= 0.03) 
	 
Testing Mode:
Used to calculate a lookup table to be used in the testing mode

	#Mandatory Options:
	 _m test
	 _w Look Up table (created from the training mode)
	 _f Fasta file of your sample
	 _n name file of your sample
	 _t Taxonomy file of your sample (created using mothur classify.seqs command)
	 _k location within the 16S rRNA gene
	 
	 #Non-mandatory Options:
	 _e mothur clustering algorithm (average, nearest or furthest, default is average)

#General Options:
	_p number of processors (default =1)
	_o Output path (MANDATORY)
 For Queries about the installaion, kindly refer to \'README.txt\'
 For Queries about the Copy rights, kindly refer to \'COPYING.txt\'

CITING [please cite the included software (Mothur)]:
(M.Mysara et al, 2015),(PD.Schloss, et al. 2009).

";
}
