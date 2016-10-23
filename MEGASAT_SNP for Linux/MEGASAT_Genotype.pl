### Separates files by primer
### TO DO:
###   Trim off primers (wobble = 2)
###   Filter weirdos (might need to process if coverage is low)
###   Look at length distribution, find modes

use strict;
use warnings;
use Parallel::ForkManager;
use SNP;

# The function to match primer sequences to reads, with a tolerance for mismatches (see definition at the bottom of this file)
sub FUZZYMATCH ($$$);
# The function to find the starting position of repeat array
sub CONTIREPEAT($$);
# The function to find fuzzy longest common substring based on hamming distance
sub lc_substr($$);
### The function to infer SNPs from sequencing errors
sub call_multinomial_snp($$$$);
# Hash table that preserve the ratio information for each locus
my %Hash_ratios = (); 
# Hash table that contains all the corresponding information for each locus
my %primerSeqs = ();
my %length_snp = ();
# array that stores the version number for the same length allele 
my %versionArray =();
# Array that contains all the loci name
my @name_MS = ();
my @MS= ();
my $suffix="-b";
my $version="_V";
my $heterozygote_limit = -3.84;
my $homozygote_limit = 3.84;
# A tab-separated primer file with column 1 = Locus name, column 2 = Forward primer, column 3 = Reverse primer, column 4 = 3' flank
# column 5 = 5' flank, column 6 = repeat array
# A header line is required
# The maximum # of mismatches as the second argument
# The minimum depth threshold as the third argument
# The directory of data set folder that contains all the input sequence read files (fastq file) as the fourth argument
# The directory to save the output folder as the fifth argument
my ($inputPrimers,$mismatches,$m_depth,$max_processors,$dataset,$saveDir,$option);
if(!(scalar @ARGV ==7 || scalar @ARGV ==6)){	
      die "Missing command line arguments!\n
#######################################################################\n
Need the directory of primer file as the first command-line argument\n
Need the maximum number of mismatches as the second command-line argument\n
Need the minimum depth threshold as the third command-line argument\n
Need the directory of data set folder that contains all the input sequence read files (fastq or fasta) as the fourth command-line argument\n
Need the directory to save the output folder as the fifth command-line argument\n"; 
} else{
	if(scalar @ARGV ==6){
		($inputPrimers,$mismatches,$m_depth,$max_processors,$dataset,$saveDir)= @ARGV;
	} elsif(scalar @ARGV ==7 && $ARGV[6]==1){
		($inputPrimers,$mismatches,$m_depth,$max_processors,$dataset,$saveDir,$option)= @ARGV;
	} else{
		die "You should type 1 in the last command line argument for not generating split files!";
	}	
	if(!(($mismatches =~ /^\d+$/) && ($m_depth =~ /^\d+$/))){
		die "The maximum number of mismatches and the minimum depth threshold should be integers";
	} elsif(! -d $saveDir){
		die "The directory to save your output is wrong!";
	}
}
# get the name of data set
my $nameofdataset = substr($dataset,rindex($dataset,'/')+1);
my $path="$saveDir/Output_SNP_$nameofdataset";
my ($S_path,$G_path,$D_path,$T_path);
# check if the folder already exists in the saving directory
if(-d $path){
	die "Output_$nameofdataset exists in $saveDir! Please save in another directory or change the existing folder name! ";
}
# create the output folder
else{
	mkdir $path, 0755;
}
if(scalar @ARGV ==6){
	$S_path = "$path/Sorted";
	mkdir $S_path, 0755;
	$D_path = "$path/Discarded";
	mkdir $D_path, 0755;
	$T_path = "$path/Trimmed";
	mkdir $T_path, 0755;
}

# disable perl output buffering
select((select(STDOUT), $|=1)[0]);
### Read the primers txt file and associate loci names with primers and flank information
open (IN,$inputPrimers) or die "cannot open primers txt file!";
my $toss = <IN>;
while (<IN>) {
	chomp;
	my @line = split /\t/,$_;
	my @Default_ratios = (0.15,0.4,0.7,0.6,0.8,0.2);
	## when the case that the primer text file has 7 columns (including ratios group)
	if (scalar @line == 7) {
		  $line[6] =~ s/\r$//;
		  ## split the ratio group by comma and save into an array
		  my @G_ratios = split(',',$line[6]);
		  @{$primerSeqs{$line[0]}} = @line[1..5];
		  ## if it does not have six ratios or the first five elements are not string, we will print the format of primers text file is wrong
		  if(scalar @G_ratios != 6 || (grep {!/^[a-zA-Z]+$/} @line[1..5])){
			  die "The format of primers txt file is wrong!There should be six ratios or the primers and flanking region should only contain string like AGCT";
		  } else{
			  ##convert the lower case letter to upper case letter in case some primers are written in lower case
			  map {$_ = uc($_)} @{$primerSeqs{$line[0]}};
   		      ### get the index list of space in the array
   		      my (@index_space) = grep {$G_ratios[$_] eq ' '}	0..@G_ratios-1;
			  ### get the index of non space in the array
			  my (@left) = grep {$G_ratios[$_] ne ' '} 0..@G_ratios-1;
			  ### check if the non space element is a ratio that is larger than 0 and smaller than 1
			  if(grep {!/^[0]+(\.[0-9][0-9]?)/} @G_ratios[@left]){
				  die "The format of primers txt file is wrong!The six ratios should be a number between 0 and 1";					   
			  } else{
   		   ## replace space in the array with default ratios
   		          foreach my $i(@index_space){$G_ratios[$i] = $Default_ratios[$i]; }
                  @{$Hash_ratios{$line[0]}} = @G_ratios[0..5];	
		      }
		  }	
	} elsif(scalar @line == 6){
		  $line[5] =~ s/\r$//;
		  if(grep {!/^[a-zA-Z]+$/} @line[1..5]){
			  die "The format of primers txt file is wrong!The primers and flanking region should only contain string like AGCT";
		  }	else{
			  @{$primerSeqs{$line[0]}} = @line[1..5];
			  ##convert the lower case letter to upper case letter in case some primers are written in lower case
			  map {$_ = uc($_)} @{$primerSeqs{$line[0]}};	
			  @{$Hash_ratios{$line[0]}} = @Default_ratios[0..5];
		  }    		  
	} else{
	      die "The format of primers txt file is wrong!The primers txt file does not contain all the primer and flanking information for some loci";
	}
}
close (IN);

### put all the loci name into the array @MS
foreach my $MSname (sort keys %primerSeqs) {
	push(@name_MS, $MSname);
    push (@MS, $MSname);
	push(@MS,$MSname.$suffix);
	}

### produce a txt file that has the minimum depth threshold and all the ratios for each locus
open (OUT, ">$path/Ratios_Threshold_$nameofdataset.txt") or die "cannot generate text file that contains rations and minimum depth threshold";
print OUT "minimum depth threshold,$m_depth\n";
print OUT "$_,@{$Hash_ratios{$_}}\n" for sort keys %Hash_ratios;
close (OUT);
	
### produce a txt file that has all the genotype information for all the individuals and loci
open (OUT1, ">$path/Genotype.txt") or die "cannot generate Genotype.txt file";
print OUT1 "Sample_idx1_idx2\t", join "\t",@MS,"\n";

### produce a txt file that has the number of discarded sequences for all the individuals and loci
open (OUT2, ">$path/Number_Discarded.txt") or die "cannot generate Number_Discarded.txt file";
print OUT2 "Sample_idx1_idx2\t", join "\t",@name_MS,"\n";

### produce a txt file that has the number of sorted sequences for all the individuals and loci
open (OUT3, ">$path/Number_Sorted.txt") or die "cannot generate Number_Sorted.txt file";
print OUT3 "Sample_idx1_idx2\t", join "\t",@name_MS,"\n";

### produce a txt file that has the number of trimmed sequences for all the individuals and loci
open (OUT4, ">$path/Number_Trimmed.txt") or die "cannot generate Number_Trimmed.txt file";
print OUT4 "Sample_idx1_idx2\t", join "\t",@name_MS,"\n";

# open the data set folder and start reading fastq file
opendir(DIR, "$dataset" ) or die "cannot open the data set folder!";
my @files = grep{/\.fastq$/ || /\.fasta$/} readdir(DIR);
closedir DIR;
## using fork to do the parallel processing
my $fork= new Parallel::ForkManager($max_processors);
select((select(STDOUT), $|=1)[0]);
my $index_file = 0;
if(scalar @files == 0){
    die "No fastq file in the data set";
}
DATA_LOOP:
foreach my $file(@files){
      $index_file++;
      my %seqsMapped = ();    ## hash table that contains sequences that are trimmed off forward primers
      my %seqsTrimmed = ();    ## hash table that contains sequences that are trimmed off forward primers and reverse primers
	  my %seqsDiscarded = ();  ## hash table that contains sequences that are discarded
      $fork->start and next DATA_LOOP;
	  open (IN, "./$dataset/$file") or die "cannot open $file!";



# For each line, we will look for matches in the following order:
#   (1) Is there a fuzzy match for any forward primer?
#   If not, throw out this line.
#   If so, delete the primer from the current line and add it to %seqsMapped. Then:
#   (2) Is there a fuzzy match for the corresponding reverse primer?
#   If so, delete the reverse primer from the current line, and 
#   (3) check the trimmed sequence to see if it contains a fuzzy match for 5'flank & 3'flank and a long continuous repeat array.
#    If so, add the trimmed sequence to %seqsTrimmed. 
#   If it cannot meet requirement (2), maybe due to the large allele, full reverse primer cannot be shown in the sequence. So
#   (4) check this line to see if it contains a fuzzy match for 5'flank $ 3'flank and the ending point of repeat array is larger
#   than the length of this line minus the length of reverse primer
#   If so, find the ending point of 3'flank and delete the remaining sequence started from the ending point,add the trimmed sequence to %seqsTrimmed.
#	(5) If not, check this line to see if it contains a fuzzy match for 5'flank and the ending point of repeat array is larger 
#   than the length of this line minus the length of 3'flank
#   If so, this line represents a large allele, add the missing 3'flank part to this line. And add the new line to %seqsTrimmed.
#	(6) if this line contains a fuzzy match for 5'flank and the ending point of repeat array is very close to the end of this line 
#   This line represents a large allele, we use 200+ the length of 3'flank to represent this allele.

while (<IN>) {
	chomp;
	my $thisLine = $_;
	if ($thisLine =~ /^[AGCTagct]+$/) {
		$thisLine = uc $thisLine;
		## value of large alleles
		my $large_allele = (length $thisLine) + 100;
		my $foundForward = undef;
		my $foundReverse = undef;		
		FINDPRIMER: foreach my $checkPrimer (keys %primerSeqs) {
			if (!defined( )) {
				die "$checkPrimer";
			 }
			my $isMatched = FUZZYMATCH($primerSeqs{$checkPrimer}[0],$thisLine,$mismatches);	
			# Is there a fuzzy forward primer match?
			if ($isMatched >= 0) {
				$foundForward = $checkPrimer;
				substr($thisLine,0,$isMatched + length $primerSeqs{$checkPrimer}[0]) = "";
				push (@{$seqsMapped{$foundForward}},$thisLine);
				last FINDPRIMER;	
			 }
		}
		if (defined($foundForward)) {
		    my ($LRIndex,$NumRepeat) = CONTIREPEAT($thisLine,$primerSeqs{$foundForward}[4]);
			my $revPrimer = $primerSeqs{$foundForward}[1];
			my $fiveFlank = ($primerSeqs{$foundForward}[2] eq 'X')?'':$primerSeqs{$foundForward}[2];
			my $threeFlank = ($primerSeqs{$foundForward}[3] eq 'X')?'A':$primerSeqs{$foundForward}[3];
			my $isRevMatched = FUZZYMATCH($revPrimer,$thisLine,$mismatches);
			my $isTfMatched = FUZZYMATCH($threeFlank,$thisLine,$mismatches);
			my $endRepeat = $LRIndex+(length $primerSeqs{$foundForward}[4]);
			my $errate_fiveflank=((length $fiveFlank)!=0)?((length $fiveFlank)-(substr($thisLine,0,length $fiveFlank) ^ $fiveFlank)=~ tr[\0][\0])/(length $fiveFlank):0;
			if (!defined($revPrimer)) {
				die $foundForward;
			}
			# Is there a fuzzy reverse primer match?
		    if ($isRevMatched > 0) {
			    $foundReverse = $revPrimer;
				my $leftPart=substr($thisLine,$isRevMatched);
				substr($thisLine,$isRevMatched) = "";
				# check the 5'flank & 3'flank and repeat array
				if($errate_fiveflank<=0.2 &&  $LRIndex>=0 && ((((length $threeFlank)>1 &&(length $threeFlank)<=5) ? substr($thisLine,-(length $threeFlank)) eq $threeFlank : $isTfMatched>=0) || $foundForward eq "BF-372")){
				      push @{$seqsTrimmed{$foundForward}},$thisLine;
				 } else{
				      $thisLine .=$leftPart;
					  push @{$seqsDiscarded{$foundForward}},$thisLine;
				    }
			 } elsif($errate_fiveflank<=0.2 && ((length $threeFlank)<=3 ? $thisLine =~ /$threeFlank/ : $isTfMatched>=0) && $LRIndex>=0 && $endRepeat <= (length $thisLine)-(length $threeFlank)) {
				   if(($endRepeat < (length $thisLine)-(length $threeFlank)-(length $revPrimer)) || ($endRepeat >= (length $thisLine)-(length $threeFlank)-4 && $foundForward ne "BF-456")){
				       if((length $threeFlank)==1){
				             substr($thisLine, $endRepeat) = "";
					         push @{$seqsTrimmed{$foundForward}},$thisLine;
				         } elsif((length $threeFlank)>20){
				             substr($thisLine,$isTfMatched+length $threeFlank)="";
				             push @{$seqsTrimmed{$foundForward}},$thisLine;
				         } else{
					         substr($thisLine, $endRepeat+ (length $threeFlank)) = "";
				             push @{$seqsTrimmed{$foundForward}},$thisLine;
				            } 
				   } else{
				       substr($thisLine, -lc_substr(substr($thisLine,-(length $revPrimer)),$revPrimer)) = "";
				       push @{$seqsTrimmed{$foundForward}},$thisLine;
				        }
			} elsif($errate_fiveflank<=0.2 && $LRIndex>=0 && $endRepeat > (length $thisLine)-(length $threeFlank)){
				      if((length $fiveFlank)>0 || substr($thisLine,0,length $primerSeqs{$foundForward}[4]) eq $primerSeqs{$foundForward}[4]){
					      if($endRepeat==(length $thisLine)|| (((length $thisLine)-$endRepeat)-(substr($thisLine,$endRepeat) ^ substr($primerSeqs{$foundForward}[4],0,(length $thisLine)-$endRepeat))=~ tr[\0][\0])<=1){
						         my $startRepeat = $LRIndex-($NumRepeat-1)*(length $primerSeqs{$foundForward}[4]);
						         my $addString="X" x ($large_allele-(length $thisLine)+$startRepeat);
							     $thisLine .=$addString;
						         push @{$seqsTrimmed{$foundForward}},$thisLine;
							} elsif((((length $thisLine)-$endRepeat)- (substr($thisLine,$endRepeat) ^ substr($threeFlank,0,(length $thisLine)-$endRepeat))=~ tr[\0][\0])<=2){
				                my $leftTF=substr($threeFlank,(length $thisLine)-$endRepeat);
				                $thisLine .=$leftTF;
				                push @{$seqsTrimmed{$foundForward}},$thisLine;
				            } else{
							    push @{$seqsDiscarded{$foundForward}},$thisLine;
							}
				        } else{
						    push @{$seqsDiscarded{$foundForward}},$thisLine;
						}
			} else{
			        push @{$seqsDiscarded{$foundForward}},$thisLine;
			}
		}
	}
}
close (IN);

# Create the prefixes for the output files
my $prefix = $file;
$prefix =~ s/_.*//;
# Print the full sequences that were associated with each primer pair
my %num_Sorted = ();
if(scalar @ARGV == 6){
	foreach my $outFile (keys %seqsMapped) {
	    open (OUT, ">$S_path/Sorted_${prefix}_$outFile.split");
		foreach my $outSeq (@{$seqsMapped{$outFile}}) {
			print OUT "$outSeq\n";
			}
		close (OUT);
	}
} else{
	foreach my $outFile (keys %seqsMapped) {
		foreach my $outSeq (@{$seqsMapped{$outFile}}) {
			++$num_Sorted{$outFile};
	    }
	}
} 

## print the discarded sequences that trimmed off forward primers but don't have reverse primers and flank
## hash table that counts the number of discarded sequence for different locus
my %num_Discarded = ();
if(scalar @ARGV == 6){
	foreach my $DiscardFile (keys %seqsDiscarded) {
	    open (OUT, ">$D_path/Discarded_${prefix}_$DiscardFile.split");
		foreach my $outSeq (@{$seqsDiscarded{$DiscardFile}}) {
			print OUT "$outSeq\n";
			++$num_Discarded{$DiscardFile};
			}
		close (OUT);
	}
} else{
	foreach my $DiscardFile (keys %seqsDiscarded) {
		foreach my $outSeq (@{$seqsDiscarded{$DiscardFile}}) {
			++$num_Discarded{$DiscardFile};
			}
	}
}

# Print the trimmed sequences associated with each primer pair, *and* count up the occurrences of different lengths
my %number = ();
# hash table that count the occurrences of different lengths for different locus
my %lengths = ();
# array that preserve the length range with version number
my %lenRange_V = ();
# array that preserve the length range in ascending
my @MSlens = ();
# array that stores the identical length sequences for each locus 
my %lengthArray =();
## hash table that counts the number of trimmed sequences for each locus in individuals
my %num_Trimmed = ();

# Print the trimmed sequences associated with each primer pair
if(scalar @ARGV == 6){
	foreach my $outPrimer (keys %seqsTrimmed) {
		open (OUT, ">$T_path/Trimmed_${prefix}_$outPrimer.split");
		foreach my $outSeq (@{$seqsTrimmed{$outPrimer}}) {
			++$num_Trimmed{$outPrimer};
			print OUT "$outSeq\n";
			my $sLen = length $outSeq;
			##++$another_lengths{$outPrimer}{$outSeq};
			push @{$lengthArray{$outPrimer}{$sLen}}, $outSeq;
		}
		close (OUT);
	}
} else{
	foreach my $outPrimer (keys %seqsTrimmed) {
		foreach my $outSeq (@{$seqsTrimmed{$outPrimer}}) {
			++$num_Trimmed{$outPrimer};
			my $sLen = length $outSeq;
			##++$another_lengths{$outPrimer}{$outSeq};
			push @{$lengthArray{$outPrimer}{$sLen}}, $outSeq;
		}
	}
}

###call function "call_multinomial_snp" to check which is SNP

foreach my $outPrimer (keys %lengthArray){ 
	foreach my $sLen (keys %{$lengthArray{$outPrimer}}){
		my @index;
		if(!exists $versionArray{$outPrimer}{$sLen}){
			$versionArray{$outPrimer}{$sLen}=0;
		}
		my %array_snp = ();
		my $size = scalar @{$lengthArray{$outPrimer}{$sLen}} -1;
		foreach my $i (0..$size){
			$array_snp{$i} = substr($lengthArray{$outPrimer}{$sLen}[$i],0,1);
		}
		for (my $col = 0; $col< $sLen; $col++){
			my @snps = ();
			my %nuc = ();
			my @array = ();
			my @search = ();
			@index = ();
			$nuc {"A"} = 0;
			$nuc {"G"} = 0;
			$nuc {"C"} = 0;
			$nuc {"T"} = 0;
			foreach my $i (0..$size){
				push (@array,substr($lengthArray{$outPrimer}{$sLen}[$i],$col,1));
			}
			foreach (@array) {
				if(exists $nuc{$_}){
				   $nuc{$_}++;
			    }
			}
			call_multinomial_snp (\@snps,$col,\%nuc,1);
			#my @sort_depth = sort { $b <=> $a } values %nuc;
			 if($snps[0]->{type} eq "snp_type_het" ){
				push(@search,$snps[0]->{rank_1});
				push(@search,$snps[0]->{rank_2});
				foreach my $search (@search){
					push (@index, grep {$array_snp{$_} eq $search} keys %array_snp);
				}
				@index = sort {$a <=> $b} @index;
				%array_snp = ();
				foreach my $i (@index){
					$array_snp{$i} = substr($lengthArray{$outPrimer}{$sLen}[$i],($col< $sLen-1)?$col+1:$sLen-1,1);
				}
			 } else{
 				push(@search,$snps[0]->{rank_1});
 				foreach my $search (@search){
 					@index = grep {$array_snp{$_} eq $search} keys %array_snp;
 				}
 				@index = sort {$a <=> $b} @index;
				%array_snp = ();
 				foreach my $i (@index){
 					$array_snp{$i} = substr($lengthArray{$outPrimer}{$sLen}[$i],($col< $sLen-1)?$col+1:$sLen-1,1);
 				}
			 } 
		 }
			foreach my $i (@index){
				my $key = $lengthArray{$outPrimer}{$sLen}[$i];
				if(exists $length_snp{$outPrimer}{$key}){
					$lengths{$outPrimer}{$length_snp{$outPrimer}{$key}}++;
					$lenRange_V{$length_snp{$outPrimer}{$key}} = 1;
			    } else{
			    	$versionArray{$outPrimer}{$sLen}++;
					$length_snp{$outPrimer}{$key} = $sLen.$version.$versionArray{$outPrimer}{$sLen};
					$lengths{$outPrimer}{$length_snp{$outPrimer}{$key}} = 1;
					$lenRange_V{$length_snp{$outPrimer}{$key}} = 1;
			    }
			}		
	}
	### produce a txtx file that has all the alleles with SNP in it and the corresponding version number
	if($index_file == scalar @files){
		open (OUT5, ">$path/Sequence_Version_$outPrimer.txt") or die "cannot generate Sequence_Version.txt file";
		print OUT5 "$length_snp{$outPrimer}{$_}\t$_\n" foreach (sort {(($length_snp{$outPrimer}{$a}=~ m/(\d+)/g)[0])+($length_snp{$outPrimer}{$a}=~ m/(\d+)/g)[1]*(10**(-3)) 
		<=> (($length_snp{$outPrimer}{$b}=~ m/(\d+)/g)[0])+($length_snp{$outPrimer}{$b}=~ m/(\d+)/g)[1]*(10**(-3))} keys %{$length_snp{$outPrimer}});
		close (OUT5);
	}	
}

# Print out a table showing the length distribution of each microsatellite locus
print OUT1 "$prefix";
print OUT2 "$prefix";

@MSlens = sort {(($a=~ m/(\d+)/g)[0])+($a=~ m/(\d+)/g)[1]*(10**(-3)) <=> (($b=~ m/(\d+)/g)[0])+($b=~ m/(\d+)/g)[1]*(10**(-3))} keys %lenRange_V;
# Print a txt file that shows the length distribution and genotype of each locus
open (OUT, ">$path/Genotype_${prefix}.txt");
print OUT "Microsatellite,", join ",",@MSlens,"sum,scores\n";

foreach my $outMS (sort keys %primerSeqs) {
	print OUT "$outMS";
	my $sum = 0;
	my @Plength= ();
	my @sorted= ();
	foreach my $i (@MSlens) {
		if (defined($lengths{$outMS}{$i})) {
			print OUT ",$lengths{$outMS}{$i}";
			$sum += $lengths{$outMS}{$i};
		} else {
			print OUT ",0";
		}
	}
	print OUT ",$sum";
  if(defined($lengths{$outMS}) && scalar keys %{$lengths{$outMS}} >=2){
	 @sorted = sort {$lengths{$outMS}{$b} <=> $lengths{$outMS}{$a}} keys %{$lengths{$outMS}};
	 my $sum=$lengths{$outMS}{$sorted[0]}+$lengths{$outMS}{$sorted[1]};
	if($sum>=$m_depth){
	    if(($sorted[0]=~ m/(\d+)/)[0] <($sorted[1]=~ m/(\d+)/)[0]){
	         my $ratio= $lengths{$outMS}{$sorted[1]}/$lengths{$outMS}{$sorted[0]};
	           if($ratio>=$Hash_ratios{$outMS}[0]){
	               if((defined($sorted[2]) && ($sorted[2]=~ m/(\d+)/)[0] > ($sorted[1]=~ m/(\d+)/)[0]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[1]) || (defined($sorted[3]) && ($sorted[3]=~ m/(\d+)/)[0] > ($sorted[1]=~ m/(\d+)/)[0]&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[1])){
		               if(($sorted[2]=~ m/(\d+)/)[0]-($sorted[1]=~ m/(\d+)/)[0]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[2]){
			               push (@Plength, ($sorted[0],$sorted[2]));
				        } elsif(defined($sorted[3]) && ($sorted[3]=~ m/(\d+)/)[0]-($sorted[1]=~ m/(\d+)/)[0]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[2]){
				           push (@Plength, ($sorted[0],$sorted[3]));
				        } else{
				            push (@Plength, ('Unscorable','Unscorable'));
				        }
		            } else{
	                   push (@Plength, ($sorted[0],$sorted[1]));
			           }
	            } else{
	                   push (@Plength, ($sorted[0],$sorted[0]));
	            }
	    }
	    if(($sorted[0]=~ m/(\d+)/)[0] > ($sorted[1]=~ m/(\d+)/)[0]){
	       if((($sorted[0]=~ m/(\d+)/)[0]-($sorted[1]=~ m/(\d+)/)[0]>=3&&$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[3]*$lengths{$outMS}{$sorted[0]})||(($sorted[0]=~ m/(\d+)/)[0]-($sorted[1]=~ m/(\d+)/)[0]<=2&&$lengths{$outMS}{$sorted[1]}>=$Hash_ratios{$outMS}[4]*$lengths{$outMS}{$sorted[0]})){
	           if(defined($sorted[2]) && ($sorted[2]=~ m/(\d+)/)[0] > ($sorted[0]=~ m/(\d+)/)[0]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[0]}>=$Hash_ratios{$outMS}[5]){
		           push (@Plength, ('Unscorable','Unscorable'));
		        } else{
		           push (@Plength, ($sorted[1],$sorted[0]));
		        }
	        } elsif(defined($sorted[2]) && ($sorted[2]=~ m/(\d+)/)[0] > ($sorted[0]=~ m/(\d+)/)[0]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[0]}>=$Hash_ratios{$outMS}[5]){
	               if(defined($sorted[3]) && ($sorted[3]=~ m/(\d+)/)[0]-($sorted[2]=~ m/(\d+)/)[0]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[2]}>=$Hash_ratios{$outMS}[2]){
	                   push (@Plength, ($sorted[0],$sorted[3]));
	                } else{
	                   push (@Plength, ($sorted[0],$sorted[2]));
		            }
	        } else{
	               push (@Plength, ($sorted[0],$sorted[0]));
	        }
	    }
		if(($sorted[0]=~ m/(\d+)/)[0] == ($sorted[1]=~ m/(\d+)/)[0]){
			push (@Plength, ($sorted[0],$sorted[1]));
		}
    } else{
        push (@Plength, (0,0));
      }
 } elsif(scalar keys %{$lengths{$outMS}}==1){
      @sorted = sort {$lengths{$outMS}{$b} <=> $lengths{$outMS}{$a}} keys %{$lengths{$outMS}};
      if($lengths{$outMS}{$sorted[0]}>=$m_depth){
         push (@Plength, ($sorted[0],$sorted[0]));
        } else{
             push (@Plength, (0,0));
        }
} else{
   push (@Plength, ('X','X'));
}
if(defined ($num_Discarded{$outMS})){
   print OUT2 "\t$num_Discarded{$outMS}";
} else{
   print OUT2 "\tX";
}
## print the number of sorted sequences to text file
if(defined ($num_Sorted{$outMS})){
     print OUT3 "\t$num_Sorted{$outMS}";
} else{
     print OUT3 "\tX";
}
## print the number of trimmed sequences to text file
if(defined ($num_Trimmed{$outMS})){
     print OUT4 "\t$num_Trimmed{$outMS}";
} else{
     print OUT4 "\tX";
}
print OUT1 "\t$Plength[0]\t$Plength[1]";
print OUT ",@Plength\n"; 
}
  print OUT1 "\n";
  print OUT2 "\n";
  print OUT3 "\n";
  print OUT4 "\n";
  close (OUT);
  ## do the exit in the child process
  $index_file = $$;
  print STDOUT "Running the child process $index_file\n";
  $fork->finish;
}
$fork->wait_all_children;
print STDOUT "Completed 100% program.\n";
close (OUT1);
close (OUT2);
close (OUT3);
close (OUT4);
### Return the start position of a near-exact match in the target sequence
### Yes, the search algorithm is inefficient.
  sub FUZZYMATCH ($$$) {
	my $primer = shift;
	my $searchSeq = shift;
	my $maxMM = shift;
	my $MM;
	my $startPos=0;
    my $pLen = length $primer;
	my $sLen = length $searchSeq;
	do {	
        $MM=$pLen - ( substr($searchSeq, $startPos,$pLen) ^ $primer ) =~ tr[\0][\0];
		$startPos++;
    } while($MM > $maxMM && $startPos <= $sLen-$pLen);		
	if ($startPos > $sLen-$pLen) { return -1; }
	else { return $startPos-1; }
  }

###Find the last index of continuous repeat array
   sub CONTIREPEAT($$){
   my $Sequence = shift;
   my $Repeat = shift;
   my $offset=0;
   my @Index = ();
   my @SortRepeat = ();
   my %lastIndex=();
   my $count=1;
   my $error=0;
   my $result=index($Sequence,$Repeat,$offset);
   while ($result != -1) {
		push (@Index, $result);
         $offset = $result + 1;
        $result = index($Sequence,$Repeat,$offset);
		}
   if(scalar @Index>=3){
   for (my $i=0; $i < scalar @Index-1; $i++){
		if($Index[$i+1]-$Index[$i]==length $Repeat){
		$count++;
		$lastIndex{$count}=$Index[$i+1];
		} elsif($Index[$i+1]-$Index[$i]==2*(length $Repeat) && $error<=2 ){
		$error++;
		$count=$count+2;
		$lastIndex{$count}=$Index[$i+1];
		} else{
		$lastIndex{$count}=$Index[$i];
		$count=1;
		}
		}
		@SortRepeat = sort {$b<=>$a} keys %lastIndex;
	if($SortRepeat[0]>=4){
     return ($lastIndex{$SortRepeat[0]},$SortRepeat[0]);
    } else{
	    return (-1,-1);
	 }
	 } else{
	 return (-1,-1);
	 }
	 }
	 
##fuzzy longest common substring based on hamming distance 
 sub lc_substr($$) {
  my $str1=shift;
  my $str2=shift;
  my %similarity=(); # difference percentage
  my @percentage=();
  my $nd;  # number of difference
  my $len1 = length $str1; 
  my $len2 = length $str2; 
 for (my $n1=0; $n1<$len1-3;$n1++){
     $nd =($len2-$n1) - (substr($str1,$n1) ^ substr($str2,0,$len2-$n1))=~ tr[\0][\0];
	 $similarity{$n1}= $nd/($len2-$n1);
   }
   @percentage = sort {$similarity{$a} <=> $similarity{$b}} keys %similarity;
  return $len1-$percentage[0];
 }
 

 ####call_multinomial_snp function to call if it is a SNP or sequencing error

 sub call_multinomial_snp($$$$){
 	my $snps = shift;
  	my $col = shift;
  	my $n = shift;
 	my $record_snps = shift;
 	my @nuc;
	my $snp = new SNP("",0,0,0,0,0,0);
	my $res;
	 
 	my $total = 0;
 	my $j = 0;
	my $l_ratio = 0;
 	foreach my $i ( sort keys %{$n}){
 		if($i ne 'N'){
 			$total = $total+$n->{$i};
 			$nuc[$j][0] = $i;
 			$nuc[$j][1] = $n->{$i};
 			$j = $j+1;
 		}
 	}
 	@nuc = sort { $b->[1] <=> $a->[1] } @nuc;
 	if($nuc[0][1] == 0){
 		if($record_snps){
 			$snp = new SNP ("snp_type_unk",$col,0,'N','-',0,0);
 			#$snp->{type} = "snp_type_unk";
 			#$snp->{col} = $col;
 			#$snp->{lratio} = 0;
 			#$snp->{rank_1} = 'N';
 			#$snp->{rank_2} = '-';
			
 			push (@{$snps}, $snp);  
 		}
 		return "snp_type_unk";
 	}
 	my $nuc_1 = $nuc[0][1];
 	my $nuc_2 = $nuc[1][1];
 	my $nuc_3 = $nuc[2][1];
 	my $nuc_4 = $nuc[3][1];
 	
	
 	$l_ratio = ($nuc_1 * log($nuc_1 / $total));

 	if($total - $nuc_1 > 0){
 		$l_ratio = $l_ratio + (($total - $nuc_1) * log(($total - $nuc_1) / (3 * $total)));
 	}
	
 	if ($nuc_1 + $nuc_2 > 0){
 		$l_ratio = $l_ratio - (($nuc_1 + $nuc_2) * log(($nuc_1 + $nuc_2) / (2 * $total)));
 	}
	    
 	if ($nuc_3 + $nuc_4 > 0){
 		$l_ratio = $l_ratio - (($nuc_3 + $nuc_4) * log(($nuc_3 + $nuc_4) / (2 * $total)));
 	}
	
 	$l_ratio = $l_ratio * 2;

 	 if($l_ratio <= $heterozygote_limit){
 		 ### This locus is a heterozygous
		 
 		 if($record_snps){
 			 $snp = new SNP("snp_type_het",$col,$l_ratio,$nuc[0][0],$nuc[1][0],0,0);
  			 #$snp->{type} = "snp_type_het";
  			 #$snp->{col} = $col;
  			 #$snp->{lratio} = $l_ratio;
  			 #$snp->{rank_1} = $nuc[0][0];
  			 #$snp->{rank_2} = $nuc[1][0];
			 
 			 push (@{$snps}, $snp);
 		 }
 		 $res = "snp_type_het";
 	 } elsif($l_ratio >= $homozygote_limit){
 		 if($record_snps){
 			 $snp = new SNP("snp_type_hom",$col,$l_ratio,$nuc[0][0],'-',0,0);
  			 #$snp->{type} = "snp_type_hom";
  			 #$snp->{col} = $col;
  			 #$snp->{lratio} = $l_ratio;
  			 #$snp->{rank_1} = $nuc[0][0];
  			 #$snp->{rank_2} = '-';
			 
 			 push (@{$snps}, $snp);
 		 }
 		 $res = "snp_type_hom";
 	 } else{
 		 if($record_snps){
 			 $snp = new SNP("snp_type_unk",$col,$l_ratio,$nuc[0][0],($nuc[1][1] > 0)? $nuc[1][0]: '-',0,0);
  			 #$snp->{type} = "snp_type_unk";
  			 #$snp->{col} = $col;
  			 #$snp->{lratio} = $l_ratio;
  			 #$snp->{rank_1} = $nuc[0][0];
  			 #$snp->{rank_2} = ($nuc[1][1] > 0)? $nuc[1][0]: '-';
			 
 			 push (@{$snps}, $snp);
 		 }
 		 $res = "snp_type_unk";
 	 }
	 
 	 return $res;
 }
 
