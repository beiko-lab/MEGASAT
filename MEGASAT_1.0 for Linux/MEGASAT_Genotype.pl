### Separates files by primer
### TO DO:
###   Trim off primers (wobble = 2)
###   Filter weirdos (might need to process if coverage is low)
###   Look at length distribution, find modes

use strict;
use warnings;

# The function to match primer sequences to reads, with a tolerance for mismatches (see definition at the bottom of this file)
sub FUZZYMATCH ($$$);
# The function to find the starting position of repeat array
sub CONTIREPEAT($$);
# The function to find fuzzy longest common substring based on hamming distance 
sub lc_substr($$);

my %primerSeqs = ();
# Array that contains all the loci name
my @MS= ();
my $suffix="-b";
# A tab-separated file with column 1 = Locus name, column 2 = Forward primer, column 3 = Reverse primer, column 4 = 3' flank
# column 5 = 5' flank, column 6 = repeat array
# A header line is required
my $inputPrimers = $ARGV[0];
# The maximum # of mismatches (if parameter not given, max mismatches = 0, i.e. exact matches are required)
my $mismatches = $ARGV[1];
# The data set folder that contains all the input sequence read files (fastq file)
my $dataset = $ARGV[2];
# The directory to save the output folder
my $saveDir = $ARGV[3];
# get the name of data set
my $nameofdataset = substr($dataset,rindex($dataset,'/')+1);
my $path="$saveDir/Output_$nameofdataset";
# create the output folder
mkdir $path, 0755;
# disable perl output buffering
select((select(STDOUT), $|=1)[0]);
### Read the primers txt file and associate loci names with primers and flank information
open (IN,$inputPrimers) or die "Cannot open primers txt file!";
my $toss = <IN>;
while (<IN>) {
	chomp;
	my @line = split /\t/,$_;
	if (scalar @line >= 6) {
		@{$primerSeqs{$line[0]}} = @line[1..4];
		$line[5] =~ s/\r$//;
		push @{$primerSeqs{$line[0]}},$line[5];
	} else{
	    die "The format of primers txt file is wrong!";
	}
}
close (IN);

### put all the loci name into the array @MS
foreach my $MSname (sort keys %primerSeqs) {
    push (@MS, $MSname);
	push(@MS,$MSname.$suffix);
	}
### produce a txt file that has all the genotype information for all the individuals and loci
open (OUT1, ">$path/Genotype.txt") or die "Cannot generate Genotype.txt file";
print OUT1 "Sample_idx1_idx2\t", join "\t",@MS,"\n";

# open the data set folder and start reading fastq file
opendir(DIR, "$dataset" ) or die "Cannot open data set folder!";
my @files = grep(/\.fastq$/,readdir(DIR));
closedir DIR;
my $index_file = 0;
if(scalar @files ==0){
    die "No fastq file in the data set";
}
foreach my $file(@files){
         $index_file++;
         my %seqsMapped = ();     ## hash table that contains sequences that are trimmed off forward primers
         my %seqsTrimmed = ();    ## hash table that contains sequences that are trimmed off forward primers and reverse primers
         open (IN, "$dataset/$file") or die "Cannot open $file!";


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
		my $foundForward = undef;
		my $foundReverse = undef;		
		foreach my $checkPrimer (keys %primerSeqs) {
			if (!defined( )) {
				die "$checkPrimer";
			 }
			my $isMatched = FUZZYMATCH($primerSeqs{$checkPrimer}[0],$thisLine,$mismatches);	
			# Is there a fuzzy forward primer match?
			if ($isMatched >= 0) {
				$foundForward = $checkPrimer;
				substr($thisLine,0,$isMatched + length $primerSeqs{$checkPrimer}[0]) = "";
				push (@{$seqsMapped{$foundForward}},$thisLine);	
			 }
		}
		if (defined($foundForward)) {
		    my ($LRIndex,$NumRepeat) = CONTIREPEAT($thisLine,$primerSeqs{$foundForward}[4]);
			my $revPrimer = $primerSeqs{$foundForward}[1];
			my $threeFlank = $primerSeqs{$foundForward}[2];
			my $fiveFlank = $primerSeqs{$foundForward}[3];
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
						         my $addString="a" x (200-(length $thisLine)+$startRepeat);
							     $thisLine .=$addString;
						         push @{$seqsTrimmed{$foundForward}},$thisLine;
							} elsif((((length $thisLine)-$endRepeat)- (substr($thisLine,$endRepeat) ^ substr($threeFlank,0,(length $thisLine)-$endRepeat))=~ tr[\0][\0])<=2){
				                my $leftTF=substr($threeFlank,(length $thisLine)-$endRepeat);
				                $thisLine .=$leftTF;
				                push @{$seqsTrimmed{$foundForward}},$thisLine;
				               } 
				        } 
				}
			}
		}
	}
close (IN);

# Create the prefixes for the output files
my $prefix = $file;
$prefix =~ s/_.*//;
# Print the full sequences that are trimmed off forward primers
foreach my $outFile (keys %seqsMapped) {
    open (OUT, ">$path/Sorted_${prefix}_$outFile.split");
	foreach my $outSeq (@{$seqsMapped{$outFile}}) {
		print OUT "$outSeq\n";
		}
	close (OUT);
}

# Print the trimmed sequences associated with each primer pair, *and* count up the occurrences of different lengths
# hash table that count the occurrences of different lengths for different locus
my %lengths = ();
# has table that count the occurrences of different lengths
my %lenRange = (); 
# array that preserve the length range
my @MSlens = ();

# Print the trimmed sequences associated with each primer pair
foreach my $outPrimer (keys %seqsTrimmed) {
	open (OUT, ">$path/Trimmed_${prefix}_$outPrimer.split");
	foreach my $outSeq (@{$seqsTrimmed{$outPrimer}}) {
		print OUT "$outSeq\n";
		my $sLen = length $outSeq;
		$lenRange{$sLen} = 1;
		++$lengths{$outPrimer}{$sLen};
	}
	close (OUT);
}


# Print out a table showing the length distribution of each microsatellite locus
print OUT1 "$prefix";

@MSlens = sort {$a <=> $b} keys %lenRange;
my $minMS = $MSlens[0];
my $maxMS = $MSlens[-1];
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
	     if($sum>=20){
	        if($sorted[0] < $sorted[1]){
	            my $ratio= $lengths{$outMS}{$sorted[1]}/$lengths{$outMS}{$sorted[0]};
	            if($ratio>=0.15){
	               if((defined($sorted[2]) && $sorted[2] > $sorted[1]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[1]}>=0.4) || (defined($sorted[3]) && $sorted[3] > $sorted[1]&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[1]}>=0.4)){
		              if($sorted[2]-$sorted[1]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[1]}>=0.7){
			              push (@Plength, ($sorted[0],$sorted[2]));
				        } elsif(defined($sorted[3]) && $sorted[3]-$sorted[1]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[1]}>=0.7){
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
	       if($sorted[0] > $sorted[1]){
	         if(($sorted[0]-$sorted[1]>=3&&$lengths{$outMS}{$sorted[1]}>=0.6*$lengths{$outMS}{$sorted[0]})||($sorted[0]-$sorted[1]<=2&&$lengths{$outMS}{$sorted[1]}>=0.8*$lengths{$outMS}{$sorted[0]})){
	              if(defined($sorted[2]) && $sorted[2] > $sorted[0]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[0]}>=0.2){
		             push (@Plength, ('Unscorable','Unscorable'));
		            } else{
		              push (@Plength, ($sorted[1],$sorted[0]));
		            }
	           } elsif(defined($sorted[2]) && $sorted[2] > $sorted[0]&& $lengths{$outMS}{$sorted[2]}/$lengths{$outMS}{$sorted[0]}>=0.2){
	                 if(defined($sorted[3]) && $sorted[3]-$sorted[2]== (length $primerSeqs{$outMS}[4])&& $lengths{$outMS}{$sorted[3]}/$lengths{$outMS}{$sorted[2]}>=0.7){
	                    push (@Plength, ($sorted[0],$sorted[3]));
	                   } else{
	                     push (@Plength, ($sorted[0],$sorted[2]));
		                }
	          } else{
	               push (@Plength, ($sorted[0],$sorted[0]));	                
	          }
	        }
      } else{
        push (@Plength, (0,0));
       }
   } elsif(scalar keys %{$lengths{$outMS}}==1){
       @sorted = sort {$lengths{$outMS}{$b} <=> $lengths{$outMS}{$a}} keys %{$lengths{$outMS}};
       if($lengths{$outMS}{$sorted[0]}>=20){
          push (@Plength, ($sorted[0],$sorted[0]));
         } else{
           push (@Plength, (0,0));
         }
   } else{
     push (@Plength, ('X','X'));
   }
   print OUT1 "\t$Plength[0]\t$Plength[1]";
   print OUT ",@Plength\n"; 
}
  print OUT1 "\n";
  close (OUT);
  my $n_task = sprintf("%.3f",$index_file/(scalar @files));
  my $per_task = $n_task*100;
  print STDOUT "Completed $per_task% program.\n"; 
}
close (OUT1);
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