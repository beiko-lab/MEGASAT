use strict;
use warnings;
## hash table to save the correct genotype
my %lociNewScore = ();
## first argument which is the scores text file
my $inputScores = $ARGV[0];
## second argument which is the original genotype text file
my $originalMatrix = $ARGV[1];
## third argument which is the directory you want to save the output
my $SaveDir = $ARGV[2];
my $prefix = "New";
## open scores text file and put the values into a hash table. Concatenate loci names and original genotypes and store
##them as the keys of hash table, and the new genotypes are the values of corresponding keys.
open (my $scores,'<',$inputScores) or die "Cannot open socres txt file";
while (my $line = <$scores>) {
	chomp $line;
	my @fields = split "," , $line;
	if (scalar @fields >= 5) {
	   $fields[4]=~ s/\r$//;
	   @{$lociNewScore{$fields[0].$fields[1].$fields[2]}} = ($fields[3],$fields[4]);
	} else{
	   die "The format of scores txt file is wrong!";
	}
}
## get the name of original genotype text file
my $originalName = substr($originalMatrix,rindex($originalMatrix,'/')+1);
my $newName = $prefix.$originalName;
## generate a new text file that contains the updated genotypes
open(my $fh, ">", "$SaveDir/$newName") or die "Cannot open new txt file";

## open the original genotype text file and read first line to get the loci names, then store them into an array
open (my $data, '<', $originalMatrix) or die "Cannot open the origianl genotype txt file"; 
my $firstLine = <$data>;
my @lociName = split "\t" ,$firstLine;
$lociName[scalar @lociName -1]=~ s/\r$//;
print $fh "$firstLine"; 
## read every line after the first line.
while (my $line = <$data>) {
	chomp $line;
	my @Genotype = split "\t" , $line;
## Windows and Unix (including OS X) use different ways to express the end of a line. That regex deletes both kinds
## ensuring it will work no matter which type of machine produced the file or which type is reading it.
    $Genotype[scalar @Genotype -1]=~ s/\r$//;
    print $fh "$Genotype[0]\t";
	for (my $i =1; $i< scalar @Genotype-2;$i+=2){
	     my $LNKey = $lociName[$i];
	     my $Oldscore = $Genotype[$i].$Genotype[$i+1];
	     if(exists($lociNewScore{$LNKey.$Oldscore})){
	         print $fh "$lociNewScore{$LNKey.$Oldscore}[0]\t";
	         print $fh "$lociNewScore{$LNKey.$Oldscore}[1]\t";
	     } else{
	         print $fh "$Genotype[$i]\t";
	         print $fh "$Genotype[$i+1]\t";
	     }
	}
### After print the last genotype, we need to change to another line. That's why we separate last genotype from other
### genotypes for each line. 
	my $lastKey = $lociName[-2];
	my $lastOldscore = $Genotype[-2].$Genotype[-1];
	if(exists($lociNewScore{$lastKey.$lastOldscore})){
	     print $fh "$lociNewScore{$lastKey.$lastOldscore}[0]\t";
	     print $fh "$lociNewScore{$lastKey.$lastOldscore}[1]\n";
	  } else{
	         print $fh "$Genotype[-2]\t";
	         print $fh "$Genotype[-1]\n";
	     }
}
close($data);
close($fh);
print STDOUT "The program is completed!";

