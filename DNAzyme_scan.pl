#!/usr/bin/perl -w
use strict;

use Getopt::Long qw(:config no_ignore_case bundling);
use List::Util qw(max min);

#=begin Program Description
#
#	DNAzyme_scan version 0.3
#   Scanning RNA sequences for putative DNAzyme target sites
#	Antonio Marco, 2005-2022
#	University of Essex
#	Basic usage:
#		perl DNAzyme_scan.pl -i input_file -o output_file
#
#	This program is free to use and distribute for non-commercial purposes
#	For commercial use please contact the author
#=cut

# INPUT/OUTPUT parameters (default values)
my $help = 0;
my $infile = "infile";
my $outfile = "outfile";
my $type = "10-23";
my $only_cds = 0;
my $arm_length = 9;
my $linker = 'ggctagctacaacga';# for 10-23
my $dnazyme_type = "10-23";


# Quick help if no parameters are provided
unless((@ARGV) || (-e $infile)){ # if no arguments AND infile does not exist
	print STDERR "\nNo arguments specified\nType \"DNAzyme_scan.pl -h\" for help\n\n";
	exit;
}

# Modify parameters according to the ARGUMENT line (user defined)
GetOptions(
	"infile=s"=>\$infile, "i=s"=>\$infile,
	"outfile=s"=>\$outfile, "o=s"=>\$outfile,
	"armlength=i"=>\$arm_length, "l=i"=>\$arm_length,
	"cds"=>\$only_cds,  "c"=>\$only_cds,
	"help"=> \$help, "h"=> \$help
);

# Print HELP screen
if($help){
        system("clear");
	print helpScreen();
	exit;
}


# PRINT RUNNING HEADER
system("clear");
print STDERR runningHeader();


open(QUERY,$infile);
my %query_transcripts = hashFromFasta(join('',<QUERY>));
close QUERY;

my $linker = 'ggctagctacaacga';# for 10-23
my $dnazyme_type = "10-23";


open(OUTFILE,">$outfile");
print OUTFILE "Transcript\tType\tPosition\tRNA\tDNAzyme\tDeltaG\n";

my $total = scalar keys %query_transcripts;
my $count = 0;
my $number_or_putative_dnazymes = 0;
foreach my $query (keys %query_transcripts){
    $count++; # How many transcripts has been explored so far
	my $total_length = length($query_transcripts{$query});

    # Detect Open Reading Frames
    my @CDS_start_end = (0,$total_length); 
    if($only_cds == 1){@CDS_start_end = start_end_cds($query_transcripts{$query});}
  
    my $last_char = " - "; # For progress bar
    while($query_transcripts{$query} =~ /[AG][UTC]/gi ){ # 10-23     # FOR EACH POTENTIAL SITE
        my $score = 100; # TO BE MODIFIED
        my $upstream_cleavage = pos($query_transcripts{$query})-1;
        my $start_left_arm = $upstream_cleavage - $arm_length;
        
        # Progress bar
        my $percent_process = ($upstream_cleavage*100)/$total_length;
        if($only_cds==1){$percent_process = (($upstream_cleavage-$CDS_start_end[0])*100)/($CDS_start_end[1]);if($percent_process>100){$percent_process=100;}}
        $last_char = rotating_char($last_char);
        print STDERR "Scanning transcript $query ($count out of $total)".rotating_char($last_char).sprintf("%.0f", $percent_process)." \%    \r";
        
        unless(($start_left_arm<0)||(($start_left_arm+($arm_length*2)+2)>$total_length)||($start_left_arm<$CDS_start_end[0])||(($start_left_arm+($arm_length*2)+2)>$CDS_start_end[1])){
        # avoid partial matches at both ends, and ignore outside CDS if parameter set
            my @why_fileterd_out = ();
            my $left_RNA = substr($query_transcripts{$query},$start_left_arm - 1,$arm_length + 1);
            my $right_RNA = substr($query_transcripts{$query},$start_left_arm+$arm_length,$arm_length);
            my $left_DNA = reverseComp($right_RNA); # ADD CONTROL FOR NO NUCLEOTIDES OR TYPOS
            my $right_DNA = reverseComp($left_RNA);
            $right_DNA =~ s/^.//;
            
            $left_RNA =~ tr/T/U/;
            $right_RNA =~ tr/T/U/;
            
            my $DNAzyme = $left_DNA .$linker.$right_DNA;
            my $DNAzyme_branches = $left_DNA.$right_DNA;
            my $target_RNAsite = $left_RNA.$right_RNA;

            # hybridization energy and difference with RNA inner energy
            my $energy_DNA_RNA = free_energy_DNA_RNA_duplex($DNAzyme_branches); # Assumes is perfectly aligned        

            # OUTPUT INFORMATION
            my $output = $query ."\t";
            $output .= $dnazyme_type."_target\t". $start_left_arm ."-". ($start_left_arm + ($arm_length*2)) ."\t";
            $output .= $left_RNA ."*". $right_RNA. "\t". $DNAzyme;
            $output .= "\t". $energy_DNA_RNA;
			$number_or_putative_dnazymes++;
			print OUTFILE $output ."\n";
        } # END avoid partial matches at both ends, discard outside ORF if set
    }
}

print STDERR "\n";
print STDERR $number_or_putative_dnazymes . " putative DNAzymes identified (in $outfile file) \n\n";
close OUTFILE;

exit 0;



sub runningHeader{
	my $output = "============================================================\n";
	$output .= "=               DNAzyme_scan version 0.3                   =\n";
	$output .= "============================================================\n";
	$output .= "- Scanning RNA sequences for putative DNAzyme target sites -\n";
	$output .= "------------------------------------------------------------\n";
	$output .= "-                       Antonio Marco                      -\n";
	$output .= "-               University of Essex, 2015-2022             -\n";
	$output .= "------------------------------------------------------------\n\n";
	return $output;
}

sub helpScreen{
	my $output = runningHeader()."\n";
	$output .= "Usage:\n   perl DNAzyme_scan.pl -i <INPUT_FILE> -o <OUTPUT_FILE> [-OPTIONS]\n";
	$output .= "\nArguments [default values in brackets]:\n";
	$output .= "\t-h / --help\t\tPrint this screen\n";
	$output .= "\t-i / --infile\t\tInput fasta file of RNA transcript sequences [infile]\n";
	$output .= "\t-o / --outfile\t\tOutput table [outfile]\n";
	$output .= "\t-c / --cds\t\tScan only putative coding sequence [disable]\n";
	$output .= "\n\n\n";
	return $output;
}

sub hashFromFasta{
	my %OUT = ( );
	my @ELEMENTS = split(">",$_[0]);
	shift(@ELEMENTS);
	for (@ELEMENTS){
		my @LINES = split("\n",$_);
		my $name = $LINES[0];
		shift(@LINES);
		my $outseq = join("",@LINES);
		$outseq =~ s/[^A-Za-z0-9-~._,;:<>!?#$%&()=]//g; # there are also ? in our sequences
		$outseq =~ tr/[a-z]/[A-Z]/;
		$OUT{$name} = $outseq;
	}
	return %OUT;
}

sub reverseComp{
	my $dna = $_[0];
	my $revcomp = reverse($dna);
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}


sub free_energy_DNA_RNA_duplex{ # Kcal/mol at 37 degrees, 1M Na
    my $energy = 0;
    for(my $i=0;$i<(length($_[0])-1);$i++){
        $energy+=nearest_neighbor_energy_DNA_RNA_duplex(substr($_[0],$i,2));
    }
    return $energy + 3.1; # initiation energy added
}

sub nearest_neighbor_energy_DNA_RNA_duplex{# Kcal/mol at 37 degrees, 1M Na
    if($_[0] eq 'AA'){
        return -1.00;
    }elsif($_[0] eq 'TT'){
        return -1.00;
    }elsif($_[0] eq 'AT'){
        return -0.88;
    }elsif($_[0] eq 'TA'){
        return -0.58;
    }elsif($_[0] eq 'CA'){
        return -1.45;
    }elsif($_[0] eq 'TG'){
        return -1.45;
    }elsif($_[0] eq 'GT'){
        return -1.44;
    }elsif($_[0] eq 'AC'){
        return -1.44;
    }elsif($_[0] eq 'CT'){
        return -1.28;
    }elsif($_[0] eq 'AG'){
        return -1.28;
    }elsif($_[0] eq 'GA'){
        return -1.30;
    }elsif($_[0] eq 'TC'){
        return -1.30;
    }elsif($_[0] eq 'CG'){
        return -2.17;
    }elsif($_[0] eq 'GC'){
        return -2.24;
    }elsif($_[0] eq 'GG'){
        return -1.84;
    }elsif($_[0] eq 'CC'){
        return -1.84;
    }else{
        return 0;
    }
}

sub complementary{
    my $temp = $_[0];
    $temp =~ tr/ACGT/TGCA/;
    if($temp eq $_[1]){
        return 1;
    }else{
        return 0;
    }
}

sub start_end_cds{ # Returns the start and end position of an Open Reading Frame within a transcript
    my $transcript = $_[0];
    for(my $start=0;$start<length($transcript);$start++){
        if((substr($transcript,$start,3) eq 'AUG')||(substr($transcript,$start,3) eq 'ATG')){
            for(my $end=$start+3;$end<length($transcript);$end+=3){
                my $codon = substr($transcript,$end,3);
                if(($codon eq 'UAG')||($codon eq 'UAA')||($codon eq 'UGA')||($codon eq 'TAG')||($codon eq 'TAA')||($codon eq 'TGA')){
                    return ($start,$end);
                }
            }
        }
    }
    return (0,0);
}

sub rotating_char{
    if($_[0] eq " - "){
        return " \\ ";
    }elsif($_[0] eq " \\ "){
        return " | ";
    }elsif($_[0] eq " | "){
        return " / ";
    }elsif($_[0] eq " / "){
        return " - ";
    }else{
        return "   "
    }
}
