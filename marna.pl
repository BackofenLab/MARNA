#!/usr/bin/perl

# MARNA: Multiple Alignment and Consensus Structure Prediction
#        of RNAs based on Sequence Structure Comparisons.
# Bioinformatics 2005, 21: 3352-3359
# (c) Sven Siebert

# used variables :
# @NAMES[$i]
# @SEQ[$i]
# @STR[$i][$j]
# @SIZE_STR[$i]

@NAMES=();
@SEQ=();
@STR=();
@SIZE_STR=();

use Getopt::Std;

# print out the help and exit
sub usage{
  print "call: marna.pl [-g,-h,-n,-s] FILENAME\n";
  print "FILENAME : file in fasta (-like) format. This format has to meet the following conditions :\n";
  print "           RNA sequences consist of names, sequences and structures which are represented\n";
  print "           in the file as :\n";
  print "              >seq1_name                   e.g.  >RNA1\n";
  print "              sequence                           GCGACGACUGUCGU\n";
  print "              structure_1                        ((((......))))\n";
  print "              ...                                .((((....)))).\n";
  print "              structure_i                        >RNA2\n";
  print "              >seq2_name                         ACGUGACUGUACGUCAGU\n";
  print "              sequence                           ((.((((.....))))))\n";
  print "              structure_1                        >RNA3\n";
  print "              ...                                ACGCGAUCUAGGCAGUGCU\n";
  print "              structure_j\n";
  print "           If no structure is given to a sequence, then they are automatically generated via\n";
  print "           the -g option.\n";
  print "options: -g structure : how to generate structure(s) where missing.\n";
  print "                         structure=0 (default) : mfe (RNAfold)\n";
  print "                         structure=1           : RNAshapes\n";
  print "                         structure=2           : RNAsubopt\n";
  print "         -h           : this help\n";
  print "         -n num       : number of suboptimal structures (needed for RNAsubopt).\n";
  print "         -p scorefile : reads alignment parameters from file (see readme.txt for description) .\n";
  exit;
}

my %options;
getopts('hn:g:p:', \%options);   # read the options with getopts
                                 # h : help
                                 # g : how to generate structures (see help)
                                 # n : number of suboptimal structures (needed for RNAsubopt)
                                 # p : score parameter file 

$options{h} && usage();     # the -h switch
my $NUM_STR=$options{n}||3 ;
my $STR_GEN=$options{g}||0 ;
my $SCORE_FILE=$options{p}||"" ;

print "MARNA: Multiple Alignment and Consensus Structure Prediction\n";
print "of RNAs based on Sequence Structure Comparisons.\n";
print "Bioinformatics 2005, 21: 3352-3359\n";
print "(c) Sven Siebert\n\n";

#check if sequence file exists
if ($ARGV[0] eq ""){
  print STDERR "An error has occurred. Filename missing.\n";
  exit 1;
}

# check if t_coffee is installed
# if you want to call t_coffee with absolute path then 1)modify line 186 and 2) comment the next 8 lines
system("whereis t_coffee | sed -e \"s/^[^/]*//\" | tr -d \"\n\" > tmp_str");
open (TMPSTR,"tmp_str");
if (-z TMPSTR){
  print STDERR "An error has occurred. t_coffee is not installed. Install it and then proceed.\n";
  exit 1;
}
close(TMPSTR);
unlink tmp_str;

print "Scores used:\n";
#check if parameterfile exists
if ($SCORE_FILE ne ""){
  if (!-e $SCORE_FILE){
    print STDERR "An error has occurred. Parameter file doesn't exist.\n";
    exit 1;
  }
  else{
    open (SCORES,$SCORE_FILE);
    print "base deletion (w_d) = ".<SCORES>;
    print "base mismatch (w_m) = ".<SCORES>;
    print "arc  removing (w_r) = ".<SCORES>;
    print "arc  breaking (w_b) = ".<SCORES>;
    print "arc  mismatch (w_am)= ".<SCORES>;
    close SCORES;
  }
}
else{
  print "base deletion (w_d) = 2.0\n";
  print "base mismatch (w_m) = 1.0\n";
  print "arc  removing (w_r) = 2.0\n";
  print "arc  breaking (w_b) = 1.5\n";
  print "arc  mismatch (w_am)= 1.8\n";
}

#check sequences
system ("./check_sequences.pl -g $STR_GEN -n $NUM_STR $ARGV[0] > $ARGV[0].formatted");
die "\n" if ($?!=0);

print "\nMissing structures are generated ";
if ($STR_GEN==0){
  print "with RNAfold to compute the minimum free energy (mfe) structures.\n";
}
elsif ($STR_GEN==1){
  print "with RNAshapes to compute ensembles of shaped structures.\n";
}
else{
  print "with RNAsubopt to compute random samples of ".$NUM_STR." structures.\n";
}
print "All sequences and structures can be viewed in the file $ARGV[0].formatted.\n";

# read data
open(DATA,"<$ARGV[0].formatted") or die "can't open file";
my $NUM_SEQ=-1;
while (defined($LINE=<DATA>)){  # read fasta file line by line
  chomp $LINE;

  if ( $LINE=~s/^\s*>// ){      # line begins with ">"
    $NUM_SEQ++;
    $i=0  ;
    $NAMES[$NUM_SEQ]=$LINE;
 }
  elsif ($LINE=~m/[ACGU]/ ){    # read sequence
    $SEQ[$NUM_SEQ]=$LINE ;
  }
  else  {                       # read structures
    $STR[$NUM_SEQ][$i++]=$LINE ;
    $SIZE_STR[$NUM_SEQ]=$i;
  }
}
close DATA;

# $NUM_SEQ is last index
open(LIB,">coffee.lib") or die "can't open file";
print LIB ($NUM_SEQ+1)."\n";

for($i=0;$i<=$NUM_SEQ;$i++){
  print LIB "seq".($i+1)."F ".(length $SEQ[$i])." ".$SEQ[$i]."\n";
}
close LIB;

print "Pairwise comparisons...\n";

# $NUM_SEQ is last index
# the one and only computation
$VIS=0;
for($i=0;$i<$NUM_SEQ;$i++){
  for($j=$i+1;$j<=$NUM_SEQ;$j++){
    $PER_CENT=200*($VIS++)/($NUM_SEQ*($NUM_SEQ+1));
    print STDERR "\r"; printf "%.0f", $PER_CENT; print "%";
    $i1=$i+1 ; $j1=$j+1;
    system("echo \"\#$i1 $j1\" >> coffee.lib");
    for($is=0;$is<$SIZE_STR[$i];$is++){
      for($js=0;$js<$SIZE_STR[$j];$js++){
	if ($SCORE_FILE ne ""){
	  system("./align -s $SCORE_FILE -t $SEQ[$i] \"$STR[$i][$is]\" $SEQ[$j] \"$STR[$j][$js]\" >> coffee.lib");
	}
	else{
	  system("./align -t $SEQ[$i] \"$STR[$i][$is]\" $SEQ[$j] \"$STR[$j][$js]\" >> coffee.lib");
	}
	die if ($?!=0);
      }
    }
  }
}
print "\rready.\n";

open(LIB,">>coffee.lib") or die "can't open file";
print LIB "CPU 0\n! SEQ_1_TO_N\n";
close LIB;

print "Build multiple alignment...";
# if t_coffee is properly installed then the next command should succeed,
# if not, then add the absolute path of t_coffee, e.g. : system("/usr/local/bin/t_coffee ....
system("t_coffee -in=Lcoffee.lib,Mclustalw_pair &> /dev/null");

print "\nready.\n\n";

# output sequence structure alignment

open(COFALN,"<coffee.aln")    or die "cannot open coffee.aln ."   ;
open(POSALN,">positions.aln") or die "cannot open positions.aln .";
@COFFEE_ALN=();
chomp(@COFFEE_ALN=<COFALN>);


print "\r";
$FIRST_RUN=0;
$LONGEST_NAME=0;
@ALIGNMENT_SEQ=();
@ALIGNMENT_STR=();

for($i=0;$i<=$NUM_SEQ;$i++){
  $k=$i+1;

  @ALIGN_LINES=();
  @ALIGN_LINES=grep /seq($k)F/, @COFFEE_ALN;
  chomp(@ALIGN_LINES);

  $ALIGN_SEQ=join "",@ALIGN_LINES;
  $ALIGN_SEQ=~s/seq($k)F//g;
  $ALIGN_SEQ=~s/\s//g;
  $ALIGN_SEQ=~tr/acgu/ACGU/;


  push(@ALIGNMENT_SEQ,$ALIGN_SEQ);

  # align structures according to sequences
  # first @SEQA = aligned sequence from t_coffee,
  #       @STRA = structure, if structure given, otherwise a sequence of stars '*'

  @SEQA=();
  @STRA=();

  @SEQA=split "",$ALIGN_SEQ;
  if ($SIZE_STR[$i]==1 ){          # display structure only if exact one structure available
    @STRA=split "",$STR[$i][0];
  }
  else{
    $STRSEQ="*"x(length $SEQ[$i]);
    @STRA=split "",$STRSEQ;
  }

  # second  insert gaps in @SEQA and @STRA and generate positions
  $C=0;
  @POS=();
  for($j=0;$j<=$#SEQA;$j++){
    if ($SEQA[$j] eq "-"){
      splice @STRA,$j,0,"-";
      @POS[$j]=-1;
    }
    else{
      $C++;
      @POS[$j]=$C;
    }
  }
  $POSSEQ=join " ",@POS;
  if($FIRST_RUN==0){                # write number of sequences and alignment length into position file
    print POSALN $NUM_SEQ." ".($#POS+1)."\n";
    $FIRST_RUN=1;
  }
  print POSALN "$POSSEQ\n";
  $ALIGN_STR=join "",@STRA;
  push(@ALIGNMENT_STR,$ALIGN_STR);

  if ((length $NAMES[$i])>$LONGEST_NAME){ $LONGEST_NAME=(length $NAMES[$i]);}
}

#real output
for($i=0;$i<=$NUM_SEQ;$i++){
  $NUM_SPACE=int ($LONGEST_NAME+4-(length $NAMES[$i]));
  $MAX_NAME_LENGTH=((length $NAMES[$i])+$NUM_SPACE);
  print $NAMES[$i]." "x$NUM_SPACE.$ALIGNMENT_SEQ[$i]."\n";
}
print "\n";
for($i=0;$i<=$NUM_SEQ;$i++){
  $NUM_SPACE=int ($LONGEST_NAME+4-(length $NAMES[$i]));
  print $NAMES[$i]." "x$NUM_SPACE.$ALIGNMENT_STR[$i]."\n";
}

close(POSALN);
close(COFALN);
print "\n";

print "Consensus:\n";
system("./consensus positions.aln $ARGV[0].formatted $MAX_NAME_LENGTH");
