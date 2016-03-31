#!/usr/bin/perl

use Getopt::Std;
use File::Find;

my $NUM_SEQ    = -1   ;  # number of RNA sequences
my $MAX_LENGTH = 1000 ;  # the length of an RNA sequence should not exceed $MAX_LENGTH

# $STR_GEN     :      how to generate structure where missing (see usage,  default is mfe)
# $NUM_STR     :      number of structures assigned to each sequence where structures are missing
# @NAMES       :      contains the sequence names from 0,...,$NUM_SEQ-1
# @SEQ         :      contains the sequences from 0,...,$NUM_SEQ-1
# @STR[$i][$j] :      contains the j-th structure for sequence $i
# @SIZE_STR[$i]:      number of structures assigned to each sequence

# print out the help and exit
sub usage
  {
    print "check_sequences.pl [-n num,-g structure] FILENAME\n";
    print "    -n num       : number of suboptimal structures (needed for RNAsubopt).\n";
    print "    -g structure : how to generate structure(s) where missing.\n";
    print "                     structure=0 (default) : mfe (RNAfold)\n";
    print "                     structure=1           : RNAshapes\n";
    print "                     structure=2           : RNAsubopt\n";
    exit;
  }

if ($ARGV[0] eq ""){
  print STDERR "filename missing.\n";
  exit 1;
}

my %options;
getopts('hn:g:', \%options);  # read the options with getopts
                                 # h : help
                                 # g : how to generate structures (see help)
                                 # n : number of suboptimal structures (needed for RNAsubopt)

$options{h} && usage();     # the -h switch
my $NUM_STR=$options{n}||3 ;
my $STR_GEN=$options{g}||0 ;

sub check_name {
  if (length $_[0]==0){                                 # name length = 0
    print STDERR "An error has occurred. No sequence name given.\n";
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n";
    exit 1 ;
  }

  for($i=0;$i<$NUM_SEQ;$i++){                           # check for double names
    if($_[0] eq $NAMES[$i]){
      print STDERR "An error has occurred. Sequence name not unique:\n>".$_[0]."\n";
      exit 1;
    }
  }
}

sub check_seq {
  if (length $_[1]==0){                                 # sequence length = 0
    print STDERR "An error has occurred. Can't recognize sequence.\n";
    print STDERR "Detected a sequence with length 0.";
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n";
    exit 1 ;
  }
  if (length $_[1]>$MAX_LENGTH){                        # sequence length > MAX_LENGTH
    print STDERR "An error has occurred. Sequence length of ".(length $_[1]);
    print STDERR " nt exceeds the maximal size of ".$MAX_LENGTH." nt.\n";
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n";
    exit 1 ;
  }
  if (grep /[^ACGU]/,$_[1]){                            # illegal characters
    print STDERR "An error has occurred. Sequence contains illegal characters.\n";
    $ERR_SEQ=$_[1];
    $ERR_SEQ=~tr/ACGU/\^/c;
    $ERR_SEQ=~tr/ACGU/ /;
    print STDERR ">".$_[0]."\n".$_[1]."\n".$ERR_SEQ."\n";
    exit 1;
  }
}

sub check_str {                                         # check indiviudal structures
  if (grep /[^\(\)\.]/,$_[2]){                          # illegal characters
    print STDERR "An error has occurred. Structure contains illegal characters.\n";
    $ERR_STR=$_[2];
    $ERR_STR=~tr/\(\)\./\^/c;
    $ERR_STR=~tr/\(\)\./ /;
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n".$ERR_STR."\n";
    exit 1;
  }

  $LEFTBRACKETS=$RIGHTBRACKETS=$_[2];                    # check for balanced brackets
  $CL=$LEFTBRACKETS=~s/\(//g;               # count left brackets
  $CR=$RIGHTBRACKETS=~s/\)//g;              # count right brackets

  if($CL>$CR){
    print STDERR "An error has occurred. The structure occupies more left binding bases ";
    print STDERR "than right binding bases. Check your structure !\n";
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n";
    exit 1;
  }
  elsif($CR>$CL){
    print STDERR "An error has occurred. The structure occupies more right binding bases ";
    print STDERR "than left binding bases. Check your structure !\n";
    print STDERR ">".$_[0]."\n".$_[1]."\n".$_[2]."\n";
    exit 1;
  }
}

sub check_prog_exists{    # gets value 0,1 or 2 depending on choice for structure generating program
  if ($_[0] eq "0"){
    system("whereis RNAfold | sed -e \"s/^[^/]*//\" | tr -d \"\n\" > tmp_str");
  }
  elsif ($_[0] eq "1"){
    system("whereis RNAshapes | sed -e \"s/^[^/]*//\" | tr -d \"\n\" > tmp_str");
  }
  else {
    system("whereis RNAsubopt | sed -e \"s/^[^/]*//\" | tr -d \"\n\" > tmp_str");
  }
  open (TMPSTR,"tmp_str");
  if (-z TMPSTR){
    if    ($_[0] eq "0"){ print STDERR "RNAfold "  ; }
    elsif ($_[0] eq "1"){ print STDERR "RNAshapes "; }
    else                { print STDERR "RNAsubopt "; }
    print STDERR "does not exist. Install it before using it.\n";
    exit 1;
  }
  $EXECUTE=<TMPSTR>;
#  print "Execute : ".$EXECUTE."\n";
  close TMPSTR;
  unlink tmp_str;
}

# read data
open(DATA,"<$ARGV[0]") or die "can't open file";

$READ_NAME="";
$READ_SEQ ="";
$READ_STRS="";

while (defined($LINE=<DATA>)){                 # read fasta file line by line
  $LINE=~s/\r//g;                              # newlines in windoof : "\r\n"
  chomp $LINE;

  if($LINE=~s/^\s*>//){                        # line begins with ">"
    if ($READ_NAME ne "" && $READ_SEQ eq ""){  # sequence name given, but not sequence
      print STDERR "An error has occurred. Sequence missing.\n";
      print STDERR ">".$READ_NAME."\n";
      exit 1;
    }
    elsif ($READ_NAME ne "" && $READ_SEQ ne ""){  # not the first run, at least the 2nd,
      $NUM_SEQ++;                              # therefore be careful with NUM_SEQ,
                                               # it is the sequence number from the last run.
      # error detection of name and sequence
      check_name($READ_NAME,$READ_SEQ,$READ_STRS);
      check_seq ($READ_NAME,$READ_SEQ,$READ_STRS);

      push(@NAMES,$READ_NAME);
      push(@SEQ,$READ_SEQ);

      # check if and how many structures are given
      if ((length $READ_STRS)==0){             # no structure given
	$LENGTH=length $READ_SEQ;
	if ($STR_GEN==0){                      # RNAfold
	  check_prog_exists(0);
	  system("echo $READ_SEQ | RNAfold | cut -c 1-$LENGTH > tmp_str");
	  unlink "rna.ps";

	  open (TMPSTR,"tmp_str");
	  while(<TMPSTR>){
	    chomp($STR[$NUM_SEQ][0]=$_);
	  }
	  close TMPSTR;
	  unlink tmp_str;
	  $SIZE_STR[$NUM_SEQ]=1;
	}

	elsif ($STR_GEN==1){                   # RNAshapes
	  $LENGTH=$LENGTH+5;
	  check_prog_exists(1);
	  system("RNAshapes -s $READ_SEQ | grep \"\.\.\.\" | sed -n -e \"s/.*\\t\\([\\(\\)\\.]*\\)\\t.*/\\1/p\" > tmp_str");
	  open (TMPSTR,"tmp_str");
	  $i=0;
	  while(<TMPSTR>){
            if ($i<$NUM_STR){ chomp($STR[$NUM_SEQ][$i++]=$_);}
	  }
	  close TMPSTR;
#	  unlink tmp_str;
	  $SIZE_STR[$NUM_SEQ]=$i;
	}

	elsif ($STR_GEN==2){                   # RNAsubopt
	  check_prog_exists(2);
	  system("echo $READ_SEQ | RNAsubopt -p $NUM_STR | cut -c 1-$LENGTH > tmp_str");
	  open (TMPSTR,"tmp_str");
	  $i=0;
	  while(<TMPSTR>){
	    chomp($STR[$NUM_SEQ][$i++]=$_);
	  }
	  close TMPSTR;
	  if ($i!=$NUM_STR){
	    print STDERR "An error has occurred. Number of produced structures(".$i.") not ";
	    print STDERR "equal to what it should be(".$NUM_STR.")\n";
	    print STDERR ">".$READ_NAME."\n".READ_SEQ."\n".READ_STR."\n";
	    exit 1;
	  }
	  unlink tmp_str;
	  $SIZE_STR[$NUM_SEQ]=$i;
	}
	else{
	  print STDERR "An error has occurred. Don't know how to generate structures.\n";
	  exit 1;
	}
      }
      else{
	$factor=(length $READ_STRS)/(length $READ_SEQ);
	if($factor == int($factor)){           # $factor structures assigned to sequence
	  for($i=0;$i<$factor;$i++){
	    $STR[$NUM_SEQ][$i]=substr($READ_STRS,$i*(length $READ_SEQ),(length $READ_SEQ));

	    # error detection of structures
	    check_str ($READ_NAME,$READ_SEQ,$STR[$NUM_SEQ][$i]);
	  }
	  $SIZE_STR[$NUM_SEQ]=$factor;
	}
	else{
	  check_str ($READ_NAME,$READ_SEQ,$STR[$NUM_SEQ][0]);
	  print STDERR "An error has occurred. Structure length(s) not valid.\n";
	  print STDERR "Name      : >".$READ_NAME."\nSequence  : ".$READ_SEQ."\nStructure : ".$READ_STRS."\n";
	  exit 1;
	}
      }
    }
    $READ_NAME = $LINE; $READ_NAME=~s/ /_/g;
    $READ_SEQ  = ""   ;
    $READ_STRS = ""   ;
  }
  elsif ($LINE=~m/^\s*\w/ ){                    # read sequence
    $LINE=~s/\s//g;
    $LINE=~tr/[a-z]/[A-Z]/;
    $LINE=~tr/tT/U/;
    $READ_SEQ.=$LINE;
  }
  else {                                       # read structures
    $LINE=~s/\s//g;
    $READ_STRS.=$LINE;
  }
}

# last run
if ($READ_NAME ne "" && $READ_SEQ eq ""){  # sequence name given, but not sequence
  print STDERR "An error has occurred. Sequence missing.\n";
  print STDERR ">".$READ_NAME."\n";
  exit 1;
}

elsif ($READ_NAME ne "" && $READ_SEQ ne "" ){
  $NUM_SEQ++;

  # error detection of name and sequence
  check_name($READ_NAME,$READ_SEQ,$READ_STRS);
  check_seq ($READ_NAME,$READ_SEQ,$READ_STRS);

  push(@NAMES,$READ_NAME);
  push(@SEQ,$READ_SEQ);

  # check if and how many structures are given

  if ((length $READ_STRS)==0){             # no structure given
    $LENGTH=length $READ_SEQ;
    if ($STR_GEN==0){                      # RNAfold
      check_prog_exists(0);
      system("echo $READ_SEQ | RNAfold | cut -c 1-$LENGTH > tmp_str");
      unlink "rna.ps";
      open (TMPSTR,"tmp_str");
      while(<TMPSTR>){
	chomp($STR[$NUM_SEQ][0]=$_);
      }
      close TMPSTR;
      unlink tmp_str;
      $SIZE_STR[$NUM_SEQ]=1;
    }
    elsif ($STR_GEN==1){                  # RNAshapes
      $LENGTH=$LENGTH+5;
      check_prog_exists(1);
      system("RNAshapes -s $READ_SEQ | grep \"\.\.\.\" | sed -n -e \"s/.*\\t\\([\\(\\)\\.]*\\)\\t.*/\\1/p\" > tmp_str");
      open (TMPSTR,"tmp_str");
      $i=0;
      while(<TMPSTR>){
	if ($i<$NUM_STR){ chomp($STR[$NUM_SEQ][$i++]=$_);}
      }
      close TMPSTR;
      unlink tmp_str;
      $SIZE_STR[$NUM_SEQ]=$i;
    }
    elsif ($STR_GEN==2){                  # RNAsubopt
      check_prog_exists(2);
      system("echo $READ_SEQ | RNAsubopt -p $NUM_STR | cut -c 1-$LENGTH > tmp_str");
      open (TMPSTR,"tmp_str");
      $i=0;
      while(<TMPSTR>){
	chomp($STR[$NUM_SEQ][$i++]=$_);
      }
      close TMPSTR;
      if ($i!=$NUM_STR){
	print STDERR "An error has occurred. Number of produced structures(".$i.") not ";
	print STDERR "equal to what it should be(".$NUM_STR.")\n";
	print STDERR ">".$READ_NAME."\n".READ_SEQ."\n".READ_STR."\n";
	exit 1;
      }
      unlink tmp_str;
      $SIZE_STR[$NUM_SEQ]=$i;
    }
    else{
      print STDERR "An error has occurred. Don't know how to generate structures.\n";
      exit 1;
    }
  }
  else{
    $factor=(length $READ_STRS)/(length $READ_SEQ);
    if($factor == int($factor)){           # $factor structures assigned to sequence
      for($i=0;$i<$factor;$i++){
	$STR[$NUM_SEQ][$i]=substr($READ_STRS,$i*(length $READ_SEQ),(length $READ_SEQ));
	
	# error detection of structures
	check_str ($READ_NAME,$READ_SEQ,$STR[$NUM_SEQ][$i]);
      }
      $SIZE_STR[$NUM_SEQ]=$factor;
    }
    else{
      print STDERR "An error has occurred. Structure length(s) not valid.\n";
      print STDERR "Name      : >".$READ_NAME."\nSequence  : ".$READ_SEQ."\nStructure : ".$READ_STRS."\n";
      exit 1;
    }
  }
}
close(DATA);

if ($NUM_SEQ<0){
  print STDERR "An error has occurred. No sequence given.\n";
  exit 1;
}
elsif ($NUM_SEQ==0){
  print STDERR "An error has occurred. At least two sequence are needed.\n";
  exit 1;
}

for($i=0;$i<=$NUM_SEQ;$i++){
  print ">".$NAMES[$i]."\n";
  print $SEQ[$i]."\n";
  for($j=0;$j<$SIZE_STR[$i];$j++){
    print $STR[$i][$j]."\n";
  }
}
