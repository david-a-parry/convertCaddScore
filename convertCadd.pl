#!/use/bin/env perl
use strict;
use warnings;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Bio::DB::Sam;
use File::Temp qw/ tempfile /;
use Getopt::Long;

my $chain = "$ENV{HOME}/ref/hg19ToHg38.over.chain";
my $fasta = "$ENV{HOME}/ref/hg38/hg38.fa";
my $liftover = "liftOver";
GetOptions
(
    'f|fasta=s' => \$fasta,
    'c|chain=s' => \$chain,
    'l|liftover=s' => \$liftover,
    'h|help',
) or usage("Syntax error\n");
usage ( "Usage: $0 cadd_scores.tsv\n" ) if not @ARGV;
die "--chain file '$chain' does not exist!\n" if not -e $chain;
die "--fasta file '$fasta' does not exist!\n" if not -e $fasta;

my $fai = Bio::DB::Sam::Fai->load($fasta);
while (my $in = shift){
    parseFile($in);
}


##################################################
sub usage{
    my $msg = shift;
    print "$msg\n" if $msg;
    print <<EOT

convertCaddScores.pl - converts cadd score coordinates from hg19/GRCh37 to hg38

Usage: $0 cadd_scores.tsv [cadd_scores_2.tsv ...] [options]

Options:
    
    -f FILE --fasta FILE
        hg38 fasta file for checking REF alleles from converted CADD scores
        Default = $ENV{HOME}/ref/hg38/hg38.fa

    -c FILE --chain FILE
        hg19ToHg38 liftover chain file (can be found at UCSC)
        Default = $ENV{HOME}/ref/hg19ToHg38.over.chain

    -l LIFTOVER --liftover LIFTOVER
        Location of liftOver executable if not in your \$PATH

    -h --help
        Show this message and exit

Warning: 
    
    This program simply converts the genomic coordinates of CADD scores from 
    hg19/GRCh37 to hg38. Some will not be able to be converted due to 
    differences between the assemblies. No changes are made to scores that do 
    convert to reflect changes between the assemblies. 

    In order to run correctly you must have a working copy of the liftOver 
    commandline program from UCSC either in your \$PATH or specified using the 
    -l/--liftover option.

Author:

    David A. Parry
    University of Edinburgh
EOT
    ;
    exit 1 if $msg;
    exit;
}

##################################################
sub parseFile{
    my $f = shift;
    my $FH = openFile($f);
    my @contiguous = ();
    my $prev_chr = '';
    my $prev_pos = -999;
LINE: while (my $line = <$FH>){
        if ($line =~ /^#/){
            if ($line =~ /^##/){
                print $line;
            }else{
                print "##convertCadd.pl liftOver\n$line";
            }
            next LINE;
        }
        chomp $line;
        next LINE if not $line;
        $line =~ s/^chr//; 
        my @split = split("\t", $line);
        if ( @contiguous < 100000 and #limit the size of a region to process 
            ( $split[0] eq $prev_chr and 
            $split[1] <= $prev_pos + 1 and 
            $split[1] >= $prev_pos )
        #same as prev coordinate or is next coordinate
        ){
            push @contiguous, \@split;
        }else{
            #not part of region;
            processRegion(\@contiguous);
            @contiguous = ();
            push @contiguous, \@split;
        }
        $prev_chr = $split[0];
        $prev_pos = $split[1];
    }
    processRegion(\@contiguous);
    close $FH;
}

##################################################
sub processRegion{
    my $cadds = shift;
    return if not @$cadds;
    my $chr   = $cadds->[0]->[0];
    my $start = $cadds->[0]->[1] - 1;#for BED format
    my $end   = $cadds->[-1]->[1];
    my $coord = "chr$chr\t$start\t$end";
    my $lifted = doLiftOver($coord); 
    my ($chrom, $reg_start, $reg_end)  = parseLiftOverResult($lifted); 
    if ( defined $chrom and ($reg_end - $reg_start) == ($end - $start) ){
        #if we got a liftOver region and it is same size 
        foreach my $c (@$cadds){
            my $pos = $c->[1] - $start + $reg_start ;#plus 1 for BED format
            if ($pos > $reg_end){
                #this shouldn't happen
                warn "$pos is greater than region end for ".join("\t", @$c) . "\n";
                next;
            }
            if (checkRefCoordinate
                (
                    chr => $chrom,
                    pos => $pos,
                    ref => $c->[2]
                )
            ){
                print join("\t", $chrom, $pos, @$c[2..5]) . "\n";
            }else{
                convertIndividually([$c]);
            }
        }
    }else{
        if (@$cadds > 100){#try splitting region in half and see if one half maps
            my $half = int( (scalar@$cadds)/2 );
            processRegion( [ @$cadds[0..$half-1] ] );
            processRegion( [ @$cadds[ $half .. $#{$cadds} ] ] );
            
        }else{
            #region might be split, try each coord individually
            convertIndividually($cadds); 
        }
    }
}

##################################################
sub doLiftOver{
    my $coord = shift;
    my $cmd =  "bash -c '$liftover -minMatch=1.0 "
               ."<(echo $coord) "
               ."$chain /dev/stdout /dev/null "
               ."2> /dev/null'";
    my $result = `$cmd`;
    return $result;
}

##################################################
sub parseLiftOverResult{
    my $result = shift;
    if ($result =~ /(\S+)\t(\d+)\t(\d+)/ ){
        return ($1, $2, $3);
    }
    return;
}

##################################################
sub convertIndividually{
    my $cadds = shift;
    my %lift = (coord => '', result => '');
CADD: foreach my $c (@$cadds){
        my $chr   = $c->[0];
        my $start = $c->[1] - 1;#for BED format
        my $end   = $c->[1];
        my $coord = "chr$chr\t$start\t$end";
        if ($coord ne $lift{coord}){
            $lift{result} = doLiftOver($coord);
            $lift{coord} = $coord;
            chomp $lift{result} ;
        }
        #print "$coord =$lift{result}\n";#debug
        my $reason = '';
        if ($lift{result} and $lift{result} =~ /(\S+)\t(\d+)\t(\d+)/ ){
            my $chr = $1;
            my $pos = $2 + 1;#+1 for BED format
            if (checkRefCoordinate
                (
                    chr => $chr,
                    pos => $pos,
                    ref => $c->[2]
                )
            ){
                print join("\t", $chr, $pos, @$c[2..5]) . "\n";
                next CADD;
            }else{
                $reason = "Non-matching REF allele";
            }
        }else{
            $reason = "liftOver failed";
        }
        my $line = join("\t", @$c);
        warn "Could not convert $line: $reason\n";
    }
}
        
##################################################
sub checkRefCoordinate{
    my %args = @_;
    my $end = $args{pos} + length($args{ref}) -1 ;
    my $seq = $fai->fetch("$args{chr}:$args{pos}-$end");
    if (uc($seq) eq uc($args{ref})){
        return 1;
    }
    return 0;
}
    
##################################################
sub openFile{
    my $f = shift;
    my $FH;
    if ($f =~ /\.gz/){
        $FH = new IO::Uncompress::Gunzip $f, MultiStream => 1 or die 
          "IO::Uncompress::Gunzip failed while opening $f for reading:\n".
          $GunzipError;
    }else{
        open ($FH, $f ) or die "Can't open $f for reading: $!\n";
    }
    return $FH;
}
    

