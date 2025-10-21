#!/usr/bin/perl
# This script reads a blast outfmt6 file with an extra column at the beginning
# indicating the sample name from stdin to create a compact binary file.
# The script works optimally when the input file is sorted/clustered by sample
# name and then by query sequence name. It writes four files in the current
# directory: out.data (binary data), samples.txt (int->sample name mapping),
# genomes.txt (int->subject id mapping).
# 
# Notes:
# - Reads are assigned an even int id because the last bit indicates whether
#   the read is read1 or read2.
# - Percent identity is stored as an int normalized to the range of uint16
# - E-value is stored in phred scale

use Digest::MD5 qw(md5);
use Encode qw(encode_utf8);

open(ofx, '>', "out.data") or die $!;

my %chroms = ();
my %refs = ();
my %samps = ();

my $q_idx = 0;
open(FAIQ,$ARGV[0]) or die "could not open query index file: $!";
while(<FAIQ>){chomp; my @a=split("\t"); $chroms{$a[0]} = $q_idx; $q_idx+=1}
$chroms{"*"} = -1;

my $s_idx = 0;
open(FAIS,$ARGV[1]) or die "could not open subject index file: $!";
while(<FAIS>){chomp; my @a=split("\t"); $refs{$a[0]} = $s_idx; $s_idx+=1}

my $smp_idx = 0;
open(FAIS,$ARGV[2]) or die "could not open samples file: $!";
while(<FAIS>){chomp; my @a=split("\t"); $samps{$a[0]} = $smp_idx; $smp_idx+=1}

while(<stdin>) {
    chomp;
    my @F = split "\t";

    # record name expected in fmt:
    # >(qname)-(mapper_flag):(mapq):(rname):(rpos):(mate_rname):(mate_rpos):(qlen)
    my @read_parts = split "-", $F[1];
    my @map_info = split ":", $read_parts[1];

    print ofx pack("vQvCvVvVCvCCCCCCVVsv",
        $samps{$F[0]}, # sample number, v - uint16 LE
        unpack("Q", pack(
            "a8", md5(encode_utf8($read_parts[0]))
        )) & (~1), # read name hash, Q - uint64 LE
        $map_info[0], # read flags, v - uint16 LE
        $map_info[1], # read mapq, C - uint8
        $chroms{$map_info[2]}, # read chrom, v - uint16 LE
        $map_info[3], # read pos, V – uint32 LE
        $chroms{$map_info[4]}, # read mate chrom, v - uint16 LE
        $map_info[5], # read mate pos,  V - uint32 LE
        $map_info[6], # read length, C - uint8
        $refs{$F[2]}, # subject sequence num, v - uint16 LE
        $F[3]/100*255, # percent identity in 8 bits, C - uint8
        $F[4], # alignment length, C - uint8
        $F[5], # mismatches, C - uint8
        $F[6], # gap opens, C - uint8
        $F[7], # query sequence start, C – uint8
        $F[8], # query sequence end, C – uint8
        $F[9], # subject sequence start, V – uint32 LE
        $F[10], # subject sequence end, V – uint32 LE
        -10*(log($F[11])/log(10)), # phred-scaled evalue, s - int16 LE
        $F[12]*10, # 10x bitscore, v - uint16 LE
    );
    $refs{$F[2]}

}

close(ofx);
