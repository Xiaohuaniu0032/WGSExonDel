#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($config, $bam, $sample, $out, $help);
GetOptions(
    'Config=s' => \$config,
    'Bam=s'    => \$bam,
    'Sample=s' => \$sample,
    'Out=s'    => \$out,
    'Help'     => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage() unless defined $config && defined $bam && defined $sample && defined $out;

$config = abs_path($config) if -e $config;
$bam = abs_path($bam) if -e $bam;
die "[ERROR] Config file not found: $config\n" unless defined $config && -s $config;
die "[ERROR] BAM file not found: $bam\n" unless defined $bam && -s $bam;
die "[ERROR] Config and BAM paths must be absolute\n" unless $config =~ m{^/} && $bam =~ m{^/};
die "[ERROR] Invalid Sample ID: $sample\n" unless $sample =~ /^[A-Za-z0-9_.-]+$/;
die "[ERROR] Output path must be absolute: $out\n" unless $out =~ m{^/};
check_bam_index($bam);

my %conf = read_config($config);
my $samtools = required_config(\%conf, 'SAMTOOLS');
my $min_mapq = required_config(\%conf, 'MIN_MAPQ');
my $exclude_duplicates = required_config(\%conf, 'EXCLUDE_DUPLICATES');
my $threads = required_config(\%conf, 'SPLIT_READ_THREADS');
validate_nonnegative_integer('MIN_MAPQ', $min_mapq);
validate_bool('EXCLUDE_DUPLICATES', $exclude_duplicates);
validate_positive_integer('SPLIT_READ_THREADS', $threads);
check_executable($samtools, 'SAMTOOLS');

make_path(dirname($out)) unless -d dirname($out);
open my $OUT, '>', $out or die "[ERROR] Cannot write output: $out: $!\n";
my @header = qw(
    SampleID Read_Name Flag Primary_Chrom Primary_Start Primary_End Primary_CIGAR
    Primary_MAPQ Primary_Strand SA_Chrom SA_Start SA_End SA_CIGAR SA_MAPQ SA_Strand SA_NM
);
print {$OUT} join("\t", @header), "\n";

my $require_flag = 1;
my $exclude_flag = 4 + 256 + 512 + 2048;
$exclude_flag += 1024 if $exclude_duplicates;
my @cmd = ($samtools, 'view', '-@', $threads, '-f', $require_flag, '-F', $exclude_flag, '-q', $min_mapq, $bam);
print STDERR "[INFO] Extracting raw SA alignments for $sample\n";
print STDERR "[CMD] ", join(' ', map { shell_quote($_) } @cmd), "\n";
open my $SAM, '-|', @cmd or die "[ERROR] Cannot run samtools view: $!\n";

my %stats = (sam_records => 0, with_sa => 0, written => 0, no_sa => 0, multi_sa => 0,
             invalid_sam => 0, invalid_sa => 0, low_sa_mapq => 0, excluded_flag => 0);
while (my $line = <$SAM>) {
    chomp $line;
    next if $line eq '';
    $stats{sam_records}++;
    my @f = split /\t/, $line, -1;
    if (@f < 11) { $stats{invalid_sam}++; next; }
    my ($qname, $flag, $rname, $pos, $mapq, $cigar) = @f[0..5];
    unless (defined $flag && $flag =~ /^\d+$/ && is_paired($flag) && !is_unmapped($flag)
            && !is_secondary($flag) && !is_supplementary($flag)) {
        $stats{excluded_flag}++; next;
    }
    if ($exclude_duplicates && is_duplicate($flag)) { $stats{excluded_flag}++; next; }
    unless (defined $qname && $qname ne '' && defined $rname && $rname ne '*' &&
            defined $pos && $pos =~ /^\d+$/ && $pos > 0 &&
            defined $mapq && $mapq =~ /^\d+$/ && $mapq >= $min_mapq &&
            valid_cigar($cigar)) {
        $stats{invalid_sam}++; next;
    }
    my $primary_ref_len = cigar_ref_length($cigar);
    unless ($primary_ref_len > 0) { $stats{invalid_sam}++; next; }

    my ($has_sa, $sa_tag) = parse_sa_tag(@f);
    unless ($has_sa) { $stats{no_sa}++; next; }
    $stats{with_sa}++;
    my @raw_entries = grep { defined $_ && $_ ne '' } split /;/, $sa_tag;
    if (@raw_entries != 1) { $stats{multi_sa}++; next; }
    my $sa = parse_one_sa_entry($raw_entries[0]);
    unless ($sa) { $stats{invalid_sa}++; next; }
    if ($sa->{mapq} < $min_mapq) { $stats{low_sa_mapq}++; next; }
    my $sa_ref_len = cigar_ref_length($sa->{cigar});
    unless ($sa_ref_len > 0) { $stats{invalid_sa}++; next; }

    my @values = (
        $sample, $qname, $flag, $rname, $pos, $pos + $primary_ref_len - 1,
        $cigar, $mapq, (is_reverse($flag) ? '-' : '+'),
        $sa->{chr}, $sa->{pos}, $sa->{pos} + $sa_ref_len - 1,
        $sa->{cigar}, $sa->{mapq}, $sa->{strand}, $sa->{nm},
    );
    print {$OUT} join("\t", map { defined $_ && $_ ne '' ? $_ : 'NA' } @values), "\n";
    $stats{written}++;
}
close $SAM or die "[ERROR] samtools view failed for $bam\n";
close $OUT or die "[ERROR] Cannot close output: $out\n";

for my $key (qw/sam_records with_sa written no_sa multi_sa invalid_sam invalid_sa low_sa_mapq excluded_flag/) {
    print STDERR "[INFO] $key: $stats{$key}\n";
}
print STDERR "[INFO] Raw SA alignment table: $out\n";
exit 0;

sub parse_sa_tag {
    my @fields = @_;
    for my $i (11 .. $#fields) {
        return (1, $1) if defined $fields[$i] && $fields[$i] =~ /^SA:Z:(.*)$/;
    }
    return (0, undef);
}
sub parse_one_sa_entry {
    my ($entry) = @_;
    my @f = split /,/, $entry, -1;
    return unless @f >= 6;
    my ($chr, $pos, $strand, $cigar, $mapq, $nm) = @f[0..5];
    return unless defined $chr && $chr ne '' && defined $pos && $pos =~ /^\d+$/ && $pos > 0
        && defined $strand && $strand =~ /^[+-]$/ && valid_cigar($cigar)
        && defined $mapq && $mapq =~ /^\d+$/;
    $nm = 'NA' unless defined $nm && $nm =~ /^\d+$/;
    return { chr=>$chr, pos=>$pos, strand=>$strand, cigar=>$cigar, mapq=>$mapq, nm=>$nm };
}
sub valid_cigar {
    my ($cigar) = @_;
    return 0 unless defined $cigar && $cigar ne '' && $cigar ne '*';
    my $copy = $cigar;
    my $consumed = 0;
    while ($copy =~ /\G(\d+)([MIDNSHP=X])/gc) { $consumed = pos($copy); }
    return $consumed == length($cigar) && $cigar =~ /\d+[M=X]/ ? 1 : 0;
}
sub cigar_ref_length {
    my ($cigar) = @_;
    my $length = 0;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($n, $op) = ($1, $2);
        $length += $n if $op =~ /^[MDN=X]$/;
    }
    return $length;
}
sub read_config {
    my ($file) = @_;
    my %cfg;
    open my $fh, '<', $file or die "[ERROR] Cannot read config $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line; $line =~ s/\r$//;
        next if $line =~ /^\s*(?:#|$)/;
        $line =~ s/\s+#.*$//;
        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;
        die "[ERROR] Duplicate config key: $1\n" if exists $cfg{$1};
        $cfg{$1} = $2;
    }
    close $fh;
    return %cfg;
}
sub required_config { my ($c,$k)=@_; die "[ERROR] Missing config: $k\n" unless defined $c->{$k} && $c->{$k} ne ''; return $c->{$k}; }
sub validate_bool { my ($n,$v)=@_; die "[ERROR] $n must be 0 or 1\n" unless $v =~ /^[01]$/; }
sub validate_positive_integer { my ($n,$v)=@_; die "[ERROR] $n must be a positive integer\n" unless $v =~ /^\d+$/ && $v > 0; }
sub validate_nonnegative_integer { my ($n,$v)=@_; die "[ERROR] $n must be a non-negative integer\n" unless $v =~ /^\d+$/; }
sub is_paired { return ($_[0] & 0x1) ? 1 : 0; }
sub is_unmapped { return ($_[0] & 0x4) ? 1 : 0; }
sub is_reverse { return ($_[0] & 0x10) ? 1 : 0; }
sub is_secondary { return ($_[0] & 0x100) ? 1 : 0; }
sub is_duplicate { return ($_[0] & 0x400) ? 1 : 0; }
sub is_supplementary { return ($_[0] & 0x800) ? 1 : 0; }
sub check_executable {
    my ($cmd,$name)=@_;
    if ($cmd =~ m{/}) { die "[ERROR] $name is not executable: $cmd\n" unless -x $cmd; }
    else { system("command -v $cmd >/dev/null 2>&1") == 0 or die "[ERROR] $name not found in PATH: $cmd\n"; }
}
sub check_bam_index {
    my ($file)=@_;
    my $a="$file.bai"; my $b=$file; $b =~ s/\.bam$/.bai/i;
    die "[ERROR] BAM index not found for $file\n" unless -s $a || -s $b;
}
sub shell_quote { my ($s)=@_; $s =~ s/'/'\\''/g; return "'$s'"; }
sub usage {
    return "Usage: perl ExtractSASplitReads.pl --Config CONF --Bam BAM --Sample ID --Out split_reads.tsv\n";
}
