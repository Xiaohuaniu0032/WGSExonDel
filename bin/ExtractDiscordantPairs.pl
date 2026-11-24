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

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $bam
    && defined $sample && defined $out;

$config = abs_path($config) if -e $config;
$bam    = abs_path($bam) if -e $bam;

die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;
die "[ERROR] BAM file not found: $bam\n"
    unless defined $bam && -s $bam;
die "[ERROR] Config and BAM paths must be absolute\n"
    unless $config =~ m{^/} && $bam =~ m{^/};
die "[ERROR] Invalid Sample ID: $sample\n"
    unless $sample =~ /^[A-Za-z0-9_.-]+$/;
die "[ERROR] Output path must be absolute\n"
    unless $out =~ m{^/};

check_bam_index($bam);

my %conf = read_config($config);
my $samtools = required_config(\%conf, 'SAMTOOLS');
my $min_mapq = required_config(\%conf, 'MIN_MAPQ');
my $exclude_duplicates = required_config(\%conf, 'EXCLUDE_DUPLICATES');
my $min_insert_size = required_config(\%conf, 'MIN_DISCORDANT_INSERT_SIZE');
my $filter_orientation = required_config(\%conf, 'FILTER_DELETION_ORIENTATION');
my $threads = required_config(\%conf, 'DISCORDANT_READ_THREADS');

validate_nonnegative_integer('MIN_MAPQ', $min_mapq);
validate_bool('EXCLUDE_DUPLICATES', $exclude_duplicates);
validate_positive_integer('MIN_DISCORDANT_INSERT_SIZE', $min_insert_size);
validate_bool('FILTER_DELETION_ORIENTATION', $filter_orientation);
validate_positive_integer('DISCORDANT_READ_THREADS', $threads);
check_executable($samtools, 'SAMTOOLS');

make_path(dirname($out)) unless -d dirname($out);
open my $OUT, '>', $out
    or die "[ERROR] Cannot write $out: $!\n";

my @output_header = qw(
    Pair_ID SampleID Read_Name Flag Chrom Read_Start Read_End Mate_Chrom
    Mate_Start TLEN Insert_Size MAPQ CIGAR Read_Strand Mate_Strand
    Pair_Orientation Left_Breakpoint Right_Breakpoint Start End Size_bp
);
print {$OUT} join("\t", @output_header), "\n"
    or die "[ERROR] Cannot write header to $out: $!\n";

my $exclude_flag = 4 + 8 + 256 + 512 + 2048;
$exclude_flag += 1024 if $exclude_duplicates;

my @cmd = (
    $samtools, 'view', '-@', $threads,
    '-f', 1,
    '-F', $exclude_flag,
    '-q', $min_mapq,
    $bam,
);

print STDERR "[INFO] Extracting and filtering Discordant read pairs for $sample\n";
print STDERR "[CMD] ", join(' ', map { shell_quote($_) } @cmd), "\n";

open my $SAM, '-|', @cmd
    or die "[ERROR] Cannot run samtools view: $!\n";

my %stats = (
    records          => 0,
    invalid          => 0,
    filtered_flag    => 0,
    same_chromosome  => 0,
    positive_tlen    => 0,
    insert_size      => 0,
    orientation      => 0,
    valid_interval   => 0,
    duplicate_pair   => 0,
    written          => 0,
);
my %seen_pair;

while (my $line = <$SAM>) {
    chomp $line;
    next if $line eq '';
    $stats{records}++;

    my @fields = split /\t/, $line, -1;
    if (@fields < 11) {
        $stats{invalid}++;
        next;
    }

    my ($read_name, $flag, $chrom, $read_start, $mapq, $cigar,
        $mate_chrom, $mate_start, $tlen) = @fields[0 .. 8];

    unless (defined $flag && $flag =~ /^\d+$/
        && is_paired($flag)
        && !is_unmapped($flag)
        && !is_mate_unmapped($flag)
        && !is_secondary($flag)
        && !is_supplementary($flag)
        && !is_qc_fail($flag)
        && !($exclude_duplicates && is_duplicate($flag))) {
        $stats{filtered_flag}++;
        next;
    }

    unless (defined $read_name && $read_name ne ''
        && defined $chrom && $chrom ne '*'
        && defined $read_start && $read_start =~ /^\d+$/ && $read_start > 0
        && defined $mapq && $mapq =~ /^\d+$/ && $mapq >= $min_mapq
        && valid_cigar($cigar)
        && defined $mate_chrom && $mate_chrom ne '*'
        && defined $mate_start && $mate_start =~ /^\d+$/ && $mate_start > 0
        && defined $tlen && $tlen =~ /^-?\d+$/) {
        $stats{invalid}++;
        next;
    }

    my $reference_length = cigar_ref_length($cigar);
    unless ($reference_length > 0) {
        $stats{invalid}++;
        next;
    }

    my $resolved_mate_chrom = $mate_chrom eq '=' ? $chrom : $mate_chrom;
    next unless $chrom eq $resolved_mate_chrom;
    $stats{same_chromosome}++;

    # Retain only the positive-TLEN alignment to represent each fragment once.
    next unless $tlen > 0;
    $stats{positive_tlen}++;

    my $insert_size = abs($tlen);
    next unless $insert_size >= $min_insert_size;
    $stats{insert_size}++;

    my $read_strand = is_reverse($flag) ? '-' : '+';
    my $mate_strand = is_mate_reverse($flag) ? '-' : '+';
    my $orientation = ($read_strand eq '+' ? 'F' : 'R')
        . ($mate_strand eq '+' ? 'F' : 'R');

    if ($filter_orientation) {
        next unless $read_strand eq '+' && $mate_strand eq '-';
    }
    $stats{orientation}++;

    my $read_end = $read_start + $reference_length - 1;
    my $start = $read_end + 1;
    my $end = $mate_start - 1;
    next unless $start >= 1 && $end >= $start;
    $stats{valid_interval}++;

    my $pair_key = join('|',
        $sample, $read_name, $chrom, $read_start, $read_end,
        $resolved_mate_chrom, $mate_start, $insert_size,
    );
    if ($seen_pair{$pair_key}++) {
        $stats{duplicate_pair}++;
        next;
    }

    my $pair_id = sprintf('DR_PAIR_%06d', $stats{written} + 1);
    my @output = (
        $pair_id,
        $sample,
        $read_name,
        $flag,
        $chrom,
        $read_start,
        $read_end,
        $resolved_mate_chrom,
        $mate_start,
        $tlen,
        $insert_size,
        $mapq,
        $cigar,
        $read_strand,
        $mate_strand,
        $orientation,
        $read_end,
        $mate_start,
        $start,
        $end,
        $end - $start + 1,
    );

    print {$OUT} join("\t", @output), "\n"
        or die "[ERROR] Cannot write Discordant pair to $out: $!\n";
    $stats{written}++;
}

close $SAM
    or die "[ERROR] samtools view failed for $bam\n";
close $OUT
    or die "[ERROR] Cannot close $out: $!\n";

for my $key (qw/
    records invalid filtered_flag same_chromosome positive_tlen insert_size
    orientation valid_interval duplicate_pair written
/) {
    print STDERR "[INFO] $key: $stats{$key}\n";
}
print STDERR "[INFO] Filtered Discordant pairs: $out\n";
exit 0;

# ------------------------------------------------------------------
# CIGAR, SAM flag, configuration, and validation helpers
# ------------------------------------------------------------------

sub valid_cigar {
    my ($cigar) = @_;
    return 0 unless defined $cigar && $cigar ne '' && $cigar ne '*';
    my $position = 0;
    while ($cigar =~ /\G(\d+)([MIDNSHP=X])/gc) {
        $position = pos($cigar);
    }
    return $position == length($cigar) && $cigar =~ /\d+[M=X]/;
}

sub cigar_ref_length {
    my ($cigar) = @_;
    my $length = 0;
    while ($cigar =~ /(\d+)([MIDNSHP=X])/g) {
        my ($size, $operation) = ($1, $2);
        $length += $size if $operation =~ /^[MDN=X]$/;
    }
    return $length;
}

sub read_config {
    my ($file) = @_;
    my %config;
    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*(?:#|$)/;
        $line =~ s/\s+#.*$//;
        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;
        die "[ERROR] Duplicate config key $1\n" if exists $config{$1};
        $config{$1} = $2;
    }
    close $fh;
    return %config;
}

sub required_config {
    my ($config_ref, $key) = @_;
    die "[ERROR] Missing config $key\n"
        unless defined $config_ref->{$key} && $config_ref->{$key} ne '';
    return $config_ref->{$key};
}

sub validate_bool {
    my ($name, $value) = @_;
    die "[ERROR] $name must be 0 or 1\n" unless $value =~ /^[01]$/;
}

sub validate_positive_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a positive integer\n"
        unless $value =~ /^\d+$/ && $value > 0;
}

sub validate_nonnegative_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a nonnegative integer\n"
        unless $value =~ /^\d+$/;
}

sub is_paired         { return $_[0] & 0x1   ? 1 : 0; }
sub is_unmapped       { return $_[0] & 0x4   ? 1 : 0; }
sub is_mate_unmapped  { return $_[0] & 0x8   ? 1 : 0; }
sub is_reverse        { return $_[0] & 0x10  ? 1 : 0; }
sub is_mate_reverse   { return $_[0] & 0x20  ? 1 : 0; }
sub is_secondary      { return $_[0] & 0x100 ? 1 : 0; }
sub is_qc_fail        { return $_[0] & 0x200 ? 1 : 0; }
sub is_duplicate      { return $_[0] & 0x400 ? 1 : 0; }
sub is_supplementary  { return $_[0] & 0x800 ? 1 : 0; }

sub check_executable {
    my ($command, $name) = @_;
    if ($command =~ m{/}) {
        die "[ERROR] $name not executable: $command\n" unless -x $command;
    }
    else {
        system("command -v $command >/dev/null 2>&1") == 0
            or die "[ERROR] $name not found: $command\n";
    }
}

sub check_bam_index {
    my ($bam_file) = @_;
    my $index_a = "$bam_file.bai";
    my $index_b = $bam_file;
    $index_b =~ s/\.bam$/.bai/i;
    die "[ERROR] BAM index missing: $bam_file\n"
        unless -s $index_a || -s $index_b;
}

sub shell_quote {
    my ($value) = @_;
    $value =~ s/'/'\\''/g;
    return "'$value'";
}

sub usage {
    return "Usage: perl ExtractDiscordantPairs.pl --Config CONF --Bam BAM --Sample ID --Out discordant_pairs.tsv\n";
}
