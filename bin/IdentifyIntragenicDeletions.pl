#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# Add read-level breakpoint/deletion annotations directly to split_reads.tsv.
# The input rows are never filtered: records that do not represent a valid
# deletion are retained with Breakpoint_Status=Filtered and Reason=Other (or
# the two explicitly classified mismatch reasons).

my ($input, $out, $help);
GetOptions(
    'Input=s' => \$input,
    'Out=s'   => \$out,
    'Help'    => \$help,
) or die usage();
if ($help) { print usage(); exit 0; }
die usage() unless defined $input && defined $out;

$input = abs_path($input) if -e $input;
die "[ERROR] Input file not found: $input\n" unless defined $input && -s $input;
die "[ERROR] Input path must be absolute\n" unless $input =~ m{^/};
die "[ERROR] Output path must be absolute\n" unless $out =~ m{^/};
make_path(dirname($out)) unless -d dirname($out);

open my $IN, '<', $input or die "[ERROR] Cannot open $input: $!\n";
my $header_line = <$IN>;
die "[ERROR] Empty input: $input\n" unless defined $header_line;
chomp $header_line;
$header_line =~ s/\r$//;
my @input_header = split /\t/, $header_line, -1;

my @derived_columns = qw(
    Primary_Segment_Role SA_Segment_Role
    Left_Breakpoint Right_Breakpoint
    Chrom Start End Size_bp Strand
    Junction_Type Breakpoint_Status Breakpoint_Reason
);
my %derived = map { $_ => 1 } @derived_columns;
my @base_indices = grep { !$derived{$input_header[$_]} } 0 .. $#input_header;
my @base_header = @input_header[@base_indices];
my %idx;
for my $i (0 .. $#base_header) {
    die "[ERROR] Duplicate input column '$base_header[$i]'\n" if exists $idx{$base_header[$i]};
    $idx{$base_header[$i]} = $i;
}

my @required = qw(
    SampleID Read_Name Flag
    Primary_Chrom Primary_Start Primary_End Primary_CIGAR Primary_MAPQ Primary_Strand
    SA_Chrom SA_Start SA_End SA_CIGAR SA_MAPQ SA_Strand SA_NM
);
for my $col (@required) {
    die "[ERROR] Missing input column: $col\n" unless exists $idx{$col};
}

# Always write to a temporary file first. This supports --Input and --Out
# pointing to the same split_reads.tsv without truncating the input.
my $tmp_out = "$out.tmp.$$";
unlink $tmp_out if -e $tmp_out;
open my $OUT, '>', $tmp_out or die "[ERROR] Cannot write temporary output $tmp_out: $!\n";
my @out_header = (@base_header, @derived_columns);
print {$OUT} join("\t", @out_header), "\n";

my %stats = (
    input_rows             => 0,
    Valid_Deletion_Gap     => 0,
    Different_Chromosome   => 0,
    Different_Strand       => 0,
    Other                  => 0,
);
my $line_no = 1;
while (my $line = <$IN>) {
    $line_no++;
    chomp $line;
    $line =~ s/\r$//;
    next if $line eq '';
    $stats{input_rows}++;

    my @raw_values = split /\t/, $line, -1;
    die "[ERROR] Column count mismatch at $input:$line_no (expected "
        . scalar(@input_header) . ", observed " . scalar(@raw_values) . ")\n"
        unless @raw_values == @input_header;
    my @values = @raw_values[@base_indices];

    my %record;
    @record{@base_header} = map { defined $_ && $_ ne '' ? $_ : 'NA' } @values;

    my $annotation = infer_breakpoint(\%record);
    $stats{$annotation->{Breakpoint_Reason}}++;

    my @annotation_values = map {
        defined $annotation->{$_} && $annotation->{$_} ne ''
            ? $annotation->{$_} : 'NA'
    } @derived_columns;

    print {$OUT} join("\t", @values, @annotation_values), "\n";
}
close $IN or die "[ERROR] Failed while reading $input\n";
close $OUT or die "[ERROR] Cannot close temporary output $tmp_out: $!\n";

rename $tmp_out, $out
    or die "[ERROR] Cannot replace output $out with $tmp_out: $!\n";

for my $key (qw/input_rows Valid_Deletion_Gap Different_Chromosome Different_Strand Other/) {
    print STDERR "[INFO] $key: $stats{$key}\n";
}
print STDERR "[INFO] Annotated split-read table: $out\n";
exit 0;

# ------------------------------------------------------------------
# Breakpoint inference
# ------------------------------------------------------------------

sub infer_breakpoint {
    my ($r) = @_;

    my %result = (
        Primary_Segment_Role => 'NA',
        SA_Segment_Role      => 'NA',
        Left_Breakpoint      => 'NA',
        Right_Breakpoint     => 'NA',
        Chrom                => 'NA',
        Start                => 'NA',
        End                  => 'NA',
        Size_bp              => 'NA',
        Strand               => 'NA',
        Junction_Type        => 'NA',
        Breakpoint_Status    => 'Filtered',
        Breakpoint_Reason    => 'Other',
    );

    my ($p_chr, $sa_chr) = @{$r}{qw/Primary_Chrom SA_Chrom/};
    return \%result unless valid_text($p_chr) && valid_text($sa_chr);
    if ($p_chr ne $sa_chr) {
        $result{Breakpoint_Reason} = 'Different_Chromosome';
        return \%result;
    }

    my ($p_strand, $sa_strand) = @{$r}{qw/Primary_Strand SA_Strand/};
    return \%result unless defined $p_strand && defined $sa_strand
        && $p_strand =~ /^[+-]$/ && $sa_strand =~ /^[+-]$/;
    if ($p_strand ne $sa_strand) {
        $result{Breakpoint_Reason} = 'Different_Strand';
        return \%result;
    }

    my ($p_start, $p_end, $sa_start, $sa_end) =
        @{$r}{qw/Primary_Start Primary_End SA_Start SA_End/};
    return \%result unless is_positive_integer($p_start)
        && is_positive_integer($p_end)
        && is_positive_integer($sa_start)
        && is_positive_integer($sa_end)
        && $p_end >= $p_start && $sa_end >= $sa_start;

    my $p_cigar  = parse_cigar($r->{Primary_CIGAR});
    my $sa_cigar = parse_cigar($r->{SA_CIGAR});
    return \%result unless $p_cigar && $sa_cigar;

    # Check that the reported alignment intervals agree with their CIGARs.
    return \%result unless $p_cigar->{reference_length} == $p_end - $p_start + 1
        && $sa_cigar->{reference_length} == $sa_end - $sa_start + 1;

    # Both SA records must describe the same query length. Hard clips are
    # included because they represent omitted bases from the original read.
    return \%result unless $p_cigar->{query_length} == $sa_cigar->{query_length};

    my $p_clip_side  = terminal_softclip_side($p_cigar->{ops});
    my $sa_clip_side = terminal_softclip_side($sa_cigar->{ops});
    return \%result unless defined $p_clip_side && defined $sa_clip_side;

    my ($left_end, $right_start);
    if ($p_clip_side eq 'trailing' && $sa_clip_side eq 'leading') {
        # Primary and SA are consecutive segments in this order.
        return \%result unless $p_end < $sa_start;
        $result{Primary_Segment_Role} = 'Left';
        $result{SA_Segment_Role} = 'Right';
        ($left_end, $right_start) = ($p_end, $sa_start);
    }
    elsif ($p_clip_side eq 'leading' && $sa_clip_side eq 'trailing') {
        # SA and Primary are consecutive segments in this order.
        return \%result unless $sa_end < $p_start;
        $result{Primary_Segment_Role} = 'Right';
        $result{SA_Segment_Role} = 'Left';
        ($left_end, $right_start) = ($sa_end, $p_start);
    }
    else {
        return \%result;
    }

    $result{Left_Breakpoint}  = $left_end;
    $result{Right_Breakpoint} = $right_start;

    my $del_start = $left_end + 1;
    my $del_end   = $right_start - 1;
    return \%result unless $del_end >= $del_start;

    $result{Chrom}             = $p_chr;
    $result{Start}             = $del_start;
    $result{End}               = $del_end;
    $result{Size_bp}           = $del_end - $del_start + 1;
    $result{Strand}            = $p_strand;
    $result{Junction_Type}     = 'Deletion_Gap';
    $result{Breakpoint_Status} = 'PASS';
    $result{Breakpoint_Reason} = 'Valid_Deletion_Gap';

    return \%result;
}

sub parse_cigar {
    my ($cigar) = @_;
    return unless defined $cigar && $cigar ne '' && $cigar ne 'NA' && $cigar ne '*';

    my @ops;
    my $consumed = 0;
    while ($cigar =~ /\G([1-9]\d*)([MIDNSHP=X])/gc) {
        my ($length, $op) = ($1, $2);
        push @ops, [$length + 0, $op];
        $consumed = pos($cigar);
    }
    return unless @ops && $consumed == length($cigar);

    # Hard clips are only valid at the two ends of a CIGAR.
    for my $i (0 .. $#ops) {
        next unless $ops[$i][1] eq 'H';
        return unless $i == 0 || $i == $#ops;
    }

    # Soft clips may only occur at the two alignment ends (after optional H).
    my $first_non_h = 0;
    $first_non_h++ while $first_non_h <= $#ops && $ops[$first_non_h][1] eq 'H';
    my $last_non_h = $#ops;
    $last_non_h-- while $last_non_h >= 0 && $ops[$last_non_h][1] eq 'H';
    for my $i (0 .. $#ops) {
        next unless $ops[$i][1] eq 'S';
        return unless $i == $first_non_h || $i == $last_non_h;
    }

    my ($reference_length, $query_length) = (0, 0);
    my $has_aligned_bases = 0;
    for my $op (@ops) {
        my ($length, $code) = @$op;
        $reference_length += $length if $code =~ /^[MDN=X]$/;
        $query_length += $length if $code =~ /^[MIS=XH]$/;
        $has_aligned_bases = 1 if $code =~ /^[M=X]$/;
    }
    return unless $has_aligned_bases && $reference_length > 0 && $query_length > 0;

    return {
        ops              => \@ops,
        reference_length => $reference_length,
        query_length     => $query_length,
    };
}

sub terminal_softclip_side {
    my ($ops_ref) = @_;
    my @ops = @$ops_ref;

    shift @ops while @ops && $ops[0][1] eq 'H';
    pop @ops while @ops && $ops[-1][1] eq 'H';
    return unless @ops;

    my $leading  = $ops[0][1] eq 'S' ? 1 : 0;
    my $trailing = $ops[-1][1] eq 'S' ? 1 : 0;
    return if $leading == $trailing; # neither end, or both ends, are soft-clipped

    return $leading ? 'leading' : 'trailing';
}

sub valid_text {
    return defined $_[0] && $_[0] ne '' && $_[0] ne 'NA' && $_[0] ne '*';
}

sub is_positive_integer {
    return defined $_[0] && $_[0] =~ /^\d+$/ && $_[0] > 0;
}

sub usage {
    return "Usage: perl IdentifyIntragenicDeletions.pl --Input split_reads.tsv --Out split_reads.tsv\n";
}
