#!/usr/bin/perl -w

# This script was originally downloaded from the ABI's de novo accessory tools package:
# http://solidsoftwaretools.com/gf/project/denovo/frs/
use Getopt::Long;

$CURRENT_VERSION = "0.1";
$PROGRAM_NAME    = $0;
$PROGRAM_NAME    =~ s|.*/||;

GetOptions('f|input_fasta=s'                 => \$inputFastaFile,
	   'q|input_qv_file=s'               => \$inputQvFile,
	   'o|output=s'                      => \$outputFile,
	   'm|min_median_qv:i'               => \$minMedianQV,
	   'x|num_colors_to_hard_trim:i'     => \$numColorsToHardTrim,
	   'l|min_read_length:i'             => \$minReadLength,
	   'n|num_consec_colors_to_trim:i'   => \$minConsecColorsToTrim,
	   'r|trim_terminal_bad_colors:i'    => \$trimTerminalBadColors,
	   's|min_single_color_QV:i'         => \$minSingleColorQV,
	   'b|max_number_bad_colors:i'       => \$maxNumberBadColors,
	   't|max_bad_colors_in_first_ten:i' => \$maxBadColorsInFirstTen,
	   'c|best_reads_fraction:f'         => \$bestReadsFraction,
	   'p|padded_length:i'               => \$paddedLength,
	   'a|pad_character:s'               => \$padCharacter,
	   'h|help'                          => \$help,
	   'v|version'                       => \$version,
	   );

if ($version) {
    print STDERR "$PROGRAM_NAME, Version $CURRENT_VERSION\n";
    exit;
}

if ($help) {
    &PrintUsage();
    exit;
}

$inputErrors = &CheckInput();
if ( defined($inputErrors) && ($inputErrors ne '') ) {
    &PrintUsage($inputErrors);
    exit;
}

open $CSFASTA_IN,  '<', $inputFastaFile  or die "Can't open input Fasta file '$inputFastaFile': $!\n";
open $QV_IN,       '<', $inputQvFile     or die "Can't open input QV file '$inputQvFile': $!\n";
open $CSFASTA_OUT, '>', $outputFile      or die "Can't open output file '$outputFile': $!\n";

# Set QV file to first data line to synchronize it with the csfasta file
$firstDataLine = 0;
while (<$QV_IN>) {
    last unless /^#/;
    $firstDataLine++;
}
close $QV_IN;
open ($QV_IN, $inputQvFile);
for ($j=0; $j<$firstDataLine; $j++) {
    $line = <$QV_IN>;
}

# Process reads

if (defined($bestReadsFraction) && ($bestReadsFraction != 0)) {
    &GetBestReads($CSFASTA_IN, $QV_IN, $CSFASTA_OUT, $bestReadsFraction);
} else {
    &FilterReads($CSFASTA_IN, $QV_IN, $CSFASTA_OUT);
}
close CSFASTA_IN;
close QV_IN;
close CSFASTA_OUT;


### SUBROUTINES ###


sub GetBestReads {
    my ($CSFASTA_IN, $QV_IN, $CSFASTA_OUT, $bestReadsFraction) = @_;
    my (%medianQvOfReadId);
    my @readQV;
    my $readQV;
    my $readID;
    my $readCount = $filteredReadCount = 0;

    print STDERR "\n!!! Retrieving top $bestReadsFraction of reads by median QV.  Other filtering is disabled.\n\n";

    while (<$QV_IN>) {
	chomp;
	if (/^>/) {
	    if (defined($readQV) && ($readQV ne '') ) {
		@readQV = split(/\s+/, $readQV);
		$medianQvOfReadId{$readID} = &CalculateMedianQV(@readQV);
	    }
	    $readID = $_;
	    $readQV = '';
	} else {
	    $readQV .= $_;
	}
    }

    @sortedReads = sort { $medianQvOfReadId{$b} <=> $medianQvOfReadId{$a} } keys %medianQvOfReadId;
    $lastReadIndex = $bestReadsFraction * scalar(@sortedReads);
    for ($i=0; $i < $lastReadIndex; $i++) {
	$readsToRetrieve{$sortedReads[$i]}++;
    }
    
    while (<$CSFASTA_IN>) {
	chomp;
	next if (/^#/);
	if (/^>/) {
	    if (defined($readSeq) && defined($readsToRetrieve{$readID})) {
		print $CSFASTA_OUT "$readID\n$readSeq\n";
		$filteredReadCount++;
	    }
	    $readID = $_;
	    $readSeq = '';
	    $readCount++;
	} else {
	    $readSeq .= $_;
	}
    }

    # Handle the last record
    if (defined($readSeq) && defined($readsToRetrieve{$readID})) {
	print $CSFASTA_OUT "$readID\n$readSeq\n";
	$filteredReadCount++;
    }
    print STDERR "$readCount reads processed, $filteredReadCount reads remain after filtering\n";
}


sub FilterReads {
    my ($CSFASTA_IN, $QV_IN, $CSFASTA_OUT) = @_;
    my $readCount = $filteredReadCount = 0;
    
    while (<$CSFASTA_IN>) {
	$line = $_;
	chomp($line);
	next if ($line =~ /\#/);
	if ($line =~ /^>/) {
	    if (defined($readSeq) && $readSeq ne '') {
		$qvReadID = <$QV_IN>;
		chomp($qvReadID);
		unless ($qvReadID eq $readID) {
		    die "QV read ID does not equal color read ID\n";
		}
		$readQV = <$QV_IN>;
		chomp($readQV);
		$trimmedRead = &TrimRead($readSeq, $readQV);
		if ($trimmedRead ne '') {
		    if (defined($paddedLength) && ($paddedLength > length($trimmedRead) - 1) ) {
			print $CSFASTA_OUT "$readID\n" . &PadRead($trimmedRead, $paddedLength, $padCharacter) . "\n";
		    } else {
			print $CSFASTA_OUT "$readID\n" . $trimmedRead . "\n";
		    }
		    $filteredReadCount++;
		}
	    }
	    $readCount++;
	    $readID = $line;
	    $readSeq = '';
	} else {
	    $readSeq .= $line;
	}
    }
    
    # Handle the last record
    if (defined($readSeq) && $readSeq ne '') {
	$qvReadID = <$QV_IN>;
	chomp($qvReadID);
	unless ($qvReadID eq $readID) {
	    die "QV read ID does not equal color read ID\n";
	}
	$readQV = <$QV_IN>;
	chomp($readQV);
	$trimmedRead = &TrimRead($readSeq, $readQV);
	if ($trimmedRead ne '') {
	    if (defined($paddedLength) && ($paddedLength > length($trimmedRead) - 1) ) {
		print $CSFASTA_OUT "$readID\n" . &PadRead($trimmedRead, $paddedLength, $padCharacter) . "\n";
	    } else {
		print $CSFASTA_OUT "$readID\n" . $trimmedRead . "\n";
	    }
	    $filteredReadCount++;
	}
    }
    print STDERR "$readCount reads processed, $filteredReadCount reads remain after filtering\n";
}


sub TrimRead {
    my ($readSeq, $readQV) = @_;
    my @readQV = split(/\s+/, $readQV);
    if (defined($numColorsToHardTrim)) {
	if ($numColorsToHardTrim < length($readSeq)) {
	    my $lastIndex = length($readSeq) - 1 - $numColorsToHardTrim;
	    $readSeq = substr($readSeq, 0, $lastIndex + 1); # account for first base
	    @readQV  = @readQV[0..$lastIndex - 1];
	} else {
	    die "num_colors_to_hard_trim value of $numColorsToHardTrim is greater than read length of " . 
		(length($readSeq) - 1) . "\n";
	}
    }
    my $readLength = scalar(@readQV);
    my $medianQV = &CalculateMedianQV(@readQV);
    if ($medianQV < $minMedianQV) {
	return('');		# discard the read if median QV is below the cutoff
    }

    my $numConsecutiveColorsBelowCutoff = 0;
    my $numTotalColorsBelowCutoff       = 0;
    my $trimmedRead                     = '';
    
    for (my $i=0; $i<$readLength; $i++) {
	if ($readQV[$i] < $minSingleColorQV) {
	    $numTotalColorsBelowCutoff++;
	    if ($numTotalColorsBelowCutoff > $maxNumberBadColors) {
		return('');	# discard if number of poor quality colors is above cutoff
	    }
	    if ( ($i < 10) && ($numTotalColorsBelowCutoff > $maxBadColorsInFirstTen) ) {
		return('');	# discard if number of poor quality colors in first ten is above cutoff
	    }
	    $numConsecutiveColorsBelowCutoff++;
	    if ( ($numConsecutiveColorsBelowCutoff >= $minConsecColorsToTrim) && ($trimmedRead eq '') ) {
		$trimmedRead = substr($readSeq, 0, ($i - $minConsecColorsToTrim) + 2);  # everything up to this point
	    }
	} else {
	    $numConsecutiveColorsBelowCutoff = 0;
	}
    }
    if ($trimmedRead eq '') {	# not trimmed
	if ($trimTerminalBadColors) {
	    $trimmedRead = &TrimTerminalBadColors($minSingleColorQV, $readSeq, @readQV);
	    if (length($trimmedRead) >= $minReadLength) {
		return ($trimmedRead);
	    } else {
		return('');
	    }
	} else {
	    return ($readSeq);
	}
    } elsif (length($trimmedRead) >= $minReadLength) {
	if ($trimTerminalBadColors) {
	    $trimmedRead = &TrimTerminalBadColors($minSingleColorQV, $trimmedRead, @readQV[0..length($trimmedRead)-2] );
	    if (length($trimmedRead) >= $minReadLength) {
		return ($trimmedRead);
	    } else {
		return('');
	    }
	} else {
	    return ($trimmedRead);
	}
    } else {
	return('');
    }
}


sub PadRead {
    my ($trimmedRead, $paddedLength, $padCharacter) = @_;
    
    my $padsToAdd = $paddedLength - length($trimmedRead) + 1; # account for the first base of the read
    for (my $i = 0; $i < $padsToAdd; $i++) {
	$trimmedRead .= $padCharacter;
    }
    return $trimmedRead;
}


sub TrimTerminalBadColors {
    my ($minSingleColorQV, $readSeq, @readQV) = @_;
    my $i;

    if (length($readSeq) < 2) {
	return('');
    } else {
	for ($i = length($readSeq); $i >= 0; $i--) {
	    last if ($readQV[$i-2] >= $minSingleColorQV);
	}
	return ( substr($readSeq, 0, $i) );
    }
}


sub CalculateMedianQV {
    my @qvArray = @_;
    my @sortedQvArray = sort {$a<=>$b} (@qvArray);
    my $medianQV;
    my $arrayLength = scalar(@qvArray);

    if ($arrayLength % 2 == 0) {
	$medianQV = ($sortedQvArray[$arrayLength/2 - 1] + $sortedQvArray[($arrayLength/2)]) / 2;
    } else {
	$medianQV = $sortedQvArray[int($arrayLength/2)];
    }
    return($medianQV);
}


sub CheckInput {
    my @errors;

    if (!defined($inputQvFile) || ($inputQvFile eq '') ) {
	push(@errors, 'Input QV file not defined');
    } elsif ( !(-e $inputQvFile)) {
	push(@errors, "Input QV file '$inputQvFile' does not exist");
    }
	
    if (!defined($inputFastaFile) || ($inputFastaFile eq '') ) {
	push(@errors, 'Input Fasta file not defined');
    } elsif (!(-e $inputFastaFile)) {
	push(@errors, "Input Fasta file '$inputFastaFile' does not exist");
    }

    if (!defined($outputFile) || ($outputFile eq '')) {
	push(@errors, 'Output file must be defined');
    }

    if (defined($minMedianQV)) {
	if ($minMedianQV < 0) {
	    push (@errors, 'min_median_qv value must be a non-negative integer');
	}
    } else {
	$minMedianQV = 1;
    }

    if (defined($minReadLength)) {
	if ($minReadLength < 0) {
	    push (@errors, 'min_read_length value must be a non-negative integer');
	}
    } else {
	$minReadLength = 1;
    }

    if (defined($minConsecColorsToTrim)) {
	if ($minConsecColorsToTrim < 0) {
	    push (@errors, 'num_consec_colors_to_trim value must be a non-negative integer');
	}
    } else {
	$minConsecColorsToTrim = 100;
    }

    if (defined($minSingleColorQV)) {
	if ($minSingleColorQV < 0) {
	    push (@errors, 'min_single_color_QV value must be a non-negative integer');
	}
    } else {
	$minSingleColorQV = 1;
    }

    if (defined($maxNumberBadColors)) {
	if ($maxNumberBadColors < 0) {
	    push (@errors, 'max_number_bad_colors value must be a non-negative integer');
	}
    } else {
	$maxNumberBadColors = 100;
    }

    if (defined($maxBadColorsInFirstTen)) {
	if ($maxBadColorsInFirstTen < 0) {
	    push (@errors, 'max_bad_colors_in_first_ten value must be a non-negative integer');
	}
    } else {
	$maxBadColorsInFirstTen = 10;
    }

    if (defined($numColorsToHardTrim)) {
	if ($numColorsToHardTrim < 0) {
	    push (@errors, 'num_colors_to_hard_trim value must be a non-negative integer');
	}
    } else {
	$numColorsToHardTrim = 0;
    }

    if (defined($bestReadsFraction)) {
	if ( ($bestReadsFraction <= 0) || ($bestReadsFraction > 1) ) {
	    push (@errors, 'best_reads_fraction value must be between 0 and 1');
	}
    }

    if (defined($trimTerminalBadColors)) {
	if (! ( ($trimTerminalBadColors == 0) || ($trimTerminalBadColors == 1) ) ) {
	    push (@errors, 'trim_terminal_bad_colors value must be 0 (OFF) or 1 (ON)');
	}
    } else {
	$trimTerminalBadColors = 1;
    }

    if (defined($paddedLength)) {
	if ($paddedLength < 1) {
	    push (@errors, '--padded_length must be a non-negative integer');
	} else {
	    if (defined($padCharacter)) {
		unless ($padCharacter =~ /^\S$/) { # single character
		    push (@errors, '--pad_character must be a single character');
		}
	    } else {
		$padCharacter = '.';
	    }
	}
    } else {
	if (defined($padCharacter)) {
	    push (@errors, '-padded_length parameter must be specified if --pad_character is defined');
	}
    }

    return join("\n", @errors);
}


sub PrintUsage {
    my $errors = shift;

    if (defined($errors)) {
	print STDERR "\n$errors\n";
    }

    print STDERR <<'END';

Usage: csfasta_quality_filter.pl [options] -f <csfasta_file> -q <qv_file> -o <output_file>

Options:
    -m or --min_median_qv <int>               : Reads with median QV below this value 
                                                will be removed from the data set (default = 1)
    -l or --min_read_length <int>             : Reads shorter than this value after trimming 
                                                will be removed from the data set (default = 1)
    -x or --num_colors_to_hard_trim <int>     : This number of colors will be trimmed from every read
                                                before any other processing (default = 0)
    -n or --num_consec_colors_to_trim <int>   : Reads will be trimmed at the beginning of the first
                                                consecutive run of this many colors with QV values
                                                below the minimum single color QV (default = 100)
    -b or --max_number_bad_colors             : Reads with more than this number of low quality colors 
                                                (below the minimum single base QV) will be removed from
						the dataset (default = 100)
    -r or --trim_terminal_bad_colors [0|1]    : Any continuous run of bad colors at the end of the 
                                                read will be trimmed (default = 1 [ON])
    -t or --max_bad_colors_in_first_ten <int> : Reads with more than this number of low quality colors 
                                                (below the minimum single base QV) among the first ten
						colors will be removed from the dataset (default = 10)
    -s or --min_single_color_QV <int>         : Minimum QV value for a single color.  Used in 
                                                conjunction with -n, -b, or -t options.  Consecutive 
						terminal colors that are lower than this value will also
						be trimmed.  Use of a very low value (e.g., the default) 
						may result in a lack of trimming or filtering (default = 1).
    -c or --best_reads_fraction <float>       : Reads are sorted by median QV, and this fraction of
                                                the best reads are returned.  Value must be between 
						0 and 1.  When this mode is used, all other filters 
						are turned off.  This mode is off by default.
    -p or --padded_length <int>               : Length to which each read should be padded.  Reads with 
                                                at least this many bases after filtering will not be padded.
						If this parameter is not used, padding is turned off 
						(default is OFF)
    -a or --pad_character <char>              : Character used to pad reads when padding is turned on
                                                (default = .)
    -v or --version                           : Prints version of the program
    -h or --help                              : Prints this usage summary

Note that defaults values are not stringent and  will result in very little filtering, so at least 
one filtering parameter (-x, -m, -l, -n, -b, -t, or -c) should be specified. When multiple parameters
are specified, they are applied in the following order:

1. num_colors_to_hard_trim
2. min_median_qv
3. max_bad_colors_in_first_ten
4. max_number_bad_colors
5. num_consec_colors_to_trim
6. trim_terminal_bad_colors
7. min_read_length

If -c option is used, no other filtering or trimming takes place.

					
END
    return;
}
