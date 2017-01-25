#! /usr/bin/perl

#  Copyright 2010 Brian Ondov
# 
#  This file is part of SOCS.
# 
#  SOCS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  SOCS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with SOCS.  If not, see <http://www.gnu.org/licenses/>.

# Cheers to RPCM bioinformatics lab team for this script:
# http://rcpcm.org/nauchnye-issledovanija/otdel-molekuljarnoj-biologii-i-genetiki/laboratorija-bioinformatiki/

use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;

my $csfastaFileName;
my $qualFileName;
my $outputDir;
my $qualityThreshold;

GetOptions('f|input_fasta=s'                 => \$csfastaFileName,
	   'q|input_qv_file=s'               => \$qualFileName,
	   'o|output=s'                      => \$outputDir,
	   't|max_bad_colors_in_first_ten=i' => \$qualityThreshold,
	   );

my $readsTotal = 0;
my $readsFiltered = 0;
my $missingQualities = 0;
my $extraQualities = 0;
my $lineQualityTag;

my ($name, $dir, $ext) = fileparse($csfastaFileName, ".csfasta");
my $outCsfastaFile = File::Spec->catfile($outputDir, $name . "." .  $qualityThreshold . ".filtered" . $ext);

my ($name, $dir, $ext) = fileparse($qualFileName, ".qual");
my $outQualFile = File::Spec->catfile($outputDir, $name . "." .  $qualityThreshold . ".filtered" . $ext);

open INFILE_CSFASTA, "<$csfastaFileName" or die $!;
open INFILE_QUAL, "<$qualFileName" or die $!;

open OUTFILE_CSFASTA, ">$outCsfastaFile" or die $!;
open OUTFILE_QUAL, ">$outQualFile" or die $!;

# eat quality file comments
#
while ($lineQualityTag = <INFILE_QUAL>)
{
	if ((substr $lineQualityTag, 0, 1) ne "#")
	{
		last;
	}
}

while (my $line = <INFILE_CSFASTA>)
{
	if ((substr $line, 0, 1) eq ">")
	{
		my $lineColors = <INFILE_CSFASTA>;
		
		while
		(
			$lineQualityTag ne "" &&
			$lineQualityTag lt $line
		)
		{
			my $lineQualities = <INFILE_QUAL>;
			$lineQualityTag = <INFILE_QUAL>;
			
			$extraQualities++;
		}
		
		if ($lineQualityTag eq $line)
		{
			my $lineQualities = <INFILE_QUAL>;
			$lineQualityTag = <INFILE_QUAL>;
			
			if
			(
				$lineColors =~ /^[ACGT][0123]+\n/ &&
				$lineQualities !~ /-/ && # no negative quality values allowed
				($qualityThreshold == 0 || averageQuality($lineQualities) >= $qualityThreshold)
			)
			{
				print OUTFILE_CSFASTA $line;
				print OUTFILE_CSFASTA $lineColors;
				print OUTFILE_QUAL $line;
				print OUTFILE_QUAL $lineQualities;
			}
			else
			{
				$readsFiltered++;
				
	#			print $lineQualities;
	#			print averageQuality($lineQualities);
	#			print "\n\n";
			}
		}
		else
		{
			$readsFiltered++;
			$missingQualities++;
		}
		
		$readsTotal++;
	}
}

close INFILE_CSFASTA;
close INFILE_QUAL;

close OUTFILE_CSFASTA;
close OUTFILE_QUAL;

if ($missingQualities > 0)
{
	print "WARNING: $missingQualities reads filtered due to missing quality scores.\n";
}

if ($extraQualities > 0)
{
	print "WARNING: $extraQualities extra quality tags ignored.\n";
}

print "Filtered $readsFiltered of $readsTotal reads.\n";

sub averageQuality
{
	my ($string) = @_;
	
	my @qualities = split / /, $string;
	
	pop @qualities; # the space on the end of each quality list will cause an extra member
	
	my $total = 0;
	
	map {$total += $_} @qualities;
	
	return $total / @qualities;
}
