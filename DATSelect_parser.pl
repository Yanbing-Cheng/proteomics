#! /usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd 'abs_path';

my ($case_files,$control_files,$protein_file,$normalize,$help);
GetOptions (
	'case=s' => \$case_files,
	'ctrl=s' => \$control_files,
	#'prot=s' => \$protein_file,
	'norm!' => \$normalize,
	'help!'  => \$help,

);

if (defined ($help) || !defined ($case_files) || !defined ($control_files) ) {
	print "Usage:\n";
	print "Example 1: perl  format_patternlab_data.pl -case case/1A,case/1B,case/1C -ctrl control/1A,control/1B,control/1C   -norm\n";
        print "Example 2: perl  format_patternlab_data.pl -case case/1 -ctrl control/1  -norm\n";
        print "Author: Yanbing Cheng, chengyanbing\@gmail.com\n";
	exit;
	
}

my @cases=split(',',$case_files);
my @controls=split(',',$control_files);

my %prot_len;
my %index_id;
my %index_annot;
my %raw_value;
my $index_count=1;

#my $name;
#open (FA,"<$protein_file") or die "can not open file $protein_file";
#while (<FA>) {
#	chomp;
#	if (/^>(\S+)/) {
#		$name=$1;
#		$prot_len{$name}=0;
#	}
#	else {
#		$prot_len{$name}+=(length($_));
#	}
#}
#close FA;

my @all_files=(@cases,@controls);

foreach my $file (@all_files) {
	my $skip=1;
	open (IN,"<$file") or die "can not open file $file";
	while (<IN>) {
		chomp;
		if (/^Locus/ || /^Unique/) {
			$skip=0;
			next;
		}
		next if ($skip==1);
		last if (/^Unfiltered/);
		my @line=split(/\t/,$_);
		if ($line[0] eq '' || $line[0]=~/^\*/) {
			next;
		}
		if ($line[0]=~/^\d+$/ && $line[1]!~/^\d+$/) {
			next;
		} 
		else {
#			print "$line[0]\n";
			if (!defined ($index_id{$line[0]})) {
				$prot_len{$line[0]}=$line[4];
				$index_id{$line[0]}=$index_count;
				$index_count++;
				$index_annot{$line[0]}=$line[8];
			}
			$raw_value{$file}{$line[0]}=$line[2];
		}
		
	}
	close IN;
}

open (IND,">./index.txt") or die "can not create output file at current dir";
open (MAT,">./sparseMatrix.txt") or die "can not create output file at current dir";

foreach my $protein (sort {$index_id{$a} <=> $index_id{$b}} keys %index_id) {
	print IND "$index_id{$protein}\t$protein\t$index_annot{$protein}\n";
}
close IND;

#print MAT "#ClassDescription	1\n#ClassDescription	-1\n";
foreach my $case (@cases) {
	my %normalized_value;
	my $norm_sum=0;
	print MAT "#".abs_path($case)."\n";
	print MAT "1";
	foreach my $p (keys %{$raw_value{$case}}) {
		if (!defined ($prot_len{$p})) {
			die "undefined length for protein:$p";
		}
		$normalized_value{$p}=$raw_value{$case}{$p}/$prot_len{$p};
		$norm_sum+=$raw_value{$case}{$p}/$prot_len{$p};
	}
	foreach my $prot (sort {$index_id{$a} <=> $index_id{$b}} keys %{$raw_value{$case}}) {
		if (defined ($normalize)) {
			print MAT " $index_id{$prot}:".($normalized_value{$prot}/$norm_sum);
		}
		else {
			print MAT " $index_id{$prot}:$raw_value{$case}{$prot}";
		}
	}
	print MAT "\n";
}

foreach my $control (@controls) {
	my %normalized_value;
        my $norm_sum=0;
        print MAT "#".abs_path($control)."\n";
        print MAT "-1";
        foreach my $p (keys %{$raw_value{$control}}) {
                $normalized_value{$p}=$raw_value{$control}{$p}/$prot_len{$p};
                $norm_sum+=$raw_value{$control}{$p}/$prot_len{$p};
        }
        foreach my $prot (sort {$index_id{$a} <=> $index_id{$b}} keys %{$raw_value{$control}}) {
                if (defined ($normalize)) {
                        print MAT " $index_id{$prot}:".($normalized_value{$prot}/$norm_sum);
                }
                else {
                        print MAT " $index_id{$prot}:$raw_value{$control}{$prot}";
                }
        }
        print MAT "\n";
}





