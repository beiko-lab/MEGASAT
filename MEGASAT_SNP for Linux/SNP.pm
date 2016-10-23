###class SNP
package SNP;

use strict;
use warnings;

sub new{
	my ($class, $type, $col, $lratio, $rank_1, $rank_2, $rank_3, $rank_4) = @_;
	my $self = {
		type => $type,
		col => $col,
		lratio => $lratio,
		rank_1 => $rank_1,
		rank_2 => $rank_2,
		rank_3 => $rank_3,
		rank_4 => $rank_4,
	};
	bless ($self, $class);
	return $self;
}
return 1;