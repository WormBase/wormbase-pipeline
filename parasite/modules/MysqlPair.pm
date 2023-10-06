use strict;
use Carp;
use ProductionMysql;
package MysqlPair;
sub new { 
  my ($class, $from, $to) = @_;
  Carp::croak "Usage: MysqlDump->new(\$from, \$to)" unless $from and $to;
  $from = ProductionMysql->new($from) unless ref $from eq 'ProductionMysql';
  $to = ProductionMysql->new($to) unless ref $to eq 'ProductionMysql';
  return bless {from => $from, to => $to}, $class;
}
sub dump {
  my ($self, $from_pattern, $to_pattern, $table) = @_;
  my $core_db_from = $self->{from}->core_db($from_pattern);
  my $core_db_to = $self->{to}->core_db($to_pattern);
  my $from_cmd = $self->{from}->{db_cmd} . " mysqldump $core_db_from $table " ;
  my $to_cmd =  $self->{to}->{db_cmd} . " $core_db_to " ;
  return system ("{ $from_cmd | $to_cmd 2>&1 1>&3 | grep -v insecure 1>&2; } 3>&1 ");
}
1;
