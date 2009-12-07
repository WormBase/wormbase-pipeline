#!/software/bin/perl -w

use LWP::Simple;

my $user;
my $spamlist;
while($user = shift){

    if (&check($user) ){
	$spamlist .= "SPAM $user\n";
	system("echo $user SPAM >> ~/forum_users");
    }
    else {
    
	open( OUTLOG, "|mailx -s \"WormBase Forum\" $user " );
	print OUTLOG "Dear user,
Thank you for registering with the WormBase Forum.  As your email address is not from an academic or other easily identifiable institution I can not automatically approve your account.  I would appreciate it if you could reply to confirm that you are a genuine user so that I can approve your account as soon as possible.  I hope you understand that we take these measures to keep the forum as free from spam and irrelevant material as possible.\nIf you do not respond within 2 weeks I will assume the account is no longer needed and remove it.\n\nBest wishes\nAnthony\n";
	close OUTLOG or die "didn't close mail properly\n\n";
	
	system("echo $user >> ~/forum_users");
    }
}
print $spamlist;
exit;

sub check {
    my $mail = shift;
    my $file = `wget -O /tmp/$mail "http://www.stopforumspam.com/api?email=$mail"`;
    open (SPAM,"</tmp/$mail");
	while(<SPAM>) {
	    if( /appears\>yes/ ) {
		close SPAM;
		`rm /tmp/$mail`;
		return 1;
	    }
	}
    close SPAM;
    `rm /tmp/$mail`;
}
