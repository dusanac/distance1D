#!/usr/bin/perl

use Tk;

#GLOB
#-------------------------
$Tolerance = 0.00; 
$ar1,$ar2 = 0;
$ab1 = 0;
$ab2 = 0;
$path1, $path2;

# B O D O V I
#---------------------------------------------------------
$NASbod=10;  	# not aligned spots
$BCSbod=10;		# broj cjepanih spotova (ostavljena opcija)
$ABDbod=2;		# abundance difference
$MWDbod=1;		# MW difference

# R E Z U L T A T
#---------------------------------------------------------
$DIST=0;
$MWD=0;
$ABD=0;
$TOLER=0;
$DISTswp=0;
$MWDswp=0;
$ABDswp=0;
$TOLERswp=0;

sub uparivanje {#args: $ref1 $ref2
	my @mwI = ();
	my @mwII = ();
	local @upareni = (); 	#mora local jer je koristi funkcija dif_par
	local @upareniI = ();	#mora local jer je koristi funkcija dif_par
	local @upareniII = ();	#mora local jer je koristi funkcija dif_par
	my @a1 = ();
	my @a2 = ();
	my @repet1 = ();
	my @repet2 = ();
	my @tmp = ();
	
	my $tmp, $a1,$a2,$p1,$p2, $razlika1,$razlika2 = 0;
	$tolerance = $Tolerance;	#KDa Vidjeti za ubuduce da li izraziti u procentima jer vece mol mase imaju
					#vece odstupanje
	@mwI = @{$_[0]};    #referencni parametri
	@mwII = @{$_[1]};
	
	#print "MW I  @mwI\n";
	#print "MW II @mwII\n";
	for (my $i=0; $i<@mwI;$i++) {
		my $upario = 0;
		
		for (my $j=0;$j<@mwII;$j++) {
			if ($mwI[$i] >= ($mwII[$j]-($tolerance*$mwI[$i])) && $mwI[$i] <= ($mwII[$j]+($tolerance*$mwI[$i])) ){#moze i obratno
				push (@upareni , (($i+1) . ':' . ($j+1)) );
				$upario = 1;
				push (@upareniI , $mwI[$i]);
				push (@upareniII, $mwII[$j]);
			}
		
		}
		unless ($upario) {push (@upareni, ($i+1) . ':' . "0");}
	}
	for (my $i=0; $i<@mwII;$i++){
		my $nasao = 0;
		for (my $j=0;$j<@upareniII;$j++){
			if ($mwII[$i]==$upareniII[$j]){$nasao=1;last;}
		}
		if (! $nasao){push (@upareni, ('0:' . ($i+1)));}
	}	 
	
	if (!test_tolerance (\@upareniI,\@upareniII) ){
		
		#print "Prevelika tolerancija pri uparivanju!\n";
		$total = 1000;
		$mwd = 100;
		$abd = 100;
		#exit 1;
	}
	

# Algoritam za poboljsanje uparivanja po sistemu stednje
# 
#------------------------------
	foreach (@upareni) {
		($a1, $a2) = split (":",$_,2);
		chomp ($a2);#remove the \n
		push (@a1,$a1);
		push (@a2,$a2);
	}
	
	@repet1 = occurence_in_arr_place (\@a1);
	@repet2 = occurence_in_arr_place (\@a2);
	
	for (my $i=0; $i<@repet1;$i+=2){
		my($p1,$p2) = split (":",@repet1[$i],2);
		if ( dif_par ($p1) < dif_par($p2) ){
				@upareni[$p2] = "0:@a2[$p2]";
				@a1[$p2] = 0;
		}
		else {
				@upareni[$p1] = "0:@a2[$p1]";
				@a1[$p1] = 0;
		}
	}
	for (my $i=0; $i<@repet2;$i+=2){
		my($p1,$p2) = split (":",@repet2[$i],2);
		
		if ( dif_par ($p1) < dif_par($p2) ){
				@upareni[$p2] = "@a1[$p2]:0";
				
		}
		else {
				@upareni[$p1] = "@a1[$p1]:0";
		}
	}
	
	sub dif_par {
		my $b1,$b2, $poz, $dif = 0;
		$poz = $_[0];
		($b1, $b2) = split (":",$upareni[$poz],2);	
		if ($b1 == 0 || $b2 ==0){return 100;}#beskonacno velika razlika
		$dif = abs ($upareniI[$b-1] - $upareniII[$b2 - 1]);
		return $dif;
	}
#	print "Inside uparivanje(): @upareni\n";
	for (my $i =0;$i<@upareni;$i++){
		push (@tmp, $upareni[$i]) unless (@upareni[$i] eq "0:0");
	}
	
	 
	@upareni = @tmp;
	
#	print "Inside uparivanje(): @upareni\n";
#	print "Inside uparivanje(): @upareniI\n";
#	print "Inside uparivanje(): @upareniII\n";

	@upareniI = del_the_same (\@upareniI);
	@upareniII = del_the_same (\@upareniII);
#-----------------------------------------------------	
	#print "Inside uparivanje() upareni : @upareni\n";
	#print "Inside uparivanje() upareni1: @upareniI\n";
	#print "Inside uparivanje() upareni2: @upareniII\n";
	
 	return (\@upareni,\@upareniI,\@upareniII);
}

sub del_the_same{#param ref of @array
				
	my @rep = ();
	my @arr = ();
	my @del = ();#mora posebno u red i 0 radi greske u interpretaciji
	my @copy = ();
	my $a1,$a2,$rep =0;
	@arr = @{$_[0]};#
	@rep = occurence_in_arr_place (\@arr);
	for (my $i =0;$i<@rep;$i+=2){#brisanje istih u nizu
		($a1,$a2) = split (":",$rep[$i],2);
		push (@del,$a2);
	}
	foreach (@del){
		@arr[$_] = 0;#pridruzuje 0 umjesto da brise
	}
	for (my $i =0;$i<@arr;$i++){
		unless ($arr[$i] == 0){
			push (@copy,$arr[$i]);
		}
			
	}
	
	return (@copy);
}


sub test_tolerance {#OK param: @ref1,@ref2
	my $upareniIref = $_[0];
	my $upareniIIref = $_[1];
	if (occurence_in_arr ($upareniIref,3) > 0 || occurence_in_arr ($upareniIIref,3) > 0){
		return 0;
	}
	else {return 1;}	#ok 1;  not ok 0
}

sub occurence_in_arr {#OK  	#searchig for the same ocurence in array and return the namber of it
							# function (\@refarray, INT); INT is 3 for occurence of 3 times. 
							# INT defoult is 2
	my $n = 0;
	my @arr = @{$_[0]};
	my $arg = 2;
	my @mjesta = ();
	if (defined $_[1]){$arg = $_[1];}
	
	for (my $i=0; $i<@arr;$i++){
		my $tmp = 0;
		for (my $j=0; $j<@arr;$j++){
			if ($i==$j){next;}#da nebi sebe ubrojao
			else {
				if ($arr[$i] == $arr[$j]){$tmp++;}
			}
		}
		if (($tmp) == ($arg-1)){$n++;}
	}
	return $n/$arg;#da ne racuna za svaki par 2 puta npr: 1>>3 3>>1
 }
 
sub occurence_in_arr_place {#OK #searchig for the same ocurence in array and return the PLACE of it
							
	my $n = 0;
	my @arr = @{$_[0]};
	my $arg = 2;
	my @mjesta = ();
	if (defined $_[1]){$arg = $_[1];}
	
	for (my $i=0; $i<@arr;$i++){
		my $tmp = 0;
		for (my $j=0; $j<@arr;$j++){
			if ($i==$j || $arr[$i] == 0 || $arr[$j]==0){next;}#da nebi sebe ubrojao ili nule
			else {
				if ($arr[$i] == $arr[$j]){
					$tmp++;
					#print "PRONASAO> $i = $j\n"; 
					push (@mjesta, ($i . ":" . $j));}
			}
		}
		
	}
	#print "mjesta @mjesta<\n";
	return (@mjesta);#da ne racuna za svaki par 2 puta npr: 1>>3 3>>1
 }
 
sub not_aligned_spot_br {#OK
	my @arr = @{$_[0]};
	my $b1,$b2 = 0;
	my $rezultat = 0;
	foreach (@arr){
		($b1,$b2) = split (":",$_,2);
		if ($b1 == 0 || $b2 == 0) {$rezultat++;}
	}
	return $rezultat;
}

sub MW_difference {
	local @upareni = (); 
	local @uz1 = (); 
	local @uz2 = (); 
	local $w1, $w2,$zbir = 0; 
	@upareni = @{$_[2]};
	@uz1 = @{$_[0]};
	@uz2 = @{$_[1]};
	foreach (@upareni) {
		($w1,$w2) = split (":",$_,2);
		unless ($w1 == 0 || $w2 == 0){
			$zbir += abs (@uz1[$w1-1]-@uz2[$w2-1]);
		}
	}
	return $zbir;
}

sub Abund_dif {#total params REF 
	my $A1, $A2, $temp, $sum1, $sum2=0;
	my @uz1 = ();
	my @uz2 = ();
	my @diff = ();
	my @diff2 = ();
	my @upareni = @{$_[2]};
	@uz1 = @{$_[0]};
	@uz2 = @{$_[1]};
	$sum1 = sum (@uz1);
	$sum2 = sum (@uz2);
	#print "uz1 @uz1\n";
	#print "uz2 @uz2\n";
	#print "sum1 $sum1\nsum2 $sum2\n";
	for (my $i = 0; $i<@upareni; $i++){
			($A1,$A2) = split (":", $upareni[$i] , 2);
			unless ($A1 == 0 || $A2 == 0){ #iskljucuje neuparene spotove
					push (@diff,  abs ( ($uz1[$A1-1]) - ($uz2[$A2-1]) ) );	
			}
	}
	#print "@diff\n";
	return sum (@diff);
	} 

sub Abund_dif2 {#sub (uzorakA1, uzorakA2, parovi) 
				# uracunava cijepane spotove 
	my $A1, $A2, $temp;
	my @diff = (); 
	my @upareni = @{$_[2]};
	my @uz1 = @{$_[0]};
	my @uz2 = @{$_[1]};
	my @upareniI = ();
	my @upareniII = ();
	
	##print "upareni su: @upareni\n";
	for (my $i = 0 ; $i < @upareni; $i++){#load sve parove u dva niza
		(@upareniI[$i],@upareniII[$i]) = split (":",$upareni[$i],2);
	}
	#print "upareniI:  @upareniI\n";
	#print "upareniII: @upareniII\n";
	#print "uz1-2 prije  : @uz1\t@uz2\n";
	my @occ1 = occurence_in_arr_place (\@upareniI);
	my @occ2 = occurence_in_arr_place (\@upareniII);
	#print "OCC1 @occ1\n";
	if (@occ1){
		#print "pokrece Occ IF\n";
		for (my $i=0; $i<@occ1; $i += 2){ #premjesta rezultat cjepanih spotova
				my ($tmp1,$tmp2) = split (":",$occ1[$i],2);
				unless ($upareniI[$tmp1] == 0 || $upareniI[$tmp2] == 0){
					@uz2[$tmp1] += @uz2[$tmp2];
					@uz2[$tmp2] = 0;
					@upareniI[$tmp2] = 0;
				}
		}
	}
	
	if (@occ2){
		for (my $i=0; $i<@occ2; $i += 2){
				my ($tmp1,$tmp2) = split (":",$occ2[$i],2);
				unless ($upareniII[$tmp1] == 0 || $upareniII[$tmp2] == 0){	
					@uz1[$tmp1] += @uz1[$tmp2];
					@uz1[$tmp2] = 0;
					@upareniII[$tmp2] = 0;
				}
		}
	}
	#print "uz1-2 poslije: @uz1\t@uz2\n";
	#print "upareniI:  @upareniI\n";
	#print "upareniII: @upareniII\n";
			
	
	for (my $i = 0; $i<@upareniI; $i++){#???	
		unless ($upareniI[$i] == 0 || $upareniII [$i] == 0){
			#print "$upareniI[$i] , $upareniII[$i] \n";
			push (@diff, abs (@uz1[$upareniI[$i]-1]-@uz2[$upareniII[$i]-1]));
		}
	}
	return sum (@diff);
}
	
sub sum {# daje zbir svih clanova u nizu;
	my @arr = ();
	my $summ;
	$summ = 0;
	@arr = @_;
	for (my $i = 0;$i<@arr; $i++){
		$summ += @arr[$i];
	}
	return $summ; 
}

sub load_file{
	my $n, $fname;
	
	my @arr1 = ();
	my @arr2 = ();#razdvoio (pravilo probleme)
	
	$fname = @_[0];
	$n=0; #counter
	
	if (open (INFILE, "<$fname")) {
		foreach (<INFILE>){
			if (m/^\#/) {next;}
			chomp;
			($arr1[$n],$arr2[$n]) = split (":", $_, 2);
			$n++;
		}
	return (\@arr1, \@arr2);
	}
	else {
		print "NO File\n";
		return 0;
		}
}

sub scoreOLD{
	my $nas,$bcs;
	my $mwd,$abd,$total;
	my @rvalue = uparivanje ($ar1, $ar2);
	
	$nas = not_aligned_spot_br ($rvalue[0]);
	$bcs = occurence_in_arr ($rvalue[1])+occurence_in_arr ($rvalue[2]);
	$mwd = MW_difference ($ar1,$ar2,@rvalue[0]);
	
	$abd = Abund_dif ($ab1,$ab2,@rvalue[0]);
	print "nas = $nas\nbcs = $bcs\nmwd = $mwd\nabd = $abd\n";
	$total = ($nas*$NASbod)+($bcs*$BCSbod)+($mwd*$MWDbod)+($abd*$ABDbod);
	print "TOTAL DIF: $total\n\n";
}

sub score{
	local $nas,$bcs =0;
	local $mwd,$abd,$total=0;
	
	my @rvalue = uparivanje ($ar1, $ar2);
	
	if ($total == 0){#test tolerance pass, date velike vrijednosti
		$nas = not_aligned_spot_br ($rvalue[0]);
		$bcs = occurence_in_arr ($rvalue[1])+occurence_in_arr ($rvalue[2]);
		$mwd = MW_difference ($ar1,$ar2,@rvalue[0]);
		$abd = Abund_dif ($ab1,$ab2,@rvalue[0]);
		$total = ($nas*$NASbod)+($bcs*$BCSbod)+($mwd*$MWDbod)+($abd*$ABDbod);
		#print "total = $total\n";
	}
	return ($total, $mwd, $abd, $Tolerance);
}

sub swap {
	my $swp =0;
	$swp = $ar1;
	$ar1 = $ar2;
	$ar2 = $swp;
	
	$swp = $ab1;
	$ab1 = $ab2;
	$ab2 = $swp;
}

sub smallest_place {# daje mjasto najmanjeg clana niza
	my @arr = ();
	my $broj, $index = 0;
	@arr = @_;
	if (@arr != 0){#test da li je aray prayna
		$broj = @arr[$i];
		$index = 0;
	}
	for (my $i = 0;$i<@arr; $i++){	
		if (@arr[$i] < $broj){
			$broj = @arr[$i];
			$index = $i;
		}
	}
	return $index; 
}

## MAIN CODE
#-----------------------------------------------------------------------

 sub main {
	 
if ((-e $path1) && (-e $path2)) {
($ar1, $ab1) = load_file ($path1);
($ar2, $ab2) = load_file ($path2);

@distAR = ();
@mwdAR = ();
@abdAR = ();
@tolAR = ();
$tempvar = 0;

for ($Tolerance=0 , $counter = 0; $Tolerance<0.15; $Tolerance += 0.001, $counter ++){
	
	($distAR[$counter], $mwdAR[$counter],$abdAR[$counter],$tolAR[$counter] ) = score ();
	if ($distAR[$counter]==1000){
		$tempvar +=1;
		last if ($tempvar>3);
	}
		
}

$place = smallest_place (@distAR);
$DIST = $distAR[$place];
$MWD = $mwdAR[$place];
$ABD = $abdAR [$place];
$TOLER = $tolAR [$place];
print "DIST = $DIST\tTOLER = $TOLER\n";

swap();
$tempvar = 0;
for ($Tolerance=0 , $counter = 0; $Tolerance<0.15; $Tolerance += 0.001, $counter ++){
	
	($distAR[$counter], $mwdAR[$counter],$abdAR[$counter],$tolAR[$counter] ) = score ();
	if ($distAR[$counter]==1000){
		$tempvar +=1;
		last if ($tempvar>3);
	}
		
}

$place = smallest_place (@distAR);
$DISTswp = $distAR[$place];
$MWDswp = $mwdAR[$place];
$ABDswp = $abdAR [$place];
$TOLERswp = $tolAR [$place];
print "DISTswp = $DISTswp\tTOLERswp = $TOLERswp\n";
if ($DIST>$DISTswp){ #swap / zamjena vrijednosti rezultata
	$tempvar = $DIST;
	$DIST = $DISTswp;
	$DISTswp = $tempvar;
	
	$tempvar = $MWD;
	$MWD = $MWDswp;
	$MWDswp = $tempvar;
	
	$tempvar = $ABD;
	$ABD = $ABDswp;
	$ABDswp = $tempvar;
	
	$tempvar = $TOLER;
	$TOLER = $TOLERswp;
	$TOLERswp = $tempvar;
}
$label3R->configure (-text => $DIST);
$label4R->configure (-text => $MWD);
$label5R->configure (-text => $ABD);
$label6R->configure (-text => $DISTswp);
}
else {
	$label3R->configure (-text => -1);
	$label4R->configure (-text => -1);
	$label5R->configure (-text => -1);
	$label6R->configure (-text => -1);
}
}

#-----------------------------------------------------------------------
#incijacija widget-a
$main = MainWindow->new();

$main->title("distance1D 0.1");

$FrameLeft = $main->Frame (-width => 150, -height => 300);
$FrameRight = $main ->Frame(-width => 150, -height => 300);

#Sledeci kod RIJESAVA PROBLEM promjene geometrije prozora pri racunu

$frameL0 = $FrameLeft->Label(  -text => "                              ")->pack;
$frameR0 = $FrameRight-> Label(-text => "                                              ")->pack;

#Generisanje potrebnih winget-a (polja, dugmadi i sl.)
$label1 = $FrameLeft->Label(-text => "Species 1: ");
$entry1 = $FrameRight->Entry (-textvariable => \$path1);

$label2 = $FrameLeft->Label(-text => "Species 2: ");
$entry2 = $FrameRight->Entry (-textvariable => \$path2);

$label3 = $FrameLeft->Label(-text => "TOTAL DIF:");
$label3R = $FrameRight->Label(-text => $DIST);#defoult value

$label4 = $FrameLeft->Label(-text => "(mwd):");
$label4R = $FrameRight->Label(-text => $MWD);

$label5 = $FrameLeft->Label(-text => "(abd):");
$label5R = $FrameRight->Label(-text => $ABD);

$label6 = $FrameLeft->Label(-text => "swDIST");
$label6R = $FrameRight->Label(-text => $DISTswp);

$buttonR = $main->Button(-text => 'Submit', -command => sub {main();});

#Geometrijske funkcije elemenata prozora
$FrameLeft -> pack (-side => 'left');
$FrameRight -> pack (-side => 'left');

$label1-> pack (-side => 'top');
$entry1-> pack (-side => 'top');

$label2-> pack (-side => 'top');
$entry2-> pack (-side => 'top');

$label3-> pack (-side => 'top');
$label3R-> pack (-side => 'top');

$label4-> pack (-side => 'top');
$label4R-> pack (-side => 'top');

$label5-> pack (-side => 'top');
$label5R-> pack (-side => 'top');

$label6-> pack (-side => 'top');
$label6R-> pack (-side => 'top');

$buttonR-> pack (-side => 'bottom');

$label3R->configure (-fg => '#A11818');


#POKRETANJE PROZORA
MainLoop();

