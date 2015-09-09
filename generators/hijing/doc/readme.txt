

	Version 2.03

	ESK parameterization of nuclear shadowing of parton distributions
	is replaced by a user implemented subroutine SHADOW. Users can replace
	this subroutine with their own parametrization

     	Versio 2.02

     	in order to fit the total cross section and dn/dy in pp or ppbar
     	collisions the pt cut-off now varies with energy roughly as
     	p0=1.3+0.09*log(sqrt(s)/10)+0.065*log(sqrt(s)/10)**2


     	Version 2.01

     	Major modification started on Feb. 22, 1999:
     	Include options for different parametrization of parton distributions.
     	IHPR2(7):(D=4)
               =1 or 2 Duke-Owns set 1 or 2
               =3 MRS D-' parametrization
               =4 GRV HO
               =5 CTEQ 2pM


     	IHPR2(6):(D=1)
               =1 the original HIJING parametrization of shadowing
               =2 Eskola-Kolhinen-Salgado parametrization

	*****************
	*		*
     	* Version 2.00	*
	*		*
	*****************

     	Modification Oct. 8, 1998. In hijing, log(ran(nseed)) occasionally
     	causes overfloat. It is modified to log(max(ran(nseed),1.0e-20)).

     	replaces all the generic random number generator ran(nseed) by the 
	portable random number generator in jetset  RLU(0).
 
	add the following common block to store the formation time and position
	of each hadron
	COMMON/HIMAIN3/VATT(130000,4)
	For directly produced or massless particles (gamma, dilepton, etc),
	we set VATT(I,4)=-10, and VATT(I,1-3)=0
	HIPR1(50)=particle formation time in rest-frame with default value 
	1.0 fm



     	Version 1.36


	(5/30/98) Some statements in HIJSFT ae modified to avoid questionable
     	branch into another IF block (new compiler on alpha machines does not
     	like this).

	Nothing important has been changed here. A few 'garbage' has been
	cleaned up here, like common block HIJJET3 for the sea quark strings
	which were originally created to implement the DPM scheme which
	later was abadoned in the final version. The lines which operate
	on these data are also deleted in the program.



	*************************
	*			*
	*     Version 1.35	*
	*			*
	*************************

	Last modification on Feb. 24, 1997.
	There are some changes in the program: subroutine HARDJET is now
	consolidated with HIJHRD. HARDJET is used to re-initiate PYTHIA
	for the triggered hard processes. Now that is done  altogether
	with other normal hard processes in modified JETINI. In the new
	version one calls JETINI every time one calls HIJHRD. In the new
	version the effect of the isospin of the nucleon on hard processes,
	especially direct photons is correctly considered.
	For A+A collisions, one has to initilize pythia
	separately for each type of collisions, pp, pn,np and nn,
	or hp and hn for hA collisions. In JETINI we use the following
	catalogue for different types of collisions:
	h+h: h+h (I_TYPE=1)
	h+A: h+p (I_TYPE=1), h+n (I_TYPE=2)
	A+h: p+h (I_TYPE=1), n+h (I_TYPE=2)
	A+A: p+p (I_TYPE=1), p+n (I_TYPE=2), n+p (I_TYPE=3), n+n (I_TYPE=4)


	****************************
	*			   *
	*       Version 1.34	   *
	*			   *
	****************************

       Last modification on January 5, 1998. Two misstakes are corrected in
       function G. A Misstake in the subroutine Parton is also corrected.
       (These are pointed out by Yasushi Nara).


       Last modifcation on April 10, 1996. To conduct final
       state radiation, PYTHIA reorganize the two scattered
       partons and their final momenta will be a little
       different. The summed total momenta of the partons
       from the final state radiation are stored in HINT1(26-29)
       and HINT1(36-39) which are little different from 
       HINT1(21-24) and HINT1(41-44).





	*************************
	*			*
	*      Version 1.33	*
	*			*
	*************************

       Last modfication  on September 11, 1995. When HIJING and
       PYTHIA are initialized, the shadowing is evaluated at
       b=0 which is the maximum. This will cause overestimate
       of shadowing for peripheral interactions. To correct this
       problem, shadowing is set to zero when initializing. Then
       use these maximum  cross section without shadowing as a
       normalization of the Monte Carlo. This however increase
       the computing time. IHNT2(16) is used to indicate whether
       the sturcture function is called for (IHNT2(16)=1) initialization
       or for (IHNT2(16)=0)normal collisions simulation

       Last modification on Aagust 28, 1994. Two bugs associate
       with the impact parameter dependence of the shadowing is
       corrected.

	*********************************
	*	HIJING version 1.31	*
	*********************************

	No major change has been made in this version. A bug 
associated with triggering direct photon is corrected. This 
problem was reported by Jim Carroll and Mike Beddo. Another 
bug in tracing decay trees reported by Matt Bloomer was also 
corrected.

	Another feature of triggered jet production (including 
direct photon) in HIJING needs to be clarified here. The energy 
and momentum stored in HINT1(21)-HINT1(25); HINT1(31)-HINT1(35); 
HINT1(41)-HINT1(45); and HINT1(51)-HINT1(55) are the momenta 
of the hard scattered partons in the lowest order calculation. 
To include the higher order correction, final state radiations 
are implemented. In this implementation, momentum and energy 
are reshuffled so that the momentum and energy you find in the 
finally produced parton or particle might not match those in 
HINT1(21)-HINT1(55). But the momentum is still balanced for 
the back-to-back jets (or jet and photon). So information in 
HINT1(21)-HINT1(55) still give the right direction of the two 
high pt jets.

	*********************************
	*	HIJING Version 1.30	*
	*********************************

	The option to trig on heavy quark production (charm IHPR2(18)=0 
or beauty IHPR2(18)=1) is added. To do this, set IHPR2(3)=3. For inclusive 
production, one should reset HIPR1(10)=0.0. One can also trig larger pt
QQbar production by giving HIPR1(10) a nonvanishing value. The mass 
of the heavy quark in the calculation of the cross section (HINT1(59)-
-HINT1(65)) is given by HIPR1(7) (the default is the charm mass D=1.5). 
We also include a separate K-factor for heavy quark and direct photon 
production by HIPR1(23)(D=2.0).

	*********************************
	*	HIJING version 1.20	*
	*********************************

Beside some modifications to correct some minor errors, there three
major changes since the first version.

	(1) An option is added to keep the information of all 
particles including the decay chains. The switch of this option 
is IHPR(21)=1 (default=0). The array KATT(130000,2) in 
COMMON/HIMAIN2/ is expanded to KATT(130000,4). K(I,3) is the line 
number of the parent particle of the current line which is produced 
via a decay. KATT(I,4) is the status number of the particle: 11=particle
which has decayed; 1=finally or directly produced particle.

	(2)Heavy flavor production processes through q+\bar{q}-->c+\bar{c}
and g+g-->c+\bar{c} have been switched on. Charm production is the default, 
B-quark option is IHPR2(18)=1 (default=0). Note that when B-quark production
is switched on, charm quark is automatically off.

	(3) A subroutine  PROFILE(XB) is added to assist the simulation 
of minimum biased events of proton-nucleus and nucleus-nucleus collisions. 
Refer to the examples in hijing.tex for the simulation.



