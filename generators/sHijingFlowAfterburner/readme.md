# Flow Afterburner for sHIJING
## Motivation
The output of the sHIJING generator is a HepMC-formatted text file with final state particles, as well as other event-level information. As written, the particles do not possess any global angular flow-like correlation. We have written an afterburner for sHIJING, which parses the HepMC output, and writes another file in the same format, but where the azimuthal angle of individual particles is rotated by a *small* amount to endow them with a $p_T$- and centrality-dependent $v_2$. In this manner, local jet correlations are also preserved.

## Mathematical Statement of the Problem
Let the azimuthal angle $\phi$ of particles at a given $p_T$ be a random variable distributed uniformly in $[0, 2\pi]$ according to 
$$f(\phi) = \frac{1}{2\pi}.$$

We want to derive a mapping $y:\phi\rightarrow\phi'$, where $\phi'$ is a new random variable distributed according to
$$g(\phi')=\frac{1}{2\pi}[1+2v_2\cos(2\phi')].$$

It can be shown that the distributions $f(\phi)$, $g(\phi')$, and the mapping $\phi(\phi')$ must satisfy the following differential equation
$$g(\phi') = f[\phi(\phi')]\frac{d\phi}{d\phi'}.$$
Thus, we would like to solve the equation above for $\phi(\phi')$. Plugging in the expressions for $f$ and $g$ results in, fortunately, a separable equation which can be solved, yielding the following transcendental equation
$$\phi'+v_2\sin(2\phi')-\phi-\frac{1}{2}=0.$$
Given the desired value of $v_2$, this equation must be solved numerically for $\phi'$.

##Parameterization of Centrality-Dependent $v_2(p_t)$
The above transcendental equation is solved numerically *for every particle* in the event record, given the particle's desired $v_2$ and the event centrality, using the secant root-finding algorithm. For a first implementation, we take centrality-dependent measurements of $v_2(p_T)$ for identified pions in Au+Au collisions at $\sqrt{s_{NN}}=200$ GeV, published by PHENIX in PRC 88, 064910. The measurements in three centrality classes, 0-20%, 20-40%, and 40-60%, are fit with a polynomial within $0 < p_T < 5$ GeV, and a linear extrapolation of the fit is used above that, matching the derivative at $p_T=5$ GeV. 

The event centrality is determined based on the impact parameter value, as output by sHIJING in the HepMC file. Instead of examining event multiplicity at forward rapidity and correlating that with impact parameter, we plotted the distribution of impact parameters in minimum bias sHIJING Au+Au events at $\sqrt{s_{NN}}=200$ GeV, and sliced it into quantiles to categorize centrality. This determines the values of the impact parameter that are used in the afterburner to categorize centrality. Events more peripheral than 60% are treated as belonging to the 40-60% category for the purposes of adding flow modulation.

## The Code
The afterburner code can be found on the sPHENIX repository, under **coresoftware/generators/sHijingFlowAfterburner/hepMCFlowAfterburner.C**. It takes two arguments, corresponding to the input HepMC file name, and the desired output file name. The code is thoroughly commented. As of July 2018, the modulation only comprises an elliptic component. Higher order flow can be added by solving the corresponding differential equation. The flow modulation, as currently written, is applied to all final-state particles. A more realistic treatment would be to apply it to parent particles and not decay products, thus better preserving jet structure. This is left for future work.
