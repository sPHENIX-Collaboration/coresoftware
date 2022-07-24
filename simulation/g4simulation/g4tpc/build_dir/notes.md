1. From the test files that I have, it appears that most tracks have hits in most of the gem pad rows. 
-> I assume that it is worth the tradoff to save the energy centroids as array< array<float, 5>, 54> for the data
   per track, as opposed to having a map or vector to hold the 90% the array<float, 5>  that exists.

2. vector< array< array<float,5>, 55>> energy centroids;
