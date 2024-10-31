


no_predecessor(R) :- reactant(R,_) ; not product(R,_).

path(G1,G2) :- reactant(G1,R), product(G2,R), G1!=G2.
path(G1,G3) :- path(G1,G2), path(G2,G3), G1!=G2, G1!=G3, G2!=G3.

ancestor(G1,G2) :- path(G1,G2), no_predecessor(G1), intermediate(G2), G1!=G2.

% #show path/2.
#show ancestor/2.