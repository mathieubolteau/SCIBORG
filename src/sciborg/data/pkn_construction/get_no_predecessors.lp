



no_predecessor(R) :- reactant(R,_) ; not product(R,_).

nb_no_predecessors(N) :- N={no_predecessor(_)}.

#show no_predecessor/1.
#show nb_no_predecessors/1.

