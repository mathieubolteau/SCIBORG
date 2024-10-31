



no_successor(R) :- product(R,_) ; not reactant(R,_).

nb_no_successors(N) :- N={no_successor(_)}.

#show no_successor/1.
#show nb_no_successors/1.

