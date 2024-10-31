



complexes(P) :- product(P,R) ; reaction(R,"PART_OF").

nb_complexes(N) :- N={complexes(_)}.

#show complexes/1.
#show nb_complexes/1.

