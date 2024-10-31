

node(N) :- reactant(N,_).
node(N) :- product(N,_).

nb_nodes(N) :- N={node(_)}.

#show node/1.
#show nb_nodes/1.