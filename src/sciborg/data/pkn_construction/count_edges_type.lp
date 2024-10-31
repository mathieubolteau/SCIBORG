

nb("ACTIVATION", N) :- N = #count{ R: reaction(R, "ACTIVATION")}.
nb("INHIBITION", N) :- N = #count{ R: reaction(R, "INHIBITION")}.
nb("PART_OF", N) :- N = #count{ R: reaction(R, "PART_OF")}.

#show nb/2.