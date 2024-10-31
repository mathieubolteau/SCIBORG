% Search of the maximal set of pairs of cells where the expression 
% of the k genes selected are the same for the two cells of the pair.
%
% INPUTS:
%   - pert(C,G,S,CL): expression S of the gene G for the cell C belonging to the class CL 
%   - input(G): gene G is an input, i.e. a node without any predecessor in the PKN
% OUTPUTS:  
%   - selgene(G): gene G is a selected gene belonging to the k genes
%   - affinity(I,J): a pair of cells where their expression for the k selected genes 
%     are the same for both cells. I and J belong to two different classes.


% Generate combinations of k genes
{selinput(G) : pert(C,G,S,CL), not intermediate(G)} = i .
{selinter(G): intermediate(G)} = k-i.

1{selinput(G)  : ancestor(G,I), input(G)} :- selinter(I).
% 1{selinter(I)  : ancestor(G,I), intermediate(I)} :- selgene(G).

% Generate the corresponding perturbation vectors
selpert(E,V,S,C) :- selinput(V), pert(E,V,S,C).
selpert(E,V,S,C) :- selinter(V), pert(E,V,S,C).

% Generate a equal/3 predicate for cells I and J, from two different classes, where their expression for the gene G is equal
equal(I,J,G) :- selpert(I,G,S1,C1), selpert(J,G,S2,C2), I!=J, S1==S2.


pot_affi(I,J) :- k={equal(I,J,_)},  selpert(I,_,_,C1), selpert(J,_,_,C2), C1<C2, I!=J.

% Generate, or not, an affinity/2 predicate when the number of countequal/3
% is equal to k for cells I and J of different classes, which guarantees that all genes are similarly expressed in both cells I and J
0{affinity(I,J)}1 :- pot_affi(I,J).

% % Count the number of input genes expressed at 1 for each affinity/2
nbInputOnes(C,N) :- N={pert(C,G,1,_) : selinput(G), input(G)}, affinity(C,_).

% % Forbid the affinity/2 where the number of input genes expressed at 1 is less than 1
:- affinity(C,_), nbInputOnes(C,N), N < 1.

:- affinity(I,J1), affinity(I,J2), J1!=J2.
:- affinity(I1,J), affinity(I2,J), I1!=I2.

:- pot_affi(I1,J), pot_affi(I2,J), affinity(I2,J), I1<I2.
:- pot_affi(I,J1), pot_affi(I,J2), affinity(I,J2), J1<J2.


selgene(G) :- selinput(G).
selgene(G) :- selinter(G).

% Maximize the number of pairs of cells
#maximize{1, I : affinity(I,_)}.

% Show selgene/1 and affinity/2 predicates
#show selgene/1.
#show affinity/2.

#show selinput/1.
#show selinter/1.
