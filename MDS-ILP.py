#! /usr/bin/env python3

import numpy as np
import gurobipy as gp
import decycling_utilities as utilities
import decycling_constraints as constraints

import igraph
import networkx as nx
import time

from gurobipy import GRB
from typing import Iterable, Union
from argparse import ArgumentParser


oracle_time = 0.0
cycles_time = 0.0
longest_path_time = 0.0

def build_decycling_oracle(b: int, k: int, l: int, DS: gp.MVar):
    V = set(range(b ** k))
    G = nx.DiGraph(
        igraph.Graph.De_Bruijn(b, k).to_networkx()
    )

    def oracle(solver: gp.Model, where):
        global oracle_time, cycles_time, longest_path_time
        oracle_start = time.time()
        if where == GRB.Callback.MIPSOL:
            X = solver.cbGetSolution(DS).nonzero()[0]
            H = nx.induced_subgraph(G, V - set(X))

            decycled = True
            cycles_start = time.time()
            cycles = nx.cycles.simple_cycles(H)
            for _, cycle in zip(range(1, 4), cycles):
                decycled = False
                solver.cbLazy(1 <= gp.quicksum(
                        DS[x].item() for x in cycle
                ))
            cycles_time += time.time() - cycles_start

            if decycled:
                longest_path_start = time.time()
                path = nx.dag_longest_path(H)

                if (L := len(path)) > l:
                    solver.cbLazy(
                        L // l <= gp.quicksum(
                            DS[x].item() for x in path
                        )
                    )
                longest_path_time += time.time() - longest_path_start

        oracle_time += time.time() - oracle_start
    return oracle


def DS_ILP_setup(
        solver: gp.Model, b: int, k: int, length: int,
        minimize_expectation: bool, partition: Union[None, Iterable[int]]
):
    solver.ModelSense = GRB.MINIMIZE
    In, Out = utilities.DBG_walk_maps(b, k)

    # DAG and path length enforcement
    DS, C, P = constraints.walk_breaking(
        b ** k, length, In, Out, solver, with_probs=minimize_expectation
    )

    # Set a higher length_bound to add more cycle constraints.
    # The cycle-breaking and PCR constraints are redundant but make optimization ~10-100x faster.
    cycles = constraints.DeBruijn_cycles(
        b, k, length_bound=k+2, X=DS
    )

    for i, cycle in enumerate(cycles):
        solver.addConstr(
            cycle, f'cycle-{i}'
        )

    # Add strong PCR and PNR constraints when no partition given
    if partition is None:
        for i, cycle in enumerate(utilities.PCR(b, k)):
            solver.addConstr(
                DS[cycle].sum() == 1, f'PCR-{i}'
            )

        for t, c in enumerate(utilities.PCCR(b, k)):

            # Special constraints for base MDS problem.
            #   Gives the solver freedom to use an ordered-set branch-and-bound heuristic.
            if (nc := len(c)) > 2:
                for I in range(nc):
                    J, K = (I + 1) % nc, (I + 2) % nc
                    solver.addSOS(
                        GRB.SOS_TYPE2, DS[[I, J, K]].tolist(), [1, 1, 1]
                    )

            # Optimize the circuit instead of the size
            #   since we know the MDS size in the base graph.
            solver.setObjective(C[-1].sum(), sense=GRB.MINIMIZE)

            return DS, C, P, None

    # Minimize the selected partitions, while enforcing the max-path constraint.
    solver.addConstr(C[-1].sum() == 0, name='nilpotent')

    # Build the partition classes and linearize the indexing
    partition_ids, linear_ids = np.unique(partition, return_inverse=True)

    p = solver.addMVar(
        shape=(partition_ids.size,), vtype=GRB.BINARY, name="partition-selector"
    )

    solver.addConstr(DS == p[linear_ids], name="partition-selection")

    # Minimize selected classes
    solver.setObjectiveN(p.sum(), index=0, priority=1)

    # Experimental secondary objectives: minimize/maximize selected vertices
    # solver.setObjectiveN(DS.sum(), index=1, priority=0)  # maximize
    # solver.setObjectiveN(-DS.sum(), index=1, priority=0)  # minimize

    for i, cycle in enumerate(utilities.PCR(b, k)):
        solver.addConstr(
            DS[cycle].sum() >= 1, f'cycle-{i}'
        )

    # Secondary objective: total decycling set size
    # solver.setObjectiveN(DS.sum(), index=1)

    # WARNING: Needs an update to accommodate for partitions.
    #   In particular, we need to slightly alter how we calculate expected path length.
    #   The difference is minor and quite subtle, yet it changes the optimal value by a large margin.

    # Set hierarchical objective, i.e., minimize size then expected length.
    # This is scaled up by a constant, so it will appear larger than it actually is.
    #
    # if minimize_expectation:
    #     solver.setObjectiveN(
    #         (np.arange(P.shape[0]) @ P).sum() / i, index=1
    #     )

    return DS, C, P, p


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('b', type=int, help='set base b alphabet')
    parser.add_argument('k', type=int, help='set order k De Bruijn graph')
    parser.add_argument('length', type=int, help='set max remaining path length')
    parser.add_argument('--sym', action='store_true', help='consider x ~ rev(not(x))')
    parser.add_argument('--out', help='write output file to given path', default=None)
    parser.add_argument('--verbose', action='store_true', help='allow solver messages when True')
    # Needs a fix before use.
    # parser.add_argument('--EPL', action='store_true', help='additionally minimize expected path length')
    args = parser.parse_args()

    b, k, l = args.b, args.k, args.length - 1

    kmer = utilities.kmer_formatter(b, k)
    In, Out = utilities.DBG_walk_maps(b, k)

    solver = gp.Model()

    if not args.verbose:
        solver.setParam(GRB.Param.LogToConsole, False)

    # solver.setParam(GRB.Param.Cuts, 3)
    solver.setParam(GRB.Param.Presolve, 2)
    solver.setParam(GRB.Param.PreCrush, 1)
    solver.setParam(GRB.Param.Symmetry, 2)
    solver.setParam(GRB.Param.PrePasses, 5)
    solver.setParam(GRB.Param.LazyConstraints, 1)

    print(
        f'\nSetting up {l+1}-walk circuit for the base {b} De Bruijn graph of order {k}…', end='\n\n'
    )

    knonical = None
    if args.sym:
        cp = utilities.Complementer(b, k)
        # Select whichever symmetries you want
        knonical = np.asarray([
            min(x, cp.rc(x)) for x in range(cp.d)
        ])

    DS, C, P, p = DS_ILP_setup(
        solver, b, k, l, False, knonical
    )

    # Do not use until the chordless cycles are fixed…
    oracle = build_decycling_oracle(b, k, l, DS)
    solver.optimize(oracle)

    print(f'\nOptimizing…')
    solver.optimize()

    MDS = DS.X.nonzero()[0]
    # parts seems unused
    # parts = None
    # if p is not None:
    #     parts = p.X.nonzero()[0]

    print(
        '\nOptimal MDS:',
        f'\tVertices selected: {DS.x.sum()}',
        sep='\n'
    )

    if p is not None:
        print(
            f'\tPartitions selected: {p.x.sum()}',
        )

    print(f'oracle {oracle_time} cycles {cycles_time} longest path {longest_path_time}')

    if args.out is not None:
        with open(args.out, 'w') as out:
            out.writelines(
                '\n'.join(kmer(x) for x in MDS)
            )
