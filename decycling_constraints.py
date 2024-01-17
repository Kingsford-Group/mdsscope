import igraph
import networkx as nx
import gurobipy as gp

from tqdm import trange
from gurobipy import GRB
from typing import Callable

# Walk operators return a node neighborhood.
# Must return a list rather a general iterable because of Gurobi's indexing API.
Walk = Callable[[int], list[int]]


def intersection_lb(x, y, count=1):
    return x[y].sum() >= count


def intersection_ub(x, y, count=1):
    return x[y].sum() <= count


def cycles(graph: nx.Graph, X: gp.MVar, length_bound: int):
    """
        Return an iterator of intersection constraints of X and G's chordless
            cycles, returning only cycles with length <= length_bound.
    """

    cycles = nx.cycles.simple_cycles(
        graph, length_bound=length_bound
    )

    for cycle in cycles:
        yield intersection_lb(X, cycle)


def DeBruijn_cycles(b: int, k: int, length_bound: int, X: gp.MVar):
    """
        Creates an intersection constraint iterator over the De Bruijn graph of base b and order k.
            These enforce that X and the chordless cycles of length <= length_bound overlap at least once.
    """

    G = nx.DiGraph(
        igraph.Graph.De_Bruijn(b, k).to_networkx()
    )

    cycles = nx.cycles.simple_cycles(
        G, length_bound=length_bound
    )

    for cycle in cycles:
        yield intersection_lb(X, cycle)


def walk_breaking(
        n: int, length_bound: int, In: Walk, Out: Walk, model: gp.Model, with_probs: bool = False
):
    """
        Builds the constraint and variables to break all cycles and paths bounded by length_bound.

        When with_probs is true, we also build a random walk circuit and random walk probability vector.

        Returns:
            1. (1, n) binary blocked vertex variables
            2. (length_bound + 1, 1) active circuit length variables
            3. (length_bound + 1, n) probability circuit variables if with_probs else None
    """

    #
    # The ILP uses a decycling set indicator vector and a path circuit C.
    #
    #      Name       |       Kind          Variables                Constraints
    #  decycling set  |      binary             bᵏ                       ~bᵏ
    #  path circuit   |    continuous        ~l x bᵏ                   ~l x bᵏ
    #

    # The walk circuit encodes Boolean matrix multiplication that we use
    #   to enforce that all paths of length <= l and all cycles are broken.
    # blocked = model.addMVar(
    #     shape=(n,), vtype=GRB.BINARY
    # )

    # boolean_circuit = ForwardCircuit(
    #     blocked, length_bound + 1, In, Out, model
    # )

    blocked, boolean_circuit = BooleanForwardCircuit(
        n, length_bound + 1, In, model
    )

    prob_circuit = None

    # if with_probs:
    #     prob_circuit = ForwardProbabilities(
    #         n, length_bound + 1, In, Out, blocked, model
    #     )

    return blocked, boolean_circuit, prob_circuit


def ForwardCircuit(blocked: gp.MVar, max_length: int, In: Walk, Out: Walk, solver: gp.Model):
    # Circuit simulating or-and adjacency matrix of
    #   the induced subgraph given by V' = V - blocked.
    N, T = blocked.size, max_length
    P = solver.addMVar(shape=(T, N), lb=0.)

    # Enforce first row to be blocked complement and
    #   also toggle off the column of every blocked vertex.
    solver.addConstr(P[-1].sum() == 0, name='nilpotence')
    solver.addConstr(P[0] == 1 - blocked, 'blocked_complement')
    solver.addConstr(P.sum(axis=0) <= T * P[0], f'circuit_blocking')

    for t in trange(T):
        for v in range(N):
            # An active node is outwardly adjacent to only active and blocked nodes.
            #   These constraints are vacuous for any P[t, v] = 0 since all variables are non-negative.
            if t < T-1:
                solver.addConstr(
                    gp.quicksum(
                        P[t+1, w] + blocked[w] for w in Out(v)
                    ) >= len(Out(v)) * P[t, v], f'out_toggle_{t}_{v}'
                )

            # A node is inactive when inwardly adjacent to only inactive nodes.
            #    The constraint is vacuous when any inwardly adjacent node is active.
            if t > 0:
                solver.addConstr(
                    P[t, v] <= gp.quicksum(
                        P[t-1, u] for u in In(v)
                    ), f'in_toggle_{t}_{v}'
                )

    return P


def BooleanForwardCircuit(N: int, T: int, In: Walk, solver: gp.Model):
    # Circuit simulating e, Ae, A²e, ..., Aⁿe, where A is the induced adjacency matrix.
    P = solver.addMVar(shape=(T+1, N), vtype=GRB.BINARY)

    for v in range(N):
        solver.addGenConstrIndicator(
            P[0, v], False,
            len(In(v)) * P[1, v] >= 2 - P[0, In(v)].sum()
        )

    for t in trange(2, T+1):
        for v in range(N):
            solver.addGenConstrIndicator(
                P[t-1, v], True,
                len(In(v)) * P[t, v] >= P[t-1, In(v)].sum()
            )

    return P[0], P[1:]


def ForwardProbabilities(
        N: int, L: int, In: Walk, Out: Walk, blocked: gp.MVar, model: gp.Model
) -> tuple[gp.MVar, gp.MVar]:

    # Path circuit variables
    P = model.addMVar(
        shape=(L + 2, N), lb=0., name=f'random-walk-circuit'
    )

    # Start your walk from the decycling set with mass 1 at each vertex, requiring an extra layer.
    # This is the same as using a probability vector, without having to track the specific size to divide.
    model.addConstr(P[0] == blocked, name='initial-layer')

    # Mass flowed to a blocked node leaves the system,
    #    with all of its length already accumulated previously.
    model.addConstr(
        P[1:].sum(axis=0) <= N * (1 - blocked), name='blocking'
    )

    # Split mass equally among all outward nodes.
    for t in trange(1, L+2):
        for v in range(N):
            in_flow = P[t-1, In(v)].sum() / len(In(v))
            model.addConstr( P[t, v] <= in_flow, f'in-flow-ub-{t}-{v}' )
            model.addConstr( P[t, v] >= in_flow - N * blocked[v], f'in-flow-lb-{t}-{v}' )

    return P
