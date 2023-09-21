# This is the code used to add endomorphism rings to Florit's isogeny database
#    generating the endomorphism ring files for every prime.

# Sage may crash several times. Rerun the script and it will continue where it left off.

# To use this download the isogeny database without endomorphism ring data from:
#    https://isogenies.enricflorit.com/
# And set the database_root in 'database_utils.sage'

# Note we require constructive Deuring correspondence code. This is very rarely used when isogeny graphs have many automorphisms, actually only for p=251
try:
    from deuring.broker import starting_curve
    from deuring.correspondence import constructive_deuring
    deuring_correspondence_available = True
except:
    print("Warning: Could not find constructive Deuring correspondence code. "
        + "This is only used in special cases, specifically only when p=19,31,103,251. Without it, these values will be skipped. "
        + "To fix this, git clone the following repository into the directory you are running Sage: https://github.com/friends-of-quaternions/deuring")
    deuring_correspondence_available = False

import os

load("database_utils.sage")


def output_file(file_path, content):
    """
    Outputs a file to a given path.

    INPUT
    - ``file_path`` -- Path to output text file.
    - ``content`` -- Text to output to the file.
    """
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w') as f:
        f.write(content)

def output_quaternion_orders(p, orders):
    """
    Given the quaternion orders for prime p, outputs them to the correct file in the isogeny database.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph.
    - ``orders`` -- Quaternion orders corresponding to endomorphism rings of vertices in isogeny graph. Should be ordered correctly.
    """
    if not os.path.isdir(database_root + str(p) + "/"):
        raise Exception("Prime " + str(p) + " not found in database.")
    file_path = database_root + str(p) + "/" + str(p) + "_endrings.txt"
    content = ""
    for r in range(0, len(orders)):
        O = orders[r]
        if r != 0:
            content += os.linesep
        content += str([b.coefficient_tuple() for b in O.basis()])
    output_file(file_path, content)

def get_endring_graph(p, L_list):
    """
    Returns L-ideal graph as sage Graph object, and list of maximal quaternion orders corresponding to vertices.

    INPUT
    - ``p`` -- Prime defining the quaternion algebra whose graph of maximal orders is returned.
    - ``L_list`` -- List of small primes. The graph includes an edge/ideal if it's norm is in this list.
    OUTPUT
    - Sage graph object representing the graph of maximal quaternion orders.
    - List of maximal quaternion orders which align with the vertices of the graph,
        these are in the unique form of the right order of the smallest connecting left-ideal to a fixed O0.
    """
    BM = BrandtModule(p)
    O0 = BM.maximal_order()
    G_adj = sum([BM.hecke_matrix(l) for l in L_list])
    G = Graph(G_adj, format='adjacency_matrix', loops=True, multiedges=True)
    orders = [I.left_order() for I in BM.right_ideals()]
    orders = [small_equivalent_ideal(connecting_ideal(O0, O)).right_order() for O in orders]
    return G, orders

def check_orders_conjugate(O1, O2):
    """
    Check two quaternion orders of the same quaternion algebra are equivalent.
    """
    O0 = O1.quaternion_algebra().maximal_order()
    I1 = connecting_ideal(O1, O0)
    I2 = connecting_ideal(O2, O0)
    return I1.is_equivalent(I2)




print("starting")
for p in database_get_primes_list():
    file_path = database_root + str(p) + "/" + str(p) + "_endrings.txt"
    if os.path.isfile(file_path):
        print("skipping " + str(p) + " as endrings.txt already exists")
        continue
    print(p)
    l = database_get_l_list(p)[0]
    G, j_invariants = database_get_isog_graph(p, [l], True, False)
    G_endring, orders = get_endring_graph(p, [l])
    G.remove_multiple_edges()
    G_endring.remove_multiple_edges()
    # Find an isomorphism - note this could be an incorrect isomorphism
    isIso, m = G.is_isomorphic(G_endring, certificate=True)
    if not isIso:
        # This should never happen
        print("Error: Graphs not isomorphic for prime " + str(p))
        continue
    # Sort the orders to align with j_invariants
    orders = [orders[ZZ(m[v])] for v in G.vertices()]
    # If the graph has additional automorphisms we need to make sure we took the right isomorphism m above
    A = G.automorphism_group()
    known_automorphisms = 1 if j_invariants[-1] in GF(p) else 2
    actual_automorphisms = A.cardinality()
    if actual_automorphisms > known_automorphisms:
        if not deuring_correspondence_available:
            print("Skipping prime " + str(p) + ", requires constructive Deuring.")
            continue
        E0, iota, O0 = starting_curve(j_invariants[0].parent())
        E0_idx = j_invariants.index(E0.j_invariant())
        # Put O0 in the correct quaternion algebra representation
        _, m = quaternion_algebra_is_isomorphic(O0.quaternion_algebra(), orders[0].quaternion_algebra(), certificate=True)
        O0 = orders[0].quaternion_algebra().quaternion_order([m(b) for b in O0.basis()])
        # Check graph automorphisms which assign O0 correctly to E0
        valid_automorphisms = []
        for automorphism in A:
            new_vertex_order = [automorphism(v) for v in G.vertices()]
            order_idx = new_vertex_order.index(E0_idx)
            O0_should_be = orders[order_idx]
            if check_orders_conjugate(O0, O0_should_be):
                valid_automorphisms.append(automorphism)
        if len(valid_automorphisms) != known_automorphisms:
            # We still have multiple options for automorphisms, use constructive deuring to find which is the correct one.
            # This only happens for p=251.
            #    Note from now on we are looking for one specific automorphism and will not accept its galois conjugate,
            #    as constructive deuring fixes a galois conjugate consistently.
            fixed_tuples = []
            while len(valid_automorphisms) > 1:
                # Pick a random vertex that it mapped to a different vertex in an automorphism
                O_index = -1
                while O_index == -1:
                    rand_auto = valid_automorphisms[randrange(0, len(valid_automorphisms))]
                    if rand_auto.cycle_tuples() != []:
                        rand_tuple = rand_auto.cycle_tuples()[randrange(0, len(rand_auto.cycle_tuples()))]
                        if rand_tuple not in fixed_tuples:
                            O_index = rand_tuple[0]
                            fixed_tuples.append(rand_tuple)
                # Check where that vertex should actually be via constructive Deuring
                O = orders[O_index]
                I = connecting_ideal(O0, O)
                E, _, _ = constructive_deuring(I, E0, iota)
                j = E.j_invariant()
                # Filter out automorphisms where O does not correspond to j or j^p
                new_valid_automorphisms = []
                for automorphism in valid_automorphisms:
                    new_vertex_order = [automorphism(v) for v in G.vertices()]
                    j_inv_new = j_invariants[new_vertex_order[O_index]]
                    if j == j_inv_new:
                        new_valid_automorphisms.append(automorphism)
                valid_automorphisms = new_valid_automorphisms
        final_automorphism = valid_automorphisms[0]
        # Reorder the list of orders to coincide with vertices/js in normal order
        new_vertex_order = [final_automorphism(v) for v in G.vertices()]
        new_orders_list = [None for O in orders]
        for O_index in range(0, len(orders)):
            new_orders_list[new_vertex_order[O_index]] = orders[O_index]
        orders = new_orders_list
    # Output the orders to file
    output_quaternion_orders(p, orders)

    # There can sometimes be a issue with memory leakage causing Sage to crash,
    #   deallocating the memory used by the largest variables seems to help
    del(G)
    del(G_endring)
    del(orders)