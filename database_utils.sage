# Set the default path to the `/isogeny_database` directory
database_root = "/mnt/c/Users/<username>/Documents/isogeny_database/"

## Reading from the database

import os

def get_file_content(path):
    """
    Returns the content of a text file.

    INPUT
    - ``path`` -- Path to the text file.
    """
    text_file = open(path, "r")
    data = text_file.read()
    text_file.close()
    return data

def database_get_primes_list():
    """
    Returns list of primes available in isogeny database.
    """
    primes = [ZZ(p) for p in os.listdir(database_root) if not p.endswith(".txt")]
    primes.sort()
    return primes

def database_get_l_list(p):
    """
    Returns a list of small primes l where the l-isogeny graph is available in the isogeny database for characteristic p.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph. Restricted to p < 30,000.
    """
    if not os.path.isdir(database_root + str(p) + "/"):
        raise Exception("Prime " + str(p) + " not found in database.")
    ells = [ZZ(f[len(str(13))+1:-4]) for f in os.listdir(database_root + str(13) + "/") if f.endswith(".npz")]
    ells.sort()
    return ells

def database_get_adj_matrix(p, l):
    """
    Returns the adjacency matrix for l-isogeny graph with characteristic p from the isogeny database.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph. Restricted to p < 30,000.
    - ``l`` -- Small prime within set {2,3,5,7,11}. The degree of isogenies in the l-isogeny graph
    """
    import scipy.sparse
    if not os.path.isdir(database_root + str(p) + "/"):
        raise Exception("Prime " + str(p) + " not found in database.")
    path = database_root + str(p) + "/" + str(p) + "_" + str(l) + ".npz"
    if not os.path.isfile(path):
        raise Exception("The " + str(l) + "-isogeny graph is not in the database for prime" + str(p) + ".")
    sparse_matrix = scipy.sparse.load_npz(path)
    sparse_matrix.toarray()
    return matrix(ZZ, sparse_matrix.toarray())

def Fp2_evaluate_isomorphism(elt, Fp2):
    """
    Maps an element of a finite field GF(p^2) into another finite field GF(p^2) where the generators are different.

    INPUT
    - ``elt`` -- element of a finite field GF(p^2).
    - ``Fp2`` -- a finite field GF(p^2) with different generator.
    """
    return elt[0] + Fp2.gens()[0]*elt[1]

def database_get_j_invariants(p, Fp2=None):
    """
    Returns list of j-invariants of vertices of isogeny graphs for a prime p in the isogeny database.

    INPUT
    - ``p`` -- A prime p to fetch the supersingular j-invariants of. Restricted to p < 30,000.
    - ``Fp2`` -- Finite field GF(p^2). Returns j-invariants as elements of this field. If None, creates a new finite field.
    OUTPUT
    - List of j-invariants as elements of Fp2. Ordered corresponding to the adjacency matrix in the isogeny database.
    """
    if not os.path.isdir(database_root + str(p) + "/"):
        raise Exception("Prime " + str(p) + " not found in database.")
    filepath = database_root + str(p) + "/" + str(p) + "_nodes.txt"
    if not os.path.isfile(filepath):
        raise Exception("The j-invariants for prime " + str(p) + "are not in the isogeny database.")
    contents = get_file_content(filepath)
    j_invs_str = contents.splitlines()
    Fp2.<z> = GF(p^2)
    j_invariants = [Fp2(j) for j in j_invs_str]
    return j_invariants if Fp2 == None else [Fp2_evaluate_isomorphism(j, Fp2) for j in j_invariants]

def comma_split_string(content):
    """
    Splits string by commas, ignoring those contained within open square or round brackets.

    INPUT
    - ``content`` -- string to split.
    """
    start_pos = 0
    arr = []
    open_brackets = 0
    for i in range(0, len(content)):
        if content[i] in '([': open_brackets += 1
        if content[i] in ')]': open_brackets -= 1
        if content[i] == ',' and open_brackets == 0:
            arr.append(content[start_pos:i])
            start_pos = i+1
    if start_pos != i+1: arr.append(content[start_pos:])
    return arr

def str_to_arr(content, dim=1):
    """
    Converts a string (an array printed out) back to an array.
    Supports arrays, tuples, and higher-dim arrays with content either strings or None.

    INPUT
    - ``content`` -- String in the format of an array or tuple, e.g. "[1,2,3]", or "('a', 'b')".
    - ``dim`` -- Dimension of the array represented by string, e.g. "[1,2]" is dimension 1, "[[1],[2]]" is dimension 2.
    """
    if content == '[]' or content == '()': return []
    split = comma_split_string(content[1:-1])
    if dim > 1:
        return [str_to_arr(c.strip(), dim-1) for c in split]
    else:
        for i in range(0, len(split)):
            stripped = split[i].strip()
            if stripped == "None":
                split[i] = None
            else:
                if stripped[0] in "(['":
                    split[i] = stripped[1:-1]
                else:
                    split[i] = stripped
        return split

def quaternion_algebra_is_isomorphic(B1, B2, certificate=False):
    """
    Are quaternion algebras B1 and B2 are isomorphic.

    INPUT
    - ``B1`` -- Domain of isomorphism if it exists.
    - ``B2`` -- Codomain of isomorphism if it exists.
    - ``certificate`` -- If True, returns a map from B1 to B2.
    OUTPUT
    - Boolean indicating if the quaternion algebras are isomorphic or not
    - (if certificate=True) Map from B1 to B2

    Source: "Finding Orientations of Supersingular Elliptic Curves and Quaternion Orders",
            https://github.com/jtcc2/finding-orientations/blob/d3668840844d1a69409e78b71b6c49954fda82ee/algorithm_5_2.py#L36
    """

    if B1.invariants()[1] != B2.invariants()[1]:
        raise NotImplementedError("Quaternion algebra isomorphism only implemented for B1.invariants()[1] == B2.invariants()[1].")
    if B1 == B2:
        if certificate: return True, (lambda elt: B2(elt))
        return True
    q_old = -B1.invariants()[0]
    i, j, k = B2.gens()
    q, p = -B2.invariants()[0], -B2.invariants()[1]
    try:
        x, y = DiagonalQuadraticForm(QQ, [1,p]).solve(q_old/q)
    except:
        if certificate: return False, None
        return False
    if certificate:
        gamma = x + j*y
        return True, (lambda elt: sum([coeff*b for coeff, b in zip(elt.coefficient_tuple(), [1, i*gamma, j, k*gamma])]))
    return True

def database_get_maximal_orders(p, B=None):
    """
    Returns list of maximal quaternion orders in the quaternion algebra ramified at p and infinity from the isogeny database.

    INPUT
    - ``p`` -- Prime p, to represent the quaternion algebra ramified at p and infinity. Restricted to p < 30,000.
    - ``B`` -- Quaternion algebra ramified at p and infinity which the maximal orders lie within. If None, creates one.
    OUTPUT
    - List of maximal quaternion orders contained in B. Ordered corresponding to the adjacency matrix in the isogeny database.
    """
    if not os.path.isdir(database_root + str(p) + "/"):
        raise Exception("Prime " + str(p) + " not found in database.")
    file_path = database_root + str(p) + "/" + str(p) + "_endrings.txt"
    if not os.path.isfile(file_path):
        raise Exception("The endomorphism rings for prime " + str(p) + " do not exist in the isogeny database.")
    B_database = QuaternionAlgebra(p)
    i,j,k = B_database.gens()
    contents = get_file_content(file_path)
    O__str_basis = contents.splitlines()
    O__basis_str = [str_to_arr(basis, dim=2) for basis in O__str_basis]
    order_bases = [[QQ(b[0])+QQ(b[1])*i+QQ(b[2])*j+QQ(b[3])*k for b in basis] for basis in O__basis_str]
    if B == None:
        return [B_database.quaternion_order(basis) for basis in order_bases]
    is_iso, m = quaternion_algebra_is_isomorphic(B_database, B, certificate=True)
    if not is_iso:
        raise Exception("Quaternion algebra not isomorphic to the quaternion algebra ramified at p and infinity.")
    return [B.quaternion_order([m(b) for b in basis]) for basis in order_bases]

def database_get_l_isogeny_graph(p, l, degrees_as_edge_weights=True):
    """
    Returns l-isogeny graph as Sage Graph object.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph. Restricted to p < 30,000.
    - ``l`` -- Small prime used as the degrees of edges of the l-isogeny graph. From the set {2,3,5,7,11}.
    - ``degrees_as_edge_weights`` -- If True the edges are given edge weights l.
    """
    A = database_get_adj_matrix(p, l)
    G = Graph(A, format='adjacency_matrix', loops=True, multiedges=True)
    if degrees_as_edge_weights:
        G = Graph([G.vertices(), [(e[0], e[1], l) for e in G.edges()]], loops=True, multiedges=True)
    return G

def database_get_isog_graph(p, L_list, include_jinvariants=True, include_endrings=True, Fp2=None, B=None, degrees_as_edge_weights=True):
    """
    Returns L-isogeny graph as Sage Graph object, optionally with j_invariants and endomorphism rings.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph. Restricted to p < 30,000.
    - ``L_list`` -- List of distinct small primes within set {2,3,5,7,11}. The isogeny degrees in the L-isogeny graph.
    - ``include_jinvariants`` -- If True outputs list of j-invariants corresponding to vertices in the isogeny graph.
    - ``include_endrings`` -- If True outputs list of endomorphism rings as maximal quaternion orders, corresponding to vertices in the isogeny graph.
    - ``Fp2`` -- Finite field GF(p^2). Returns j-invariants as elements of this field. If None, creates a new finite field.
    - ``B`` -- Quaternion algebra ramified at p and infinity which the maximal orders lie within. If None, creates one.
    - ``degrees_as_edge_weights`` -- If True the output graph object will have isogeny degrees as edge weights.
    OUTPUT
    - Sage Graph object representing the L-isogeny graph.
    - If include_jinvariants=True, outputs List of j-invariants as elements of Fp2. Ordered corresponding to graph vertices.
    - If include_endrings=True, outputs List of maximal quaternion orders contained in B. Ordered corresponding to graph vertices.
    """
    G = database_get_l_isogeny_graph(p, L_list[0], degrees_as_edge_weights)
    for i in range(1, len(L_list)):
        G.add_edges(database_get_l_isogeny_graph(p, L_list[i], degrees_as_edge_weights).edges())
    if not include_jinvariants and not include_endrings:
        return G
    if include_jinvariants:
        j_invariants = database_get_j_invariants(p)
    if include_endrings:
        endrings = database_get_maximal_orders(p, B)
    if include_jinvariants and not include_endrings: return G, j_invariants
    if not include_jinvariants and include_endrings: return G, endrings
    return G, j_invariants, endrings

def small_equivalent_ideal(I):
    """
    Returns a left-ideal J of smaller norm in the left ideal class.

    Source: "Finding Orientations of Supersingular Elliptic Curves and Quaternion Orders",
            https://github.com/jtcc2/finding-orientations/blob/d3668840844d1a69409e78b71b6c49954fda82ee/ideals.py#L4C1-L12C16
    """
    _,mn = I.quadratic_form().__pari__().qfminim(None,None,1)
    el = sum(ZZ(c)*g for c,g in zip(mn, I.basis()))
    y = el.conjugate() / I.norm()
    I *= y
    return I

def connecting_ideal(O0, O1):
    """
    Returns the connecting ideal between two quaternion orders of the same quaternion algebra.
    """
    I = O0 * O1
    I = I * denominator(I.norm())
    return I

##############################################################################################
## Failover and helper functions

def get_isog_graph(p, L_list, include_jinvariants=True, include_endrings=True, Fp2=None, B=None, degrees_as_edge_weights=True):
    """
    For any p and L, returns L-isogeny graph as Sage Graph object, optionally with j_invariants and endomorphism rings. Gets graph from the database if it exists,
    otherwise computes the additional required information.

    INPUT
    - ``p`` -- Prime characteristic of base field of the curves in the isogeny graph. No restriction.
    - ``L_list`` -- List of distinct small primes. The isogeny degrees in the L-isogeny graph. No restriction.
    - ``include_jinvariants`` -- If True outputs list of j-invariants corresponding to vertices in the isogeny graph.
    - ``include_endrings`` -- If True outputs list of endomorphism rings as maximal quaternion orders, corresponding to vertices in the isogeny graph.
    - ``Fp2`` -- Finite field GF(p^2). Returns j-invariants as elements of this field. If None, creates a new finite field.
    - ``B`` -- Quaternion algebra ramified at p and infinity which the maximal orders lie within. If None, creates one.
    - ``degrees_as_edge_weights`` -- If True the output graph object will have isogeny degrees as edge weights.
    OUTPUT
    - Sage Graph object representing the L-isogeny graph.
    - If include_jinvariants=True, outputs List of j-invariants as elements of Fp2. Ordered corresponding to graph vertices.
    - If include_endrings=True, outputs List of maximal quaternion orders contained in B. Ordered corresponding to graph vertices.
    """
    def get_start_curve_j(p, Fp2):
        """Gives j-invariant of supersingular elliptic curve over Fp2"""
        q = next(q for q in Primes() if q%4 == 3 and kronecker_symbol(-q,p) == -1)
        K = QuadraticField(-q)
        H = K.hilbert_class_polynomial()
        return H.change_ring(Fp2).any_root()
    
    L_list_remaining = L_list
    G = None
    j_invs = None
    endrings = None
    endrings_BM = None
    BM = None
    shown_warning = False
    if p <= 29989:
        database_ls = database_get_l_list(p)
        L_list_remaining = []
        L_list_database = []
        for l in L_list:
            if l in database_ls:
                L_list_database.append(l)
            else:
                L_list_remaining.append(l)
        if len(L_list_remaining) == 0:
            return database_get_isog_graph(p, L_list_database, include_jinvariants, include_endrings, Fp2, B, degrees_as_edge_weights=degrees_as_edge_weights)
        G, j_invs, endrings = database_get_isog_graph(p, L_list_database, True, True, Fp2, None, degrees_as_edge_weights=degrees_as_edge_weights)
    else:
        if shown_warning == False:
            # TODO: we can fix this issue by being more careful with special curves
            shown_warning = True
            print("UserWarning: Computing graph data from outside isogeny database. Note there might be inconsistencies with loops/multi-edges at special curves with extra automorphisms.")
        if Fp2 == None:
            Fp2.<z> = GF(p^2)
        E_start = EllipticCurve(Fp2, j=get_start_curve_j(p, Fp2))
        l = L_list_remaining[0]
        G = E_start.isogeny_ell_graph(l, directed=False, label_by_j=True)
        j_invs = [Fp2(v) for v in G.vertices()]
        edgelabel = l if degrees_as_edge_weights else None
        new_edges = [(G.vertices().index(e[0]), G.vertices().index(e[1]), edgelabel) for e in G.edges()]
        G = Graph([list(range(0, len(G.vertices()))), new_edges], loops=True, multiedges=True)
        L_list_remaining.remove(l)
        if len(L_list_remaining) == 0:
            if not include_jinvariants and not include_endrings: return G
            if include_jinvariants and not include_endrings: return G, j_invs
        # Compute Brandt Module graph, isomorphism from the G to it, and use this to sort maximal orders
        Gcopy = copy(G)
        Gcopy.remove_multiple_edges()
        BM = BrandtModule(p)
        O0 = BM.maximal_order()
        endrings_BM = [small_equivalent_ideal(connecting_ideal(O0, I.left_order())).right_order() for I in BM.right_ideals()] # wrong ordering
        G_endring = Graph(BM.hecke_matrix(l), format='adjacency_matrix', loops=True, multiedges=True)
        G_endring.remove_multiple_edges()
        isIso, m = Gcopy.is_isomorphic(G_endring, certificate=True)
        if not isIso:
            raise Exception("Brandt module graph not isomorphic. This should never happen.")
        endrings = [endrings_BM[ZZ(m[v])] for v in G.vertices()] # corrected ordering
        # We ignore issue of graph automorphisms as for large graphs outside isogeny database this is very rare
        known_automorphisms = 1 if j_invs[-1] in GF(p) else 2
        actual_automorphisms = Gcopy.automorphism_group().cardinality()
        if actual_automorphisms > known_automorphisms:
            # TODO: As in `compute_database.sage` we could fix this issue using constructive Deuring correspondence
            print("UserWarning: Computing graph data from outside isogeny database. The parameters you provided are a rare case when additional graph automorphisms exist. This means endomorphism rings could be ordered incorrectly. Use Deuring correspondence to check this.")
        if len(L_list_remaining) == 0:
            if include_endrings and B != None:
                if B != endrings[0].quaternion_algebra():
                    _, m = quaternion_algebra_is_isomorphic(endrings[0].quaternion_algebra(), B, certificate=True)
                    endrings = [B.quaternion_order([m(b) for b in O.basis()]) for O in endrings]
            if not include_jinvariants and include_endrings: return G, endrings
            if include_jinvariants and include_endrings: return G, j_invs, endrings
    
    # Need to add in missing edges for larger prime degree isogenies that are not in the database
    if len(L_list_remaining) > 0:
        if shown_warning == False:
            shown_warning = True
            print("UserWarning: Computing graph data from outside isogeny database. Note there might be inconsistencies with loops/multi-edges at special curves with extra automorphisms.")
        new_G_edges = []
        for l in L_list_remaining:
            if BM == None:
                BM = BrandtModule(p)
            if endrings_BM == None:
                O0 = BM.maximal_order()
                endrings_BM = [small_equivalent_ideal(connecting_ideal(O0, I.left_order())).right_order() for I in BM.right_ideals()]
            G_endring = Graph(BM.hecke_matrix(l), format='adjacency_matrix', loops=True, multiedges=True)
            # Want to add edges of G_endring to G, noting different degree, and ordering of orders
            edgelabel = l if degrees_as_edge_weights else None
            new_G_edges.extend([(endrings.index(endrings_BM[e[0]]), endrings.index(endrings_BM[e[1]]), edgelabel) for e in G_endring.edges()])
        new_G_edges.extend(G.edges())
        G = Graph([G.vertices(), new_G_edges], loops=True, multiedges=True)
    if include_endrings and B != None:
        if B != endrings[0].quaternion_algebra():
            _, m = quaternion_algebra_is_isomorphic(endrings[0].quaternion_algebra(), B, certificate=True)
            endrings = [B.quaternion_order([m(b) for b in O.basis()]) for O in endrings]
    if not include_jinvariants and not include_endrings: return G
    if include_jinvariants and not include_endrings: return G, j_invs
    if not include_jinvariants and include_endrings: return G, endrings
    return G, j_invs, endrings

def j_to_endring(p, j, B=None):
    """
    Returns the endomorphism ring as a maximal quaternion order corresponding to a given j-invariant.

    INPUT
    - ``p`` - The prime characteristic.
    - ``j`` - The j-invariant to find the endomorphism ring of.
    """
    _, j_invs, orders = get_isog_graph(p, [2], Fp2=j.parent(), B=B)
    return orders[j_invs.index(j)]

def curve_to_endring(E, B=None):
    """
    Returns the endomorphism ring as a maximal quaternion order corresponding to a given curve E.

    INPUT
    - ``E`` - The curve to find the endomorphism ring of.
    """
    p = E.base_field().characteristic()
    return j_to_endring(p, E.j_invariant(), B)

def endring_to_j(O, Fp2=None):
    """
    Returns j-invariant corresponding to given endomorphism ring.

    INPUT
    - ``O`` -- Maximal order in quaternion algebra ramified at p and infinity.
    - ``Fp2`` -- Finite field Fp2. If None creates a new Fp2 object.
    OUTPUT
    - j-invariant as element of Fp2, where curves of that j-invariant have endomorphism ring isomorphic to O.
    """
    _, j_invariants, endrings = get_isog_graph(p, [2], Fp2=Fp2)
    O0 = endrings[0].quaternion_algebra().maximal_order()
    if O.quaternion_algebra() != endrings[0].quaternion_algebra():
        _, m = quaternion_algebra_is_isomorphic(O.quaternion_algebra(), endrings[0].quaternion_algebra(), certificate=True)
        O = endrings[0].quaternion_algebra().quaternion_order([m(b) for b in O.basis()])
    J = small_equivalent_ideal(connecting_ideal(O0, O))
    O = J.right_order()
    return j_invariants[endrings.index(O)]

def endring_to_curve(O, Fp2=None):
    """
    Returns a supersingular elliptic curve with given maximal quaternion order as it's endomorphism ring.

    INPUT
    - ``O`` -- Maximal quaternion order representing the endomorphism ring of the curve.
    """
    j = endring_to_j(O, Fp2)
    return Elliptic_curve_from_j(j)

def get_all_supersingular_j_invariants(p, Fp2=None):
    """
    Returns all supersingular j-invariants in Fp2 for a prime p.

    INPUT
    - ``p`` -- Prime characteristic.
    - ``Fp2`` -- Finite field of degree 2, characteristic p, to return j-invariants within. If None, creates new object.
    """
    if p <= 29989:
        return database_get_j_invariants(p, Fp2)
    _, j_invs = get_isog_graph(p, [2], include_jinvariants=True, include_endrings=False, Fp2=Fp2, B=None, degrees_as_edge_weights=False)
    return j_invs

def get_all_maximal_orders(p, B=None):
    """
    Returns all maximal orders of quaternion algebra ramified at p and infinity.

    INPUT
    - ``p`` -- Prime p.
    - ``B`` -- Representation of the quaternion algebra ramified at p and infinity. Returns orders in this representation. If None, creates new object.
    """
    if p <= 29989:
        return database_get_maximal_orders(p, B)
    _, endrings = get_isog_graph(p, [2], include_jinvariants=False, include_endrings=True, Fp2=None, B=B, degrees_as_edge_weights=False)
    return endrings

