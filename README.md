# Isogeny Database with Endomorphism Rings

This is a database of isogeny graphs of supersingular elltipic curves over $\mathbb{F}_{p^2}$ for primes $p \leq 29989$, $2 \leq l \leq 11$, including the endomorphism rings of every curve as a maximal order in a quaternion algebra. We include some helpful code for working with small prime examples in [Sagemath](https://www.sagemath.org/), including Deuring correspondence functions like `curve_to_endring(E)`, and `endring_to_curve(O)`.

The database is an extension of Florit and Finol's [Isogeny Database](https://isogenies.enricflorit.com/) to include endomorphism rings.

I am in the process of writing a document explaining how these computations were made, and hope to upload it soon.

**Note:** Rational maps of endomorphisms are not included, only quaternion bases. Also, a choice of Galois conjugate has been made. If you use the wrong conjugate, you'll need to conjugate the assignment of endomorphism rings in the graphs.


## Using it in Sagemath

Clone the repository with:
```
git clone https://github.com/jtcc2/isogeny-database-with-end-rings.git
```

Ensure the `database_utils.sage` is available to Sage, for example in the Sage installation directory or your current working directory.

Edit the `database_utils.sage` and set the `database_root` at the top of the file. This should be the path to the `isogeny_database` directory. For example, if on Windows and using Sage through WSL with the database in your documents this should be:
```
database_root = "/mnt/c/Users/<username>/Documents/isogeny_database/"
```

In Sage, import the functions with `load("database_utils.sage")`. For example:
```python
load("database_utils.sage")

p = 5569
G, js, endrings = get_isog_graph(p, [2,3], True, True)
```


See the functions available below.

## Documentation

The following functions are available in `database_utils.sage`:  
- `get_isog_graph(p, L_list, include_jinvariants=True, include_endrings=True, Fp2=None, B=None, degrees_as_edge_weights=True)`  
Returns an $L$-isogeny graph together with $j$-invariants and endomorphism rings, for example:  
`G, js, endrings = get_isog_graph(5569, [2,3], True, True)`  
They are ordered such that for index `t`, graph vertex `G.vertices()[t]` represents curves with $j$-invariant `js[t]` and endomorphism ring `endrings[t]`.  
- `get_all_supersingular_j_invariants(p)` returns all supersingular $j$-invariants for a given prime.  
- `get_all_maximal_orders(p)` returns all maximal quaternion orders for a given prime.  
- Functions `curve_to_endring(E)` and `endring_to_curve(O)` work as you'd expect for switching between curves and their endomorphism rings.
- Similarily for working with $j$-invariants instead of curves use `j_to_endring(p, j)` and `endring_to_j(O)`.  

Note that if the methods above are given $p$ or $L$ that is not contained in the database, they will compute the required isogeny graph with endomorphism rings. The computed graph is not stored internally and will be recomputed on every function call.

## Database Format

If you want to work with the database files directly without using Sage, it follows the same directory structure as Florit's Isogeny Database. Each prime $p$ has a directory, within which there is a `{p}_nodes.txt` file listing all supersingular $j$-invariants for that prime, and `{p}_{l}.npz` files which store the adjacency matrices of $l$-isogeny graphs in a compressed format.

The difference is the addition of endomorphism ring files. For each prime you'll find a `{p}_endrings.txt` file. Each line of this file represents an endomorphism ring. Lines are ordered so the ith endomorphism ring corresponds to the ith $j$-invariant in the `{p}_nodes.txt` file. An endomorphism ring is represented as a array of basis elements of a maximal quaternion order, each element being listed in terms of coefficients of $1, i,j,k$. The quaternion algebra representation is the same as Sagemath's function `B.<i,j,k>=QuaternionAlgebra(p)`. All orders are stored using the basis of the right order of the equivalent connecting left-ideal with smallest norm from `B.maximal_order()`.

## Recomputing Database

Code to recompute Florit's isogeny database is available [here](https://github.com/gfinol/IsogenyGraph). Alternatively there are other options for computing isogeny graphs on elliptic curves, such as [Sage's built-in algorithm](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_field.html#sage.schemes.elliptic_curves.ell_field.EllipticCurve_field.isogeny_ell_graph) or using a faster [Pari implementation](https://github.com/JamesRickards-Canada/Isogeny).

Assuming you are starting from Florit's database, the code to compute endomorphism rings and attach them to elliptic curves is given in the `compute_database.sage` file. It uses the same `database_root` as `database_utils.sage`. Also there a few primes where [constructive Deuring correspondence](https://github.com/friends-of-quaternions/deuring) code is needed.

## License

The underlying isogeny database is released under the [Open Data Commons Attribution License](https://opendatacommons.org/licenses/by/1-0/), with details [here](https://isogenies.enricflorit.com/about.html).