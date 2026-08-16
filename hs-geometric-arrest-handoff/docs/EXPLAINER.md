# What this project is about, in plain terms

## The puzzle

Pour sand into a jar and shake it. At some point it stops flowing and behaves
like a solid, even though nothing has crystallised — the grains are still
arranged messily. The same thing happens to microscopic plastic spheres
suspended in water: crowd them enough and the suspension stops flowing and
becomes a glass. Nobody fully understands why this happens at the particular
density it happens at.

De Graaf's 2D paper makes an unusual proposal. Instead of looking for a local
structural motif that appears near the transition, he argues that arrest happens
at a density picked out by *pure geometry*.

## The 2D idea

Give every particle its own territory: the region of space closer to it than to
any other particle. That is a Voronoi cell, and the cells tile the plane
completely.

Now here is the neat bit. Suppose the disks are all the same size, and suppose
the arrangement is special in a particular way: every cell has the disk perfectly
inscribed in it, touching all its edges. Then the packing fraction of the disks
is *exactly equal* to a shape number of the cell — how circle-like the cell is,
$q = 4\pi A/P^2$. Density becomes readable off shape.

There is a specific tiling of the plane by pentagons — the "floret pentagonal
tiling", which you get by taking a close-packed triangular lattice and removing
one particle in every seven — where this number comes out to
$\sqrt3\pi/7 \approx 0.777$. And that is very close to where experiments see
2D colloids fall out of equilibrium, and where 2D granular disks jam. His claim
is that this is not a coincidence: the disordered fluid becomes "entangled" with
this ordered geometric reference state, and that is what arrests it.

His 3D section is a sketch. He builds some candidate structures out of FCC, BCC
and diamond arrangements, gets some numbers, and says the connection should be
looked at more closely.

## What I did — part one: finish the geometry

**The 3D version of the key identity.** The 2D statement relies on a small
geometric fact about polygons with an inscribed circle. The same argument works
for polyhedra with an inscribed sphere, and gives

$$\eta = \frac{36\pi V^2}{S^3}$$

where $\eta$ is the packing fraction, $V$ the cell volume, $S$ its surface area.
The right-hand side is a standard measure of how sphere-like a shape is. So the
correct 3D quantity to measure is this one — not the number of faces, and not
the bond-orientational order parameters people usually reach for.

**The rule this gives you.** The identity only holds if the sphere touches every
single face of its cell. Translated: *every Voronoi neighbour must be a
touching neighbour*. If even one face of the cell belongs to a particle that
isn't touching, the whole argument breaks. That is a sharp admissibility test,
and applying it immediately throws out two of de Graaf's candidates — BCC and
diamond both fail it.

**Finding all the allowed structures.** This looked like a hard search problem
and I initially attacked it as one, with a numerical optimiser. That failed,
repeatedly. The fix was to notice that the quantities being equated are *linear*
in a standard set of lattice parameters (Selling parameters). So it isn't a
search at all — it's a system of linear equations, which you can solve exactly
and completely. The answer:

> There are exactly **three** lattices in 3D whose Voronoi cells satisfy the
> criterion: simple cubic ($\eta = 0.5236$), simple hexagonal ($\eta = 0.6046$),
> and FCC ($\eta = 0.7405$).

And no fourteen-faced solution can exist, which means BCC's failure isn't bad
luck — it's structural.

**A family nobody had noticed — now with a proof.** De Graaf's structures come
from removing particles from a close-packed lattice in a regular pattern. In 3D
the options are exactly: remove 1 in 4, 1 in 5, 1 in 7, or 1 in 13, at packing
fractions 0.5554, 0.5924, 0.6347, 0.6835. This is now *derived*, not just
observed: a simple double-counting of the holes gives the density formula, and
a spectral argument (the number of deleted neighbours per particle can be at
most the magnitude of the most negative eigenvalue of the contact lattice)
proves the list is complete.

The rule behind it: the removal pattern works exactly when each remaining
particle sits next to the same whole number of holes. Run the same argument in 2D and the spectral
bound allows exactly three rungs: 1 in 3, 1 in 4, and 1 in 7. Two of them are
de Graaf's structures — his honeycomb (caging onset) and his floret pentagonal
tiling, whose depleted lattice is known as the maple-leaf lattice. He presents
them as separate constructions; they are rungs of one ladder. The third rung is
the kagome lattice at density 0.680175, a new falsifiable prediction. (An
earlier draft of this project also claimed a fourth rung at 0.4535; the proof
shows that one cannot exist, and I have retracted it.)

**Four things in his paper that are wrong.** All checkable in a few lines of
code:

1. He says a vacancy in a BCC lattice has 4 nearest neighbours of one type. It
   has 8. His resulting number is wrong, and with the right count the volume
   bookkeeping balances exactly — with his count it doesn't.
2. He miscounts how many tetrahedra sit next to an octahedron in one honeycomb
   (2 instead of 8), so his cell only fills three-quarters of space. His claimed
   agreement with a known granular packing density disappears once fixed.
3. He writes that the honeycomb comes from removing two particles in three. It's
   one in three — his own formula in the same paragraph confirms it.
4. He says no candidate structure falls in the range where the 3D glass
   transition is reported (0.58–0.64). Three do.

That last one sounds like good news for his hypothesis. It isn't, and this
matters: the allowed values are packed roughly 0.04 apart across the whole
relevant range, while the experimental glass transition is quoted anywhere from
0.588 to 0.639. Hitting one of them by chance is nearly guaranteed. In 2D his
0.777 was isolated, which is why the match was impressive. In 3D the same kind of
match proves nothing. So the test has to be about *mechanism*, not about matching
a number.

## What I did — part two: test it

I built a hard-sphere Monte Carlo simulation (spheres that can't overlap, moved
around randomly — a reasonable stand-in for how colloids actually jiggle), plus
the analysis machinery to compute Voronoi cells and shape numbers. Both were
validated: the simulation reproduces the known hard-sphere equation of state to
0.04%, and the analysis reproduces the exact geometry numbers to eight decimal
places.

Then I ran a pilot: two very different mixtures of sphere sizes, densities from
0.40 to 0.62.

**Result 1: the arrest density couldn't be pinned down.** The standard method
(extrapolate where diffusion goes to zero) gave 0.612 ± 0.012 for one mixture and
0.626 ± 0.012 for the other. But those error bars are fake — just changing which
data points you fit moves the answer from 0.599 to 0.884. The real uncertainty is
ten times the quoted one. Which means this method can't distinguish between
candidates spaced 0.04 apart, no matter how long you run it.

**Result 2: the most direct test of the mechanism failed — twice.** If the fluid
is being pulled toward one of these special structures, more and more of its
Voronoi neighbours should be *touching* neighbours. But in a fluid of hard
spheres, exact touching essentially never happens, so measured strictly this
number is zero at every density. Measured loosely — counting neighbours within
some small gap — it does grow, but I found the gap distribution has only one
characteristic scale. That means the measurement contains no information beyond a
single average number, and that number changes perfectly smoothly through the
transition. So the observable isn't just showing nothing; it's the wrong kind of
observable.

**Result 3: one thing did work.** Inflate every sphere by a small amount and ask
when the "touching" network first spans the whole box. Call that inflation
$\varepsilon^*$. For the two very different mixtures, $\varepsilon^*$ agrees to
about 4% at every density. That independence-from-composition is exactly the
fingerprint de Graaf uses to argue something is geometric rather than an accident
of particle sizes — and it's the only place in the whole pilot where it shows up.

I couldn't tell whether that curve has a *feature* at any particular density,
because I only had three configurations per point and the noise swamps it.

## Where this leaves things

The geometry is finished, and it's stronger than what was in the paper: an exact
identity, a complete classification, a unifying family, and four corrections.

The physics is open. But the shape of the remaining question changed a lot. The
expensive part of the original plan — long dynamical simulations to find where
diffusion vanishes — turns out to be the part that doesn't work. The part that
does work needs no dynamics at all, just well-equilibrated snapshots. That makes
the decisive experiment roughly an order of magnitude cheaper than planned, and
it reduces to one question:

> Does the composition-independent $\varepsilon^*(\eta)$ curve have a real
> feature, and does it sit at one of the six allowed values?

A "no" would be as useful as a "yes". And there's a deeper worry that no amount
of simulation resolves: a geometric ground state is a single arrangement, so it
has no entropy, and it isn't obvious why a warm jiggling fluid should care about
it at all. De Graaf raises this objection against his own 2D result and doesn't
answer it. In 2D there's at least a known mechanism (Mermin–Wagner) that destroys
long-range order at finite temperature and could plausibly leave a ghost of the
structure behind. In 3D there isn't.


## Two late results that change the picture

**The special states don't have to be crystals.** The counting argument never
uses periodicity. Exhaustively listing every valid hole pattern on a small
patch shows most of them are *not* regular lattices — and yet every one of them
is exactly "tangential" with exactly the special density. So the geometric
reference state is a whole family, including disordered members. This matters
because the strongest objection to the whole idea was: a single perfect crystal
has no entropy, so why would a warm, jiggling fluid care about it? If the family
is large, that objection dissolves.

**HCP muddies the water usefully.** The same construction might also work on
the other close-packed stacking (HCP), at the *same densities* but with
differently shaped cells. If so, measuring a density can never tell you which
structure the fluid is near — you have to measure the cell shapes and their
network. That settles a design question for the simulations: topology
classification isn't optional.
