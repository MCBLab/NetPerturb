# Slide 3 — Introduction (concept image), landing on slide 4's question

Visual: `Fig/concept_grn_knockout.png`
Next slide: "Can we run a fully in-silico experiment to determine good candidates
using untreated scRNAseq data? How can we make it fast and affordable?"

---

Single-cell data gives us cell identities — and for each one we can build a gene
regulatory network: a wiring diagram of which genes influence which, inferred from
how they vary together across that identity's cells.

Once you have that wiring, you can ask a question of it. Take a drug target, cut
its wires, and see what happens to the network.

In one cell type the target is a hub: cutting it breaks the network apart, and the
genes at the end of its strongest edges move the most. In another, the same target
sits at the periphery — the wiring holds, almost nothing moves.

That difference is the measurement. A high score means the identity's network
leaned on that target; a low score means it didn't. And the genes that moved are
the knockout's readout.

The important part is that the score belongs to the **cell identity**, not to the
gene. The same target can carry one population's network and barely touch
another's — which is exactly what you want to know before choosing a candidate.

*[to slide 4]* So the question becomes: how far can we take this?

---

## Delivery

- 45–55 seconds at a normal pace.
- Say "cut its wires" while pointing at the red edges in the concept image.
- Let the last line hang before advancing — slide 4's question then answers itself.
