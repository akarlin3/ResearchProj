#!/usr/bin/env python3
"""Generate the deterministic known-answer mini corpus for S0/S1 validation.

This builds ``tests/data/mini_corpus.fasta`` so the front of the local narrowing
pipeline (S0 dereplicate -> S1 tokenize) can be exercised end-to-end without any
real download or cloud round-trip. Re-running it reproduces the byte-identical
FASTA (control sequences are embedded verbatim; the decoy is drawn from an RNG
seeded with the single global ``random_seed`` from config/proteus.yaml).

CORPUS COMPOSITION (10 input records)
-------------------------------------
Controls (6) — real chain-A sequences from the locked control PDBs
(controls/references.csv), so the 3Di S1 emits is biologically meaningful:
    IsPETase   6EQE   positive PET hydrolase
    LCC_WT     4EB0   reference PET hydrolase (leaf-branch compost cutinase)
    CalB       1TCA   negative — lid-bearing lipase (alpha/beta-hydrolase)
    AChE       1EA5   negative — Torpedo acetylcholinesterase (deep gorge)
    CRL        1CRL   negative — Candida rugosa lipase 1 (lid-bearing)
    Est2       1EVQ   negative — Alicyclobacillus carboxylesterase

Near-duplicate variants (2) — IsPETase with a handful of point mutations.
At >=98% identity they sit far above the s0_dereplicate.min_seq_id (0.95)
threshold and MUST collapse into the same cluster as IsPETase:
    IsPETase_var1   IsPETase + 3 point substitutions
    IsPETase_var2   IsPETase + 4 point substitutions

Decoys (2) — clearly non-hydrolase sequences that share homology with nothing
in the set. They pass through S0/S1 untouched as their own representatives
(they are only dropped later, at the S2 fold-class triage — NOT here):
    decoy_allalpha   B-domain of Staphylococcal protein A (3-helix bundle, all-alpha)
    decoy_random     a random-but-valid amino-acid sequence (seeded)

EXPECTED S0 OUTCOME (the known answer the test asserts)
-------------------------------------------------------
    inputs:          10 sequences
    one cluster:     {IsPETase, IsPETase_var1, IsPETase_var2} -> 1 representative
    singletons:      LCC_WT, CalB, AChE, CRL, Est2, decoy_allalpha, decoy_random
    representatives: 8   (1 collapsed cluster + 7 singletons)

i.e. EXACTLY the two near-duplicate variants collapse; every distinct sequence
survives. No homology gate is applied — the negatives and decoys are retained.
"""
from __future__ import annotations

import os
import random
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "mini_corpus.fasta")

# Single global seed (mirrors config/proteus.yaml `random_seed`). Kept inline so
# this generator stays dependency-free (no PyYAML needed to rebuild the fixture).
RANDOM_SEED = 1729

# --------------------------------------------------------------------------- #
# Real control sequences (chain A, from RCSB; accessions per controls/references.csv)
# --------------------------------------------------------------------------- #
ISPETASE = (  # 6EQE — includes the expression construct's C-terminal His tag
    "MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYP"
    "TNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSS"
    "SPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNS"
    "SALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDF"
    "RTANCSLEHHHHHH"
)
LCC_WT = (  # 4EB0
    "SNPYQRGPNPTRSALTADGPFSVATYTVSRLSVSGFGGGVIYYPTGTSLTFGGIAMSPGYTADASSLAWLG"
    "RRLASHGFVVLVINTNSRFDYPDSRASQLSAALNYLRTSSPSAVRARLDANRLAVAGHSMGGGGTLRIAEQ"
    "NPSLKAAVPLTPWHTDKTFNTSVPVLIVGAEADTVAPVSQHAIPFYQNLPSTTPKVYVELDNASHFAPNSN"
    "NAAISVYTISWMKLWVDNDTRYRQFLCNVNDPALSDFRTNNRHCQ"
)
CALB = (  # 1TCA — Candida antarctica lipase B
    "LPSGSDPAFSQPKSVLDAGLTCQGASPSSVSKPILLVPGTGTTGPQSFDSNWIPLSTQLGYTPCWISPPPF"
    "MLNDTQVNTEYMVNAITALYAGSGNNKLPVLTWSQGGLVAQWGLTFFPSIRSKVDRLMAFAPDYKGTVLAG"
    "PLDALAVSAPSVWQQTTGSALTTALRNAGGLTQIVPTTNLYSATDEIVQPQVSNSPLDSSYLFNGKNVQAQ"
    "AVCGPLFVIDHAGSLTSQFSYVVGRSALRSTTGQARSADYGITDCNPLPANDLTPEQKVAAAALLAPAAAA"
    "IVAGPKQNCEPDLMPYARPFAVGKRTCSGIVTP"
)
ACHE = (  # 1EA5 — Torpedo californica acetylcholinesterase
    "DDHSELLVNTKSGKVMGTRVPVLSSHISAFLGIPFAEPPVGNMRFRRPEPKKPWSGVWNASTYPNNCQQYV"
    "DEQFPGFSGSEMWNPNREMSEDCLYLNIWVPSPRPKSTTVMVWIYGGGFYSGSSTLDVYNGKYLAYTEEVV"
    "LVSLSYRVGAFGFLALHGSQEAPGNVGLLDQRMALQWVHDNIQFFGGDPKTVTIFGESAGGASVGMHILSP"
    "GSRDLFRRAILQSGSPNCPWASVSVAEGRRRAVELGRNLNCNLNSDEELIHCLREKKPQELIDVEWNVLPF"
    "DSIFRFSFVPVIDGEFFPTSLESMLNSGNFKKTQILLGVNKDEGSFFLLYGAPGFSKDSESKISREDFMSG"
    "VKLSVPHANDLGLDAVTLQYTDWMDDNNGIKNRDGLDDIVGDHNVICPLMHFVNKYTKFGNGTYLYFFNHR"
    "ASNLVWPEWMGVIHGYEIEFVFGLPLVKELNYTAEEEALSRRIMHYWATFAKTGNPNEPHSQESKWPLFTT"
    "KEQKFIDLNTEPMKVHQRLRVQMCVFWNQFLPKLLNATAC"
)
CRL = (  # 1CRL — Candida rugosa lipase 1
    "APTATLANGDTITGLNAIINEAFLGIPFAEPPVGNLRFKDPVPYSGSLDGQKFTSYGPSCMQQNPEGTYEE"
    "NLPKAALDLVMQSKVFEAVSPSSEDCLTINVVRPPGTKAGANLPVMLWIFGGGFEVGGTSTFPPAQMITKS"
    "IAMGKPIIHVSVNYRVSSWGFLAGDEIKAEGSANAGLKDQRLGMQWVADNIAAFGGDPTKVTIFGESAGSM"
    "SVMCHILWNDGDNTYKGKPLFRAGIMQSGAMVPSDAVDGIYGNEIFDLLASNAGCGSASDKLACLRGVSSD"
    "TLEDATNNTPGFLAYSSLRLSYLPRPDGVNITDDMYALVREGKYANIPVIIGDQNDEGTFFGTSSLNVTTD"
    "AQAREYFKQSFVHASDAEIDTLMTAYPGDITQGSPFDTGILNALTPQFKRISAVLGDLGFTLARRYFLNHY"
    "TGGTKYSFLSKQLSGLPVLGTFHSNDIVFQDYLLGSGSLIYNNAFIAFATDLDPNTAGLLVKWPEYTSSSQ"
    "SGNNLMMINALGLYTGKDNFRTAGYDALFSNPPSFFV"
)
EST2 = (  # 1EVQ — Alicyclobacillus acidocaldarius esterase-2
    "MPLDPVIQQVLDQLNRMPAPDYKHLSAQQFRSQQSLFPPVKKEPVAEVREFDMDLPGRTLKVRMYRPEGVE"
    "PPYPALVYYHGGGWVVGDLETHDPVCRVLAKDGRAVVFSVDYRLAPEHKFPAAVEDAYDALQWIAERAADF"
    "HLDPARIAVGGDSAGGNLAAVTSILAKERGGPALAFQLLIYPSTGYDPAHPPASIEENAEGYLLTGGMMLW"
    "FRDQYLNSLEELTHPWFSPVLYPDLSGLPPAYIATAQYDPLRDVGKLYAEALNKAGVKVEIENFEDLIHGF"
    "AQFYSLSPGATKALVRIAEKLRDALA"
)

# B-domain of Staphylococcal protein A — a 3-helix bundle, all-alpha, no
# alpha/beta-hydrolase fold and no catalytic Ser. A clean fold-class decoy.
DECOY_ALLALPHA = "VDNKFNKEQQNAFYEILHLPNLNEEQRNGFIQSLKDDPSQSANLLAEAKKLNDAQAPK"

# Standard 20 amino acids for the synthetic decoy.
_AA20 = "ACDEFGHIKLMNPQRSTVWY"


def _mutate(seq: str, subs: list[tuple[int, str]]) -> str:
    """Apply 0-based (index, new_residue) point substitutions; assert each is a
    real change so a 'variant' can never accidentally equal the parent."""
    chars = list(seq)
    for idx, new in subs:
        assert 0 <= idx < len(chars), f"mutation index {idx} out of range"
        assert chars[idx] != new, f"substitution at {idx} is a no-op ({new})"
        chars[idx] = new
    return "".join(chars)


def _random_sequence(length: int, rng: random.Random) -> str:
    return "".join(rng.choice(_AA20) for _ in range(length))


def build_records() -> list[tuple[str, str]]:
    rng = random.Random(RANDOM_SEED)

    # IsPETase near-duplicate variants: a few conservative point substitutions in
    # the mature region (well clear of the C-terminal His tag). >=98% identity ->
    # collapse with the parent at min_seq_id 0.95.
    var1 = _mutate(ISPETASE, [(40, "S"), (95, "A"), (150, "S")])
    var2 = _mutate(ISPETASE, [(40, "S"), (60, "G"), (120, "V"), (200, "I")])

    # Deterministic non-hydrolase decoy (length 120, seeded).
    decoy_random = _random_sequence(120, rng)

    return [
        ("IsPETase", ISPETASE),
        ("LCC_WT", LCC_WT),
        ("CalB", CALB),
        ("AChE", ACHE),
        ("CRL", CRL),
        ("Est2", EST2),
        ("IsPETase_var1", var1),
        ("IsPETase_var2", var2),
        ("decoy_allalpha", DECOY_ALLALPHA),
        ("decoy_random", decoy_random),
    ]


def write_fasta(records: list[tuple[str, str]], path: str = OUT) -> None:
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for line in textwrap.wrap(seq, 60):
                fh.write(line + "\n")


if __name__ == "__main__":
    recs = build_records()
    write_fasta(recs)
    print(f"wrote {len(recs)} records -> {os.path.relpath(OUT)}")
    print("expected S0: 10 inputs -> 8 representatives "
          "(IsPETase + 2 variants collapse; 7 singletons survive)")
