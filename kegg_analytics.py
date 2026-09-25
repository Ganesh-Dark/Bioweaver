import sys
sys.dont_write_bytecode = True
import os
import requests

# Base directory of this script (used to locate local Kegg_files/)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
KEGG_FILES = os.path.join(BASE_DIR, "Kegg_files")


# ─────────────────────────────────────────────
# INTERNAL HELPERS
# ─────────────────────────────────────────────

def _link_api(target, source):
    """
    Calls https://rest.kegg.jp/link/<target>/<source> and returns a list of
    (source_id_bare, target_id_bare) tuples.
    """
    # target: KEGG database to link TO  (e.g. 'rn', 'ec', 'cpd', 'drug')
    # source: KEGG ID or database to link FROM (e.g. 'map00010', 'hsa')
    url = f"https://rest.kegg.jp/link/{target}/{source}"
    print(f"[Analytics] KEGG Link API: /link/{target}/{source}")
    res = requests.get(url, timeout=15)
    if res.status_code != 200 or not res.text.strip():
        return []
    pairs = []
    for line in res.text.strip().split("\n"):
        if "\t" in line:
            a, b = line.split("\t")
            pairs.append((a.split(":")[-1], b.split(":")[-1]))
    return pairs


def _load_tsv_lookup(filename):
    """
    Loads a two-column TSV file into a dict: {id: name_string}.
    """
    # filename: Name of the file inside Kegg_files/ (e.g. 'Reactions.tsv')
    path = os.path.join(KEGG_FILES, filename)
    lookup = {}
    if not os.path.exists(path):
        print(f"[Analytics] WARNING: {filename} not found.")
        return lookup
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t", 1)
            if len(parts) == 2:
                lookup[parts[0]] = parts[1]
    return lookup


def _parse_equation(equation_str):
    """
    Splits a KEGG equation string into (substrates_list, products_list,
    is_reversible) using <=> or =>.
    Stoichiometric numbers and leading/trailing whitespace are preserved.
    """
    # equation_str: Raw equation string from Reactions.tsv
    if "<=>" in equation_str:
        left, right = equation_str.split("<=>", 1)
        reversible = True
    elif "=>" in equation_str:
        left, right = equation_str.split("=>", 1)
        reversible = False
    else:
        return [], [], None

    subs = [s.strip() for s in left.split("+") if s.strip()]
    prods = [p.strip() for p in right.split("+") if p.strip()]
    return subs, prods, reversible


def _get_pathway_reactions(pathway_id):
    """
    Returns a list of dicts for every reaction in the given pathway. Each dict has:
      id, name, equation, substrates, products, reversible, ec_numbers.
    Both the reaction list and EC mappings are built from the KEGG Link API +
    local Reactions.tsv, so zero LLM tokens are used.
    """
    # pathway_id: KEGG map ID such as map00010

    # 1. Get all reaction IDs for this pathway
    rn_pairs = _link_api("rn", pathway_id)
    rn_ids = set(p[1] for p in rn_pairs)
    if not rn_ids:
        print(f"[Analytics] No reactions found for {pathway_id}")
        return []

    # 2. Bulk-fetch EC assignments for the specific reactions
    rn_to_ecs = {}
    try:
        print("[Analytics] Fetching EC mappings for reactions...")
        ec_rn_res = requests.get("https://rest.kegg.jp/link/ec/rn", timeout=15)
        if ec_rn_res.status_code == 200:
            for line in ec_rn_res.text.strip().split("\n"):
                if "\t" in line:
                    rn_part, ec_part = line.split("\t")
                    rn_id = rn_part.replace("rn:", "")
                    ec_id = ec_part.replace("ec:", "")
                    if rn_id in rn_ids:
                        rn_to_ecs.setdefault(rn_id, []).append(ec_id)
    except Exception as e:
        print(f"[Analytics] Failed to fetch EC mappings: {e}")

    # 3. Cross-reference with local Reactions.tsv
    rxn_lookup = _load_tsv_lookup("Reactions.tsv")
    results = []
    for rn_id in rn_ids:
        raw = rxn_lookup.get(rn_id, "")
        equation = ""
        name = rn_id

        # Split name from equation
        if "<=>" in raw or "=>" in raw:
            chunks = raw.split("; ")
            for i in range(len(chunks) - 1, -1, -1):
                if "<=>" in chunks[i] or "=>" in chunks[i]:
                    equation = chunks[i]
                    name = "; ".join(chunks[:i]) or rn_id
                    break
        else:
            name = raw

        subs, prods, rev = _parse_equation(equation)

        results.append({
            "id": rn_id,
            "name": name.split("; ")[0] if name else rn_id,
            "equation": equation,
            "substrates": subs,
            "products": prods,
            "reversible": rev,
            "ec_numbers": rn_to_ecs.get(rn_id, [])
        })

    print(f"[Analytics] Loaded {len(results)} reactions for {pathway_id}")
    return results


def _get_pathway_compounds(pathway_id):
    """
    Returns a set of bare compound IDs (e.g. 'C00022') for the given pathway
    via the KEGG Link API.
    """
    # pathway_id: KEGG map ID such as map00010
    pairs = _link_api("cpd", pathway_id)
    return set(p[1] for p in pairs)


def _get_pathway_ecs(pathway_id):
    """
    Returns a dict {ec_id: ec_name} for the given pathway via KEGG Link API
    + local enzyme.tsv.
    """
    # pathway_id: KEGG map ID such as map00010
    pairs = _link_api("ec", pathway_id)
    ec_ids = set(p[1] for p in pairs)
    ec_lookup = _load_tsv_lookup("enzyme.tsv")
    result = {}
    for ec_id in ec_ids:
        raw_name = ec_lookup.get(ec_id, "")
        result[ec_id] = raw_name.split(";")[0].strip() if raw_name else ec_id
    return result


def _compound_name(cpd_id, lookup=None):
    """
    Returns the primary name for a compound ID from local TSV. Falls back to
    the ID itself if not found.
    """
    # cpd_id: KEGG compound ID like C00022
    # lookup: Optional pre-loaded dict to avoid re-reading the file
    if lookup is None:
        lookup = _load_tsv_lookup("kegg_compounds_synonyms.tsv")
    raw = lookup.get(cpd_id, "")
    return raw.split(";")[0].strip() if raw else cpd_id


# ─────────────────────────────────────────────
# PUBLIC ANALYTICS FUNCTIONS
# ─────────────────────────────────────────────

def find_isozymes(pathway_id):
    """
    Returns only those reactions in the pathway that are catalyzed by MORE than
    one EC number (isozymes / redundant enzymes).
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] find_isozymes: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    isozyme_rxns = [r for r in reactions if len(r["ec_numbers"]) > 1]
    print(f"[Analytics] Found {len(isozyme_rxns)} reactions with multiple ECs.")
    return {
        "pathway_id": pathway_id,
        "analysis": "isozymes",
        "total_reactions": len(reactions),
        "results": isozyme_rxns
    }


def find_irreversible_reactions(pathway_id):
    """
    Returns only the irreversible reactions (=> only, not <=>) in the pathway.
    These are often rate-limiting or regulatory steps.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] find_irreversible_reactions: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    irreversible = [r for r in reactions if r.get("reversible") is False]
    print(f"[Analytics] Found {len(irreversible)} irreversible reactions.")
    return {
        "pathway_id": pathway_id,
        "analysis": "irreversible_reactions",
        "total_reactions": len(reactions),
        "results": irreversible
    }


def find_reversible_reactions(pathway_id):
    """
    Returns only the reversible reactions (<=>) in the pathway.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] find_reversible_reactions: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    reversible = [r for r in reactions if r.get("reversible") is True]
    print(f"[Analytics] Found {len(reversible)} reversible reactions.")
    return {
        "pathway_id": pathway_id,
        "analysis": "reversible_reactions",
        "total_reactions": len(reactions),
        "results": reversible
    }


def find_orphan_reactions(pathway_id):
    """
    Returns reactions that have no EC number assigned. These are 'orphan'
    reactions – chemistry that is known to occur but whose enzyme is not yet
    characterised.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] find_orphan_reactions: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    orphans = [r for r in reactions if not r["ec_numbers"]]
    print(f"[Analytics] Found {len(orphans)} orphan reactions.")
    return {
        "pathway_id": pathway_id,
        "analysis": "orphan_reactions",
        "total_reactions": len(reactions),
        "results": orphans
    }


def get_pathway_coverage(pathway_id):
    """
    Returns a summary showing what percentage of reactions in the pathway have
    at least one EC number assigned.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] get_pathway_coverage: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    total = len(reactions)
    with_ec = len([r for r in reactions if r["ec_numbers"]])
    without_ec = total - with_ec
    coverage_pct = round((with_ec / total) * 100, 1) if total > 0 else 0
    print(f"[Analytics] Coverage: {with_ec}/{total} ({coverage_pct}%)")
    return {
        "pathway_id": pathway_id,
        "analysis": "pathway_coverage",
        "total_reactions": total,
        "reactions_with_ec": with_ec,
        "orphan_reactions": without_ec,
        "coverage_percent": coverage_pct
    }


def get_ec_class_breakdown(pathway_id):
    """
    Groups all EC numbers in the pathway by their top-level EC class:
      Class 1 = Oxidoreductases
      Class 2 = Transferases
      Class 3 = Hydrolases
      Class 4 = Lyases
      Class 5 = Isomerases
      Class 6 = Ligases
      Class 7 = Translocases
    This reveals the dominant reaction type in a pathway at a glance.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    EC_CLASSES = {
        "1": "Oxidoreductases",
        "2": "Transferases",
        "3": "Hydrolases",
        "4": "Lyases",
        "5": "Isomerases",
        "6": "Ligases",
        "7": "Translocases"
    }
    print(f"\n[Analytics] get_ec_class_breakdown: {pathway_id}")
    ecs = _get_pathway_ecs(pathway_id)
    breakdown = {}
    for ec_id, ec_name in ecs.items():
        # EC numbers start with the class digit, e.g. '2.7.1.1'
        cls = ec_id.split(".")[0]
        label = EC_CLASSES.get(cls, f"Class {cls}")
        breakdown.setdefault(label, []).append({"ec": ec_id, "name": ec_name})

    # Sort classes numerically
    sorted_breakdown = dict(sorted(breakdown.items()))
    print(f"[Analytics] Found {len(ecs)} ECs across {len(sorted_breakdown)} classes.")
    return {
        "pathway_id": pathway_id,
        "analysis": "ec_class_breakdown",
        "total_ec_numbers": len(ecs),
        "breakdown_by_class": sorted_breakdown
    }


def find_hub_metabolites(pathway_id, top_n=10):
    """
    Identifies the most frequently occurring compounds (hub metabolites) across
    all reactions in the pathway. A compound appearing in many reactions acts
    as a central "metabolic hub" (e.g. ATP, NAD+, CoA).
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    # top_n: How many top hub compounds to return
    print(f"\n[Analytics] find_hub_metabolites: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    cpd_lookup = _load_tsv_lookup("kegg_compounds_synonyms.tsv")

    # Remove stoichiometric numbers (e.g. "2 ATP" -> "ATP")
    freq = {}
    for r in reactions:
        all_cpds = r["substrates"] + r["products"]
        for cpd_str in all_cpds:
            # Drop leading stoichiometric number
            tokens = cpd_str.strip().split()
            cpd_name = " ".join(tokens[1:]) if tokens and tokens[0].replace("n", "").replace("+", "").isdigit() else cpd_str.strip()
            freq[cpd_name] = freq.get(cpd_name, 0) + 1

    # Sort by frequency, take top N
    sorted_cpds = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:top_n]
    results = [{"compound": name, "appearances": count} for name, count in sorted_cpds]
    print(f"[Analytics] Top {top_n} hub metabolites found.")
    return {
        "pathway_id": pathway_id,
        "analysis": "hub_metabolites",
        "total_reactions_scanned": len(reactions),
        "results": results
    }


def find_pathway_compound_intersection(pathway_ids):
    """
    Finds all compounds that are shared between 2 or more KEGG pathways.
    """
    if not pathway_ids or len(pathway_ids) < 2:
        return {"error": "Provide at least 2 pathway IDs for intersection."}
        
    print(f"\n[Analytics] find_pathway_compound_intersection: {pathway_ids}")
    
    # Initialize intersection with the first pathway
    shared = _get_pathway_compounds(pathway_ids[0])
    
    # Intersect with all other pathways
    for pid in pathway_ids[1:]:
        shared &= _get_pathway_compounds(pid)

    cpd_lookup = _load_tsv_lookup("kegg_compounds_synonyms.tsv")
    results = [{"id": c, "name": _compound_name(c, cpd_lookup)} for c in sorted(shared)]
    print(f"[Analytics] Shared compounds: {len(results)}")
    return {
        "pathways": pathway_ids,
        "analysis": "compound_intersection",
        "shared_count": len(results),
        "results": results
    }


def find_pathway_ec_intersection(pathway_ids):
    """
    Finds all enzymes (EC numbers) shared between 2 or more KEGG pathways.
    """
    if not pathway_ids or len(pathway_ids) < 2:
        return {"error": "Provide at least 2 pathway IDs for intersection."}
        
    print(f"\n[Analytics] find_pathway_ec_intersection: {pathway_ids}")
    
    # We need the full dict for the final names, so keep the first one
    base_ecs = _get_pathway_ecs(pathway_ids[0])
    shared_ids = set(base_ecs.keys())
    
    for pid in pathway_ids[1:]:
        ecs = _get_pathway_ecs(pid)
        shared_ids &= set(ecs.keys())
        
    results = [{"ec": ec, "name": base_ecs[ec]} for ec in sorted(shared_ids)]
    print(f"[Analytics] Shared ECs: {len(results)}")
    return {
        "pathways": pathway_ids,
        "analysis": "ec_intersection",
        "shared_count": len(results),
        "results": results
    }


def find_pathway_unique_compounds(pathway_id_1, pathway_id_2):
    """
    Finds compounds UNIQUE to pathway_id_1 (not present in pathway_id_2).
    Useful for identifying what makes one pathway distinct from another.
    """
    # pathway_id_1: The pathway whose unique compounds we want
    # pathway_id_2: The reference pathway to compare against
    print(f"\n[Analytics] find_pathway_unique_compounds: unique to {pathway_id_1}")
    cpds1 = _get_pathway_compounds(pathway_id_1)
    cpds2 = _get_pathway_compounds(pathway_id_2)
    unique = cpds1 - cpds2
    cpd_lookup = _load_tsv_lookup("kegg_compounds_synonyms.tsv")
    results = [{"id": c, "name": _compound_name(c, cpd_lookup)} for c in sorted(unique)]
    print(f"[Analytics] Unique compounds in {pathway_id_1}: {len(results)}")
    return {
        "pathway_id": pathway_id_1,
        "compared_against": pathway_id_2,
        "analysis": "unique_compounds",
        "unique_count": len(results),
        "results": results
    }


def trace_compound_in_pathway(pathway_id, compound_name_query):
    """
    Finds all reactions in a pathway where a given compound appears as a
    substrate OR a product. Shows the metabolic 'role' of that compound in
    the pathway.
    """
    # pathway_id: KEGG map ID
    # compound_name_query: Partial or full name to match (case-insensitive)
    print(f"\n[Analytics] trace_compound_in_pathway: '{compound_name_query}' in {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    query_lower = compound_name_query.lower()

    as_substrate = []
    as_product = []

    for r in reactions:
        for s in r["substrates"]:
            if query_lower in s.lower():
                as_substrate.append({"reaction_id": r["id"], "equation": r["equation"], "ec_numbers": r["ec_numbers"]})
                break
        for p in r["products"]:
            if query_lower in p.lower():
                as_product.append({"reaction_id": r["id"], "equation": r["equation"], "ec_numbers": r["ec_numbers"]})
                break

    print(f"[Analytics] '{compound_name_query}' is substrate in {len(as_substrate)} rxns, product in {len(as_product)} rxns.")
    return {
        "pathway_id": pathway_id,
        "compound_query": compound_name_query,
        "analysis": "compound_trace",
        "consumed_in": as_substrate,
        "produced_in": as_product
    }


def find_drug_targets_in_pathway(pathway_id):
    """
    Cross-references the enzymes (EC numbers) in the given pathway with the
    KEGG Drug database to find which enzymes are known drug targets.
    Uses the KEGG Link API: /link/drug/pathway.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010')
    print(f"\n[Analytics] find_drug_targets_in_pathway: {pathway_id}")

    # Get drugs linked to this pathway
    drug_pairs = _link_api("drug", pathway_id)
    if not drug_pairs:
        return {
            "pathway_id": pathway_id,
            "analysis": "drug_targets",
            "results": [],
            "note": "No drugs found linked to this pathway in KEGG."
        }

    drug_ids = [p[1] for p in drug_pairs]
    drug_lookup = _load_tsv_lookup("kegg_diseases.tsv")  # fallback name lookup

    # Fetch drug names via KEGG list API for up to 10 drugs (API limit)
    results = []
    for drug_id in drug_ids[:30]:
        drug_name = drug_id  # fallback
        res = requests.get(f"https://rest.kegg.jp/list/{drug_id}", timeout=8)
        if res.status_code == 200 and res.text.strip():
            parts = res.text.strip().split("\t", 1)
            if len(parts) == 2:
                drug_name = parts[1].split(";")[0].strip()
        results.append({"drug_id": drug_id, "name": drug_name})

    print(f"[Analytics] Found {len(results)} drug targets for {pathway_id}.")
    return {
        "pathway_id": pathway_id,
        "analysis": "drug_targets",
        "total_drugs": len(drug_ids),
        "results": results
    }


def summarize_pathway(pathway_id):
    """
    Returns a full statistical summary of a pathway in one call:
    - Total reactions, reversible vs irreversible
    - Orphan reactions (no EC)
    - Coverage percentage (reactions with EC)
    - EC class breakdown
    - Top 5 hub metabolites
    This gives a complete snapshot of the pathway without the LLM processing anything.
    """
    # pathway_id: KEGG map ID (e.g. 'map00010' for Glycolysis)
    print(f"\n[Analytics] summarize_pathway: {pathway_id}")
    reactions = _get_pathway_reactions(pathway_id)
    total = len(reactions)
    if total == 0:
        return {"error": f"No reactions found for pathway {pathway_id}"}

    reversible = [r for r in reactions if r.get("reversible") is True]
    irreversible = [r for r in reactions if r.get("reversible") is False]
    orphans = [r for r in reactions if not r["ec_numbers"]]
    with_ec = total - len(orphans)
    isozyme_rxns = [r for r in reactions if len(r["ec_numbers"]) > 1]
    coverage = round((with_ec / total) * 100, 1)

    # EC class count
    EC_CLASSES = {
        "1": "Oxidoreductases", "2": "Transferases", "3": "Hydrolases",
        "4": "Lyases", "5": "Isomerases", "6": "Ligases", "7": "Translocases"
    }
    ec_classes = {}
    for r in reactions:
        for ec in r["ec_numbers"]:
            cls = ec.split(".")[0]
            label = EC_CLASSES.get(cls, f"Class {cls}")
            ec_classes[label] = ec_classes.get(label, 0) + 1

    # Top 5 hub metabolites
    freq = {}
    for r in reactions:
        for cpd_str in r["substrates"] + r["products"]:
            tokens = cpd_str.strip().split()
            name = " ".join(tokens[1:]) if tokens and tokens[0].replace("n", "").replace("+", "").isdigit() else cpd_str.strip()
            freq[name] = freq.get(name, 0) + 1
    top5 = sorted(freq.items(), key=lambda x: x[1], reverse=True)[:5]

    return {
        "pathway_id": pathway_id,
        "analysis": "full_summary",
        "total_reactions": total,
        "reversible": len(reversible),
        "irreversible": len(irreversible),
        "reactions_with_ec": with_ec,
        "orphan_reactions": len(orphans),
        "coverage_percent": coverage,
        "reactions_with_multiple_ecs": len(isozyme_rxns),
        "ec_class_counts": ec_classes,
        "top_5_hub_metabolites": [{"compound": n, "appearances": c} for n, c in top5]
    }


# ─────────────────────────────────────────────
# Detailed Script Workings:
# 1. _link_api() is the single KEGG network call. All public functions call
#    this helper to get cross-reference IDs from the KEGG Link REST API.
# 2. _load_tsv_lookup() reads local Kegg_files TSV files into a fast dict
#    so that all subsequent lookups are purely in-memory with no network call.
# 3. _get_pathway_reactions() is the core aggregator. It combines _link_api()
#    results for reactions AND ec numbers, then cross-references with local
#    Reactions.tsv in a single O(n) pass.
# 4. All public functions (find_isozymes, find_irreversible_reactions, etc.)
#    filter the pre-built reaction list using standard Python list comprehensions.
#    No LLM tokens are consumed anywhere in this file.
# 5. trace_compound_in_pathway() does a simple case-insensitive substring match
#    on the equation substrate/product strings to locate the compound.
# 6. summarize_pathway() aggregates all statistics in a single pass over the
#    reaction list, making it the most token-efficient way to describe a pathway.
