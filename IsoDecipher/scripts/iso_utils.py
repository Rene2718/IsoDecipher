import gffutils
import pandas as pd

def explain_utr_status(tx_input, db):
    if isinstance(db,str):
        db = gffutils.FeatureDB(db)

    try:
        tx = db[tx_input]
    except Exception:
        matches = [
            t for t in db.features_of_type("transcript")
            if t.attributes.get("transcript_name", [""])[0] == tx_input
        ]
        if not matches:
            print(f"❌ No transcript found for '{tx_input}' (checked both ID and transcript_name)")
            return None
        if len(matches) > 1:
            print(f"⚠️ Multiple transcripts matched '{tx_input}'; using the first one ({matches[0].id})")
        tx = matches[0]

    gene = tx.attributes.get('gene_name',[""])[0]
    tx_id = tx.id
    tx_name = tx.attributes.get("transcript_name",[""])[0]
    strand = tx.strand
    tags = tx.attributes.get("tag",[])
    exons = list(db.children(tx, featuretype="exon"))
    cds = list(db.children(tx, featuretype="CDS"))
    print(f"\n🧬 Transcript: {tx_id} ({tx_name}) | Gene: {gene} | Strand: {strand} ")
    if tags: 
        print(f"Tags:{', '.join(tags)}")
    if not exons:
        print("❌ No exons found — malformed transcript.")
        return {"utr_len": None, "status": "no_exons"}
    # Determine last exon depending on strand
    if strand == "+":
        last_exon = max(exons, key=lambda e: e.end)
        strand_info = "+"
    elif strand == "-":
        last_exon = min(exons, key=lambda e: e.start)
        strand_info = "-"
    else:
        last_exon = max(exons, key=lambda e: e.end)
        strand_info = "unknown"

        print("⚠️  Missing or unknown strand; using genomic order as fallback.")
    print(f"Last exon: {last_exon.start:,d}-{last_exon.end:,d}")
    
    if not cds:
        print("❌ No CDS entries → utr_status = 'missing_CDS'")
        return {"utr_len": None, "status": "missing_CDS"}
    cds_starts = [c.start for c in cds]
    cds_ends =[c.end for c in cds]
    cds_start = min(cds_starts)
    cds_end = max(cds_ends)
    print(f"CDS range: {cds_start:,d}-{cds_end:,d}")
    
    # Compute UTR length depending on strand

    if strand=="+":
        utr_len =last_exon.end-cds_end
    elif strand == "-":
        utr_len = cds_start-last_exon.start
    else:
        utr_len = None

    # Status
    if utr_len is None:
        status = "missing_strand"
        color= "\x1b[31m"
    elif utr_len <= 0:
        status = "invalid_or_negative"
        color = "\x1b[33m"  # yellow
    else:
        status = "valid"
        color = "\x1b[32m"  # green
    
    print(f"Compute UTR length:{utr_len if utr_len is not None else 'N/A'} bp ")
    print(f"{color}UTR status: {status}\x1b[0m")
    
    if status == "invalid_or_negative":
        print("⚠️  CDS overlaps or extends beyond the last exon — no distinct 3′ UTR region.")
    elif status == "missing_strand":
        print("⚠️  Transcript has no strand info — cannot determine UTR direction.")
    elif status == "valid":
        print("✅  Valid UTR detected.")
    return {
        "gene": gene,
        "tx_id":tx_id,
        "tx_name": tx_name,
        "strand": strand or "unknown",
        "utr_len": utr_len,
        "status": status,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "last_exon_start": last_exon.start,
        "last_exon_end": last_exon.end,
    }
