"""Run main.main() for one scoring variant and dump the final-score rows (no Excel/Word)."""
import os, pickle, sys
from pathlib import Path
sys.path.insert(0, ".")
label, cfg = sys.argv[1], sys.argv[2]
weight = float(open(Path(cfg).parent / f"weight_{label}.txt").read())
import analysis.interaction_scoring as isc
isc.V2_COMPONENT_WEIGHTS["coexpression_gse64349"] = weight
import output.excel as excel_mod, output.word_report as word_mod
from output.report_v2 import build_workbook_sheets
out = Path(cfg).parent / f"rows_{label}.pkl"
def fake_excel(config, blast_classification, output_path, interaction_result=None, word_report_filename=None, provenance=None):
    sheets = build_workbook_sheets(config, blast_classification, interaction_result)
    pickle.dump(sheets["final_score_rows"], open(out, "wb"))
    Path(output_path).write_text("dumped"); return Path(output_path)
excel_mod.write_classification_workbook = fake_excel
word_mod.write_word_report = lambda *a, **k: Path("x")
import main
main.main(["--config", cfg])
