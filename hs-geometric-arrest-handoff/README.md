# hs-geometric-arrest — handoff bundle

Complete handoff for building the geometric-ground-state hard-sphere project
with Claude Code. Everything in here was produced and verified in the research
sessions of August 2026; nothing is placeholder.

Start: open this directory in Claude Code. It reads CLAUDE.md automatically;
then say "Execute T0 from CLAUDE_CODE_TASKS.md".

Contents:
- CLAUDE.md ................ agent constitution (read order, rules, supervision)
- CLAUDE_CODE_TASKS.md ..... ordered tasks T0-T12, gates, landmines, file status
- IMPLEMENTATION_SPEC.md ... normative spec (modules, gates, campaign)
- REFERENCE_VALUES.json .... golden numbers, immutable, three-way verified
- docs/
    paper.tex .............. manuscript (revtex4-2; article-class preview incl.)
    paper_preview.pdf ...... compiled preview
    EXPLAINER.md ........... plain-language account
    CHATGPT_REVIEW.md ...... verified verdict on the external review + retraction
    ALTERNATIVE_ROUTES.md .. positive control, sensitivity calibration, 2D-first plan
    RESULTS_pilot.md ....... pilot results incl. negative findings
    archive/ ............... superseded first spec (history only)
- prototypes/ .............. research-phase code; correct output, port cleanly
- data/pilot3d, pilot2d .... actual pilot configurations for validating ports

Known open items before submission (also in the docs):
- paper.tex affiliation is a placeholder; bibliography must be completed
  (add Donev-Torquato-Stillinger cond-mat/0408550 for the K=4 prior art).
- The 2D replication (T6) is the gate on all 3D physics claims.
