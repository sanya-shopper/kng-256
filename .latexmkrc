# After every bibtex run, regroup the produced .bbl into thematic clusters.
$bibtex = 'bibtex %O %B && python3 tools/clusterbbl.py %B.bbl';
