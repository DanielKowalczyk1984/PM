#!/usr/bin/env python3
import re
import pydot
import glob

list_files = glob.glob("wt*_4_cuts_*.txt")
list_files.sort()

for files in list_files:
    cuts = open(files, 'r')
    pattern_cuts_list = r"(.*)_cuts_(\d+)\.txt"
    match_list = re.findall(pattern_cuts_list,files)
    instance = match_list[0][0]
    zero_half_iteration = match_list[0][1]
    j = 0

    for it in cuts.readlines():
        file_path = instance + "_graph.dot"
        pattern_cuts = r"([-+])?\sx(\d+)"
        m_cuts = re.findall(pattern_cuts, it)
        i = 0
        file_handler = open(file_path, 'r')
        name_file = instance + "_cuts_%s_represent_%03d" %(zero_half_iteration,j)
        file_new = open(name_file + ".dot", 'w')
        for line in file_handler:
            color = "red" if m_cuts[i][0] == '-' else "green"
            a = int(m_cuts[i][1])
            b = str(a)
            pattern = r'(\d+) ->\s(\d+)\[label\s=\s+\"(%s)\"(\s?.+)\]' % (b)
            m = re.search(pattern, line)
            if m:
                file_new.write("%s -> %s [label=\"%s\" %s color=%s];\n" %
                           (m.group(1), m.group(2), m.group(3), m.group(4), color))
                i = i + 1
                if i == len(m_cuts):
                    i = 0
            else:
                file_new.write(line)

        j = j + 1

        file_handler.close()
        file_new.close()
        (graph,) = pydot.graph_from_dot_file(name_file + ".dot")
        graph.write_pdf(name_file + ".pdf")
