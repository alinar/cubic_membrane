#! /usr/bin/env python
import membrane.membrane as mem
m=mem.Membrane(2,220)
m.form = "G"
m.Make()
m.PrintInfo()
m.WriteOnFile("g.pdb")
