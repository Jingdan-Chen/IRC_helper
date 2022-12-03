# IRC_helper

---

This program is for processing IRC out file produced by Gaussian 16, fitting IRC curve, and generating EDA/single point computation file for Q-CHEM.

## Files and Folders

### working.ipynb

---

Notebook file as the panel of interactive program

### funarg.py

---

functions and arguments to be imported to working.ipynb

### origin_file

---

Here contains g16 IRC out files to be read

**Recommended IRC generation key words**:

*#p irc=(calcfc,maxpoints=30,lqa,stepsize=8) {ps: the "p" following "#" doesn't matter}* 

*--Link1--*

*#p irc=(restart,maxpoints=60,lqa,stepsize=15)*

## more_scripts

---

Help you to use generated .gjf/.xyz to generate .inp files for ALMO-EDA(Q-CHEM)

## irc_analysis_out

---

Output .csv and .xyz files of "working.ipynb"