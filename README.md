# IRC_helper

---

This program is for processing IRC out file produced by Gaussian 16, fitting IRC curve, and generating EDA/single point computation file for Q-CHEM.

Wei-Feng Zheng, **Jingdan Chen**, Xiaotian Qi*, Zhongxing Huang*. Asymmetric decarboxylative protonation enabled by an anchoring group that enhances noncovalent interactions. ()

## Files and Folders

### working.ipynb

---

Notebook file as the panel of interactive program

### funarg.py

---

functions and arguments to be imported to working.ipynb

### origin_file

---

Here contains g16 IRC out file(s) to be read
Hint: 1) To get a complete IRC containing ts information, please do not read Hessian from check file;
2) The code for multi-IRCfiles is not robust, please do not try it to ensure the outcome is correct


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
