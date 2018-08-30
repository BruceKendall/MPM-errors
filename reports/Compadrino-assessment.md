# Assessment of Compadrino MPM error analysis

- It appears that the compadrinos were using the development version of the database rather than 2.01: the analysis of 10.1046/j.1523-1739.1998.97054.x used SpeciesAuthor "Acinonyx_jubatus_2" and says "Matrix in paper hast just 8 classes, compadrino used 36 (maybe from additional source?)"; but in my version SpeciesAuthor is "Acinonyx_jubatus" and the model has 8 classes.

# Soren (from top)
## 10.1046/j.1523-1739.1998.97054.x
- The doi doesn't work for some reason, but I found the paper.
- It's a postbreeding Leslie matrix. Soren correctly classified it as post
- The F's have parent survival, but Soren coded SurvInRep as "unknown"
- Text and table 1 says age at first reproduction is 24, but model sets it at 30, so ReproWithMat should by "no;" Soren gets this correct.

## 10.1111/j.1365-2664.2010.01846.x
- The paper has no matrix, and there is no evidence that they constructed one (no analysis of MPM). Unclear where the Comadre matrix came from
- However, Soren provided an analysis.

## 10.1007/BF01237655
- prebreeding Leslie matrix, with offspring included in F; Soren gets this right
- Soren incorrectly says "no" for ReproWithMat

## 10.2307/1941572
- Soren correctly says "no matrix in paper" but does an analysis anyway.

## 10.2193/0022-541X(2006)70[1094:GASOBB]2.0.CO;2
- Soren correctly says "no matrix in paper" but does an analysis anyway.

## 10.2307/1941949
- Soren correctly says prebreeding but incorrectly says parent survival
- S incorrectly says "no" for reprowithmat
- S says "NA" for growth transition but it should be "observed" (based on mean of one-year transitions)


# Jakob (from top)
## 10.1007/s10646-010-0507-y
- Matrix not shown but must exist; presumably obtained from author. 
- Clearest to look at D. magna "spinosad" treatment as it has lowest survival
- postbreeding leslie model; j gets this right
- F = m; survinfec is none. J incorrectly says "parent"
- J correctly says no to rwm and none to ls

## 10.2307/2937107
- postbreeding, with no sif. J says Na and unknown
- j correctly says rwm is no
- mat rule seems to be "observed" J says "other"


