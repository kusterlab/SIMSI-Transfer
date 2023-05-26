# MSFragger input handling for SIMSI

### msms.txt vs psm.tsv files
<br>

General remarks msms.txt
- msms.txt contains the results of all runs

<br><br>
General remarks psm.tsv
- https://fragpipe.nesvilab.org/docs/tutorial_fragpipe_outputs.html#psmtsv
- psm.tsv is generated for each different 'experiment' in FragPipe run
- "spectrum" column identifier, consisting of...
  - raw file name w/o extension
  - MS2 spectrum number
  - MS2 spectrum number
  - Charge

Columns
- msms.txt <br>
  'Raw file', 'Scan number', 'Scan index', 'Sequence', 'Length',
  'Missed cleavages', 'Modifications', 'Modified sequence',
  'Oxidation (M) Probabilities', 'Oxidation (M) Score diffs',
  'Acetyl (Protein N-term)', 'Oxidation (M)', 'Proteins', 'Gene Names',
  'Protein Names', 'Charge', 'Fragmentation', 'Mass analyzer', 'Type',
  'Scan event number', 'Isotope index', 'm/z', 'Mass', 'Mass error [ppm]',
  'Mass error [Da]', 'Simple mass error [ppm]', 'Retention time', 'PEP',
  'Score', 'Delta score', 'Score diff', 'Localization prob',
  'Combinatorics', 'PIF', 'Fraction of total spectrum',
  'Base peak fraction', 'Precursor full scan number',
  'Precursor Intensity', 'Precursor apex fraction',
  'Precursor apex offset', 'Precursor apex offset time', 'Matches',
  'Intensities', 'Mass deviations [Da]', 'Mass deviations [ppm]',
  'Masses', 'Number of matches', 'Intensity coverage', 'Peak coverage',
  'Neutral loss level', 'ETD identification type', 'Reverse',
  'All scores', 'All sequences', 'All modified sequences',
  'Reporter intensity corrected 1', 'Reporter intensity corrected 2',
  'Reporter intensity corrected 3', 'Reporter intensity corrected 4',
  'Reporter intensity corrected 5', 'Reporter intensity corrected 6',
  'Reporter intensity corrected 7', 'Reporter intensity corrected 8',
  'Reporter intensity corrected 9', 'Reporter intensity corrected 10',
  'Reporter intensity corrected 11', 'Reporter intensity 1',
  'Reporter intensity 2', 'Reporter intensity 3', 'Reporter intensity 4',
  'Reporter intensity 5', 'Reporter intensity 6', 'Reporter intensity 7',
  'Reporter intensity 8', 'Reporter intensity 9', 'Reporter intensity 10',
  'Reporter intensity 11', 'Reporter mass deviation [mDa] 1',
  'Reporter mass deviation [mDa] 2', 'Reporter mass deviation [mDa] 3',
  'Reporter mass deviation [mDa] 4', 'Reporter mass deviation [mDa] 5',
  'Reporter mass deviation [mDa] 6', 'Reporter mass deviation [mDa] 7',
  'Reporter mass deviation [mDa] 8', 'Reporter mass deviation [mDa] 9',
  'Reporter mass deviation [mDa] 10', 'Reporter mass deviation [mDa] 11',
  'Reporter PIF', 'Reporter fraction', 'id', 'Protein group IDs',
  'Peptide ID', 'Mod. peptide ID', 'Evidence ID',
  'Oxidation (M) site IDs'
  <br><br>
- psm.tsv <br>
  'Spectrum', 'Spectrum File', 'Peptide', 'Modified Peptide', 'Prev AA',
  'Next AA', 'Peptide Length', 'Charge', 'Retention', 'Observed Mass',
  'Calibrated Observed Mass', 'Observed M/Z', 'Calibrated Observed M/Z',
  'Calculated Peptide Mass', 'Calculated M/Z', 'Delta Mass',
  'Expectation', 'Hyperscore', 'Nextscore', 'PeptideProphet Probability',
  'Number of Enzymatic Termini', 'Number of Missed Cleavages',
  'Protein Start', 'Protein End', 'Intensity', 'Assigned Modifications',
  'Observed Modifications', 'Purity', 'Is Unique', 'Protein',
  'Protein ID', 'Entry Name', 'Gene', 'Protein Description',
  'Mapped Genes', 'Mapped Proteins', 'Quan Usage', '30pm', '300pm', '1nm',
  '3nm', '10nm', '30nm', '100nm', '300nm', '1um', '10um', 'Ref'
